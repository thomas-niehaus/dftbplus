#!/usr/bin/env python3
"""Run the two DFTB+ calculations required for a virial S1 energy."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


HARTREE_TO_EV = 27.2113845
MARKER = "@STATE_SETTINGS@"
INPUT_NAME = "dftb_in.hsd"
RESULTS_NAME = "results.tag"
LOG_NAME = "dftb.log"
ENERGY_TAGS = ("extrapolated0_energy", "mermin_energy", "total_energy")

HELP = f"""
The template must contain exactly one line with

    {MARKER}

inside Hamiltonian = DFTB. SpinPolarisation and SpinConstants are common to
both calculations; the driver sets UnpairedElectrons to 0 for the ground state
and to 2 for ROKS. A minimal H/C template has this form:

    Geometry = GenFormat {{
        <<< geo.gen
    }}

    Hamiltonian = DFTB {{
        SCC = Yes
        SCCTolerance = 1.0e-8

        MaxAngularMomentum = {{
            H = "s"
            C = "p"
        }}

        SlaterKosterFiles = Type2FileNames {{
            Prefix = "/absolute/path/to/mio-1-1/"
            Separator = "-"
            Suffix = ".skf"
        }}

        {MARKER}

        SpinConstants = {{
            C = {{
                -0.04559 -0.02930
                -0.02930 -0.02755
            }}
            H = {{-0.07925}}
            ShellResolvedSpin = Yes
        }}
    }}

    Options = {{
        WriteResultsTag = Yes
    }}

    ParserOptions = {{
        ParserVersion = 14
    }}

The generated inputs are always overwritten. Both calculations are always run.
Files referenced by the template must be accessible from both calculation
directories. For example, use an absolute Slater-Koster prefix and place a
copy of geo.gen in both directories.

Example:

    virial_s1.py --template template.hsd \\
        --hhubbard H=0.68353 --hhubbard C=0.49748 \\
        --executable /path/to/dftb+
"""


class DriverError(RuntimeError):
    """Expected input, execution, or result error."""


def parse_hhubbard(text: str) -> tuple[str, float]:
    try:
        element, raw_value = text.split("=", 1)
        element = element.strip()
        value = float(raw_value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"invalid Hubbard parameter '{text}'; expected ELEMENT=VALUE"
        ) from error
    if not element:
        raise argparse.ArgumentTypeError(f"missing element in '{text}'")
    if value <= 0.0:
        raise argparse.ArgumentTypeError(
            f"Hartree Hubbard parameter must be positive in '{text}'"
        )
    return element, value


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a closed-shell ground state and a two-open-orbital ROKS "
            "calculation, then evaluate E(S1) = E(ROKS) - E(S0) + K_if."
        ),
        epilog=HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--template",
        type=Path,
        required=True,
        help=f"input template containing one {MARKER} line",
    )
    parser.add_argument(
        "-x",
        "--executable",
        default="dftb+",
        help="DFTB+ executable path or command name (default: dftb+)",
    )
    parser.add_argument(
        "--ground-dir",
        type=Path,
        default=Path("ground"),
        help="ground-state directory (default: ground)",
    )
    parser.add_argument(
        "--triplet-dir",
        type=Path,
        default=Path("triplet"),
        help="ROKS directory (default: triplet)",
    )
    parser.add_argument(
        "--hhubbard",
        metavar="ELEMENT=VALUE",
        action="append",
        type=parse_hhubbard,
        required=True,
        help="Hartree-only Hubbard parameter; repeat for every species",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("s1_results.json"),
        help="combined JSON result (default: s1_results.json)",
    )
    return parser.parse_args()


def collect_hhubbards(items: list[tuple[str, float]]) -> dict[str, float]:
    values: dict[str, float] = {}
    for element, value in items:
        element = element[0].upper() + element[1:].lower()
        if element in values:
            raise DriverError(f"duplicate Hubbard parameter for {element}")
        values[element] = value
    return values


def format_virial_block(hhubbards: dict[str, float]) -> str:
    lines = ["RoksVirial = {", "    HHubbard = {"]
    for element, value in hhubbards.items():
        lines.append(f"        {element} = {value:.12g}")
    lines.extend(("    }", "}"))
    return "\n".join(lines)


def state_settings(roks: bool, hhubbards: dict[str, float]) -> str:
    if not roks:
        return """\
Roks = No
SpinPolarisation = Colinear {
    UnpairedElectrons = 0
}"""

    return f"""\
Roks = Yes
RoksMaxIterations = 50
RoksTolerance = 1.0e-8
{format_virial_block(hhubbards)}
SpinPolarisation = Colinear {{
    UnpairedElectrons = 2
}}"""


def render_input(template: str, settings: str) -> str:
    marker_lines = [
        line for line in template.splitlines() if line.strip() == MARKER
    ]
    if len(marker_lines) != 1:
        raise DriverError(
            f"template must contain exactly one {MARKER} line; "
            f"found {len(marker_lines)}"
        )
    marker_line = marker_lines[0]
    indent = marker_line[: len(marker_line) - len(marker_line.lstrip())]
    settings = settings.replace("\n", "\n" + indent)
    return template.replace(marker_line, indent + settings, 1)


def prepare_directory(
    directory: Path, template: str, settings: str
) -> Path:
    directory = directory.expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / INPUT_NAME).write_text(
        render_input(template, settings), encoding="utf-8"
    )
    print(f"Wrote {directory / INPUT_NAME}")
    return directory


def run_dftb(executable: str, directory: Path) -> None:
    results_path = directory / RESULTS_NAME
    results_path.unlink(missing_ok=True)
    log_path = directory / LOG_NAME
    print(f"Running {directory.name}: {executable}")
    try:
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(
                [executable],
                cwd=directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
    except FileNotFoundError as error:
        raise DriverError(f"executable not found: {executable}") from error
    if result.returncode != 0:
        raise DriverError(
            f"DFTB+ failed in {directory} with exit code "
            f"{result.returncode}; see {log_path}"
        )
    if not results_path.is_file():
        raise DriverError(
            f"{results_path} was not written; enable WriteResultsTag"
        )


def read_scalar_tags(path: Path) -> dict[str, float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    tags: dict[str, float] = {}
    for index, line in enumerate(lines[:-1]):
        fields = line.split(":")
        if len(fields) >= 3 and fields[1] == "real" and fields[2] == "0":
            try:
                tags[fields[0].strip()] = float(
                    lines[index + 1].strip().replace("D", "E").replace("d", "e")
                )
            except ValueError as error:
                raise DriverError(
                    f"invalid scalar value for '{fields[0].strip()}' in {path}"
                ) from error
    return tags


def select_energy_tag(
    ground: dict[str, float], roks: dict[str, float]
) -> str:
    for name in ENERGY_TAGS:
        if name in ground and name in roks:
            return name
    raise DriverError(
        "ground and ROKS results have no common total-energy tag"
    )


def calculate_result(
    ground_path: Path, roks_path: Path, hhubbards: dict[str, float]
) -> dict[str, object]:
    ground = read_scalar_tags(ground_path)
    roks = read_scalar_tags(roks_path)
    energy_tag = select_energy_tag(ground, roks)
    try:
        virial_h = roks["roks_virial_energy"]
    except KeyError as error:
        raise DriverError(
            f"'roks_virial_energy' is missing from {roks_path}"
        ) from error

    ground_h = ground[energy_tag]
    roks_h = roks[energy_tag]
    triplet_excitation_h = roks_h - ground_h
    s1_h = triplet_excitation_h + virial_h
    return {
        "energy_tag": energy_tag,
        "ground_energy_h": ground_h,
        "roks_energy_h": roks_h,
        "triplet_excitation_ev": triplet_excitation_h * HARTREE_TO_EV,
        "virial_kif_ev": virial_h * HARTREE_TO_EV,
        "s1_excitation_ev": s1_h * HARTREE_TO_EV,
        "hartree_hubbard": hhubbards,
    }


def print_result(result: dict[str, object]) -> None:
    print("\nBecke virial S1 calculation")
    print("===========================")
    print(f"Ground energy:      {result['ground_energy_h']: .12f} H")
    print(f"ROKS energy:        {result['roks_energy_h']: .12f} H")
    print(f"Triplet excitation: {result['triplet_excitation_ev']: .8f} eV")
    print(f"Virial K_if:        {result['virial_kif_ev']: .8f} eV")
    print(f"S1 excitation:      {result['s1_excitation_ev']: .8f} eV")


def main() -> int:
    args = parse_arguments()
    template_path = args.template.expanduser().resolve()
    if not template_path.is_file():
        raise DriverError(f"template does not exist: {template_path}")

    ground_dir = args.ground_dir.expanduser().resolve()
    roks_dir = args.triplet_dir.expanduser().resolve()
    if ground_dir == roks_dir:
        raise DriverError("ground and ROKS directories must differ")

    template = template_path.read_text(encoding="utf-8")
    hhubbards = collect_hhubbards(args.hhubbard)
    ground_dir = prepare_directory(
        ground_dir, template, state_settings(False, hhubbards)
    )
    roks_dir = prepare_directory(
        roks_dir, template, state_settings(True, hhubbards)
    )

    run_dftb(args.executable, ground_dir)
    run_dftb(args.executable, roks_dir)

    result = calculate_result(
        ground_dir / RESULTS_NAME,
        roks_dir / RESULTS_NAME,
        hhubbards,
    )
    print_result(result)

    output_path = args.output.expanduser().resolve()
    output_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"\nWrote {output_path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except DriverError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
