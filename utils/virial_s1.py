#!/usr/bin/env python3
"""Run ground-state and ROKS-DFTB calculations for a Becke virial S1 energy."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


HARTREE_TO_EV = 27.2113845
ENERGY_TAGS = ("extrapolated0_energy", "mermin_energy", "total_energy")
STATE_SETTINGS_MARKER = "@STATE_SETTINGS@"
TRIPLET_ONLY_BEGIN = "@TRIPLET_ONLY_BEGIN@"
TRIPLET_ONLY_END = "@TRIPLET_ONLY_END@"
TAG_HEADER = re.compile(
    r"^\s*(?P<name>[^:\s]+)\s*:"
    r"(?P<type>real|integer|logical|complex):"
    r"(?P<rank>\d+):(?P<shape>[0-9,]*)\s*$"
)


class DriverError(RuntimeError):
    """An expected input, execution, or results error."""


def parse_hhubbard(value: str) -> tuple[str, float]:
    """Parse an ELEMENT=VALUE Hartree-Hubbard command-line argument."""
    try:
        element, raw_value = value.split("=", maxsplit=1)
        element = element.strip()
        hhubbard = float(raw_value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"invalid HHubbard '{value}'; expected ELEMENT=VALUE"
        ) from error
    if not element or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", element):
        raise argparse.ArgumentTypeError(
            f"invalid species name in HHubbard '{value}'"
        )
    if hhubbard <= 0.0:
        raise argparse.ArgumentTypeError(
            f"HHubbard must be positive in '{value}'"
        )
    return element, hhubbard


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a closed-shell ground state and a ROKS triplet, then "
            "calculate S1 = E_triplet - E_ground + K_if."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Successful execution
--------------------
1. Build a DFTB+ executable containing the ROKS virial implementation. The
   executable must write the scalar tag 'roks_virial_energy' for the ROKS run.

2. Prepare all files referenced by dftb_in.hsd, for example geo.gen and any
   local Slater-Koster files, so they are accessible from both calculation
   directories. Absolute Slater-Koster paths are usually the simplest choice.

3. Either prepare separate inputs manually:

       ground/{'dftb_in.hsd'}
       triplet/{'dftb_in.hsd'}

   or start from the supplied example 'utils/virial_template.hsd'. A template
   must contain exactly one line

       {STATE_SETTINGS_MARKER}

   inside 'Hamiltonian = DFTB {{ ... }}'. The driver replaces that line with
   the closed-shell settings for ground and with ROKS, spin-polarisation, and
   RoksVirial.HHubbard settings for triplet. Input needed only by the triplet,
   such as SpinConstants, must be enclosed by marker lines:

       {TRIPLET_ONLY_BEGIN}
       SpinConstants = {{ ... }}
       {TRIPLET_ONLY_END}

   The enclosed text is removed from the ground input and retained in the
   triplet input.

4. Ensure the common part of the input enables results.tag:

       Options {{
           WriteResultsTag = Yes
       }}

   For a generated triplet input, supply one Hartree-only Hubbard parameter
   for every chemical species with repeated --hhubbard arguments.

5. Run both calculations and combine their results:

       virial_s1.py --template template.hsd \\
           --hhubbard H=0.68353 --hhubbard C=0.49748 \\
           --executable /path/to/dftb+

   The default calculation directories are 'ground' and 'triplet'. Existing
   dftb_in.hsd files are protected; pass --overwrite-inputs to replace them.

6. Inspect ground/dftb.log and triplet/dftb.log if DFTB+ fails. On success,
   the driver reads the common total-energy tag from both results.tag files
   and 'roks_virial_energy' from the triplet result. It evaluates

       E(S1) = E(ROKS) - E(ground) + K_if

   and prints the triplet excitation, K_if, and S1 excitation in eV. The full
   machine-readable result is written to s1_results.json by default.

Reuse completed calculations
----------------------------
To calculate the final result without running DFTB+ again:

    virial_s1.py --ground-dir ground --triplet-dir triplet --reuse-results

Manual triplet inputs
---------------------
By default, a manual triplet input must contain both 'RoksVirial' and
'HHubbard'. This guards against accidentally using the ordinary SCC gamma
interaction. Use --allow-scc-gamma only when that older interaction is
deliberately required for a comparison.
""",
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
        help="ground-state calculation directory (default: ground)",
    )
    parser.add_argument(
        "--triplet-dir",
        type=Path,
        default=Path("triplet"),
        help="ROKS triplet calculation directory (default: triplet)",
    )
    parser.add_argument(
        "--input-name",
        default="dftb_in.hsd",
        help="input filename in both directories (default: dftb_in.hsd)",
    )
    parser.add_argument(
        "--template",
        type=Path,
        help=(
            "generate both inputs from a template containing exactly one "
            f"{STATE_SETTINGS_MARKER} line"
        ),
    )
    parser.add_argument(
        "--hhubbard",
        metavar="ELEMENT=VALUE",
        action="append",
        type=parse_hhubbard,
        default=[],
        help=(
            "Hartree-only Hubbard value in atomic units; repeat for every "
            "species, for example --hhubbard H=0.68353"
        ),
    )
    parser.add_argument(
        "--allow-scc-gamma",
        action="store_true",
        help=(
            "allow an ordinary SCC-gamma virial integral; by default the "
            "triplet input must contain RoksVirial.HHubbard"
        ),
    )
    parser.add_argument(
        "--unpaired-electrons",
        type=int,
        default=2,
        help="unpaired electrons in the generated ROKS input (default: 2)",
    )
    parser.add_argument(
        "--overwrite-inputs",
        action="store_true",
        help="allow template mode to replace existing input files",
    )
    parser.add_argument(
        "--results-name",
        default="results.tag",
        help="tagged results filename (default: results.tag)",
    )
    parser.add_argument(
        "--log-name",
        default="dftb.log",
        help="captured process-output filename (default: dftb.log)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("s1_results.json"),
        help="combined JSON output path (default: s1_results.json)",
    )
    parser.add_argument(
        "--reuse-results",
        action="store_true",
        help="do not run DFTB+; combine existing results.tag files",
    )
    return parser.parse_args()


def unique_hhubbards(
    items: list[tuple[str, float]],
) -> dict[str, float]:
    values: dict[str, float] = {}
    for element, value in items:
        key = element[0].upper() + element[1:].lower()
        if key in values:
            raise DriverError(f"Duplicate --hhubbard value for species {key}")
        values[key] = value
    return values


def format_hhubbard_block(values: dict[str, float]) -> str:
    if not values:
        raise DriverError(
            "Template mode requires at least one --hhubbard ELEMENT=VALUE "
            "unless --allow-scc-gamma is used"
        )
    lines = ["RoksVirial = {", "    HHubbard = {"]
    for element, value in values.items():
        lines.append(f"        {element} = {value:.12g}")
    lines.extend(("    }", "}"))
    return "\n".join(lines)


def generate_inputs(
    template_path: Path,
    ground_dir: Path,
    triplet_dir: Path,
    input_name: str,
    unpaired_electrons: int,
    hhubbards: dict[str, float],
    allow_scc_gamma: bool,
    overwrite: bool,
) -> None:
    template_path = template_path.expanduser().resolve()
    if not template_path.is_file():
        raise DriverError(f"Template file does not exist: {template_path}")
    if unpaired_electrons < 1:
        raise DriverError("--unpaired-electrons must be positive")

    template = template_path.read_text(encoding="utf-8")
    marker_pattern = re.compile(
        rf"^(?P<indent>[ \t]*){re.escape(STATE_SETTINGS_MARKER)}[ \t]*$",
        re.MULTILINE,
    )
    matches = list(marker_pattern.finditer(template))
    if len(matches) != 1:
        raise DriverError(
            f"Template must contain exactly one line containing only "
            f"{STATE_SETTINGS_MARKER}; found {len(matches)}"
        )

    triplet_block_pattern = re.compile(
        rf"^[ \t]*{re.escape(TRIPLET_ONLY_BEGIN)}[ \t]*\n"
        rf"(?P<body>.*?)"
        rf"^[ \t]*{re.escape(TRIPLET_ONLY_END)}[ \t]*(?:\n|$)",
        re.MULTILINE | re.DOTALL,
    )
    triplet_blocks = list(triplet_block_pattern.finditer(template))
    begin_count = template.count(TRIPLET_ONLY_BEGIN)
    end_count = template.count(TRIPLET_ONLY_END)
    if begin_count != end_count or begin_count != len(triplet_blocks):
        raise DriverError(
            "Every @TRIPLET_ONLY_BEGIN@ line must have one matching "
            "@TRIPLET_ONLY_END@ line"
        )

    ground_settings = "Roks = No"
    triplet_lines = [
        "Roks = Yes",
        "RoksMaxIterations = 50",
        "RoksTolerance = 1.0e-8",
    ]
    if hhubbards:
        triplet_lines.append(format_hhubbard_block(hhubbards))
    elif not allow_scc_gamma:
        format_hhubbard_block(hhubbards)  # raises the detailed error
    triplet_lines.extend(
        (
            "SpinPolarisation = Colinear {",
            f"    UnpairedElectrons = {unpaired_electrons}",
            "}",
        )
    )
    triplet_settings = "\n".join(triplet_lines)

    jobs = []
    for directory, settings, keep_triplet_blocks in (
        (ground_dir, ground_settings, False),
        (triplet_dir, triplet_settings, True),
    ):
        directory = directory.expanduser().resolve()
        input_path = directory / input_name
        if input_path.exists() and not overwrite:
            raise DriverError(
                f"Refusing to replace {input_path}; use --overwrite-inputs"
            )
        jobs.append((directory, input_path, settings, keep_triplet_blocks))

    indent = matches[0].group("indent")
    for directory, input_path, settings, keep_triplet_blocks in jobs:
        directory.mkdir(parents=True, exist_ok=True)
        indented_settings = settings.replace("\n", "\n" + indent)
        rendered = marker_pattern.sub(
            indent + indented_settings, template, count=1
        )
        rendered = triplet_block_pattern.sub(
            lambda match: match.group("body") if keep_triplet_blocks else "",
            rendered,
        )
        input_path.write_text(rendered, encoding="utf-8")
        print(f"Generated {input_path}")


def resolve_executable(command: str) -> str:
    candidate = Path(command).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        candidate = candidate.resolve()
        if not candidate.is_file():
            raise DriverError(f"Executable does not exist: {candidate}")
        if not candidate.stat().st_mode & 0o111:
            raise DriverError(f"File is not executable: {candidate}")
        return str(candidate)
    resolved = shutil.which(command)
    if resolved is None:
        raise DriverError(f"Executable not found on PATH: {command}")
    return resolved


def validate_calculation_directory(
    directory: Path, input_name: str
) -> Path:
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise DriverError(f"Calculation directory does not exist: {directory}")
    input_path = directory / input_name
    if not input_path.is_file():
        raise DriverError(f"Required input does not exist: {input_path}")
    return directory


def verify_triplet_kernel(
    input_path: Path, allow_scc_gamma: bool
) -> str:
    text = re.sub(
        r"#.*$",
        "",
        input_path.read_text(encoding="utf-8"),
        flags=re.MULTILINE,
    )
    pure_hartree = (
        re.search(r"\bRoksVirial\b", text, re.IGNORECASE) is not None
        and re.search(r"\bHHubbard\b", text, re.IGNORECASE) is not None
    )
    if pure_hartree:
        return "pure_hartree_custom_hhubbard"
    if not allow_scc_gamma:
        raise DriverError(
            f"{input_path} does not contain RoksVirial.HHubbard; "
            "add Hartree-only parameters or use --allow-scc-gamma"
        )
    return "ordinary_scc_gamma"


def run_calculation(
    executable: str,
    directory: Path,
    log_name: str,
    results_name: str,
) -> None:
    log_path = directory / log_name
    results_path = directory / results_name
    old_mtime = results_path.stat().st_mtime_ns if results_path.exists() else None

    print(f"Running {directory.name}: {executable}")
    with log_path.open("w", encoding="utf-8") as log_file:
        completed = subprocess.run(
            [executable],
            cwd=directory,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise DriverError(
            f"DFTB+ failed in {directory} with exit code "
            f"{completed.returncode}; see {log_path}"
        )
    if not results_path.is_file():
        raise DriverError(
            f"DFTB+ did not create {results_path}; enable WriteResultsTag"
        )
    if old_mtime is not None and results_path.stat().st_mtime_ns == old_mtime:
        raise DriverError(f"Results file was not updated: {results_path}")


def parse_tagged_file(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    entries: dict[str, Any] = {}
    index = 0
    while index < len(lines):
        match = TAG_HEADER.match(lines[index])
        if match is None:
            index += 1
            continue
        name = match.group("name")
        value_type = match.group("type")
        rank = int(match.group("rank"))
        shape = tuple(
            int(item)
            for item in match.group("shape").split(",")
            if item
        )
        count = 1
        for extent in shape:
            count *= extent
        if rank == 0:
            count = 1

        tokens: list[str] = []
        index += 1
        while index < len(lines) and len(tokens) < count:
            if TAG_HEADER.match(lines[index]):
                break
            tokens.extend(lines[index].split())
            index += 1
        if len(tokens) < count:
            raise DriverError(f"Incomplete tag '{name}' in {path}")

        raw = tokens[:count]
        if value_type == "real":
            values: list[Any] = [
                float(value.replace("D", "E").replace("d", "e"))
                for value in raw
            ]
        elif value_type == "integer":
            values = [int(value) for value in raw]
        elif value_type == "logical":
            values = [
                value.lower() in ("t", "true", ".true.") for value in raw
            ]
        else:
            continue
        entries[name] = values[0] if rank == 0 else values
    return entries


def require_scalar(
    entries: dict[str, Any], name: str, path: Path
) -> float:
    if name not in entries:
        raise DriverError(f"Required tag '{name}' is missing from {path}")
    value = entries[name]
    if isinstance(value, list):
        raise DriverError(f"Tag '{name}' is not scalar in {path}")
    return float(value)


def combine_results(
    ground_path: Path,
    triplet_path: Path,
    kernel: str,
    hhubbards: dict[str, float],
) -> dict[str, Any]:
    ground = parse_tagged_file(ground_path)
    triplet = parse_tagged_file(triplet_path)
    energy_tag = next(
        (
            name
            for name in ENERGY_TAGS
            if name in ground and name in triplet
        ),
        None,
    )
    if energy_tag is None:
        raise DriverError(
            "No common total-energy tag; tried " + ", ".join(ENERGY_TAGS)
        )

    ground_h = require_scalar(ground, energy_tag, ground_path)
    triplet_h = require_scalar(triplet, energy_tag, triplet_path)
    virial_h = require_scalar(
        triplet, "roks_virial_energy", triplet_path
    )
    if virial_h < -1.0e-12:
        raise DriverError(
            f"ROKS virial integral is negative: {virial_h:.12e} H"
        )

    triplet_exc_h = triplet_h - ground_h
    s1_h = triplet_exc_h + virial_h
    result: dict[str, Any] = {
        "method": "Becke virial exciton model with ROKS-DFTB",
        "virial_kernel": kernel,
        "energy_tag": energy_tag,
        "ground_energy_h": ground_h,
        "triplet_energy_h": triplet_h,
        "triplet_excitation_h": triplet_exc_h,
        "triplet_excitation_ev": triplet_exc_h * HARTREE_TO_EV,
        "virial_kif_h": virial_h,
        "virial_kif_ev": virial_h * HARTREE_TO_EV,
        "s1_excitation_h": s1_h,
        "s1_excitation_ev": s1_h * HARTREE_TO_EV,
        "hartree_to_ev": HARTREE_TO_EV,
        "ground_results": str(ground_path),
        "triplet_results": str(triplet_path),
    }
    if hhubbards:
        result["hartree_hubbard"] = hhubbards
    if "roks_transq_sum" in triplet:
        transition_charge_sum = require_scalar(
            triplet, "roks_transq_sum", triplet_path
        )
        if abs(transition_charge_sum) > 1.0e-8:
            raise DriverError(
                "ROKS transition charges are not neutral: "
                f"sum(q_A) = {transition_charge_sum:.12e}"
            )
        result["transition_charge_sum"] = transition_charge_sum
    if "roks_open_orbitals" in triplet:
        result["open_shell_orbitals"] = triplet["roks_open_orbitals"]
    return result


def print_report(result: dict[str, Any]) -> None:
    print("\nBecke virial S1 calculation")
    print("===========================")
    print(f"Virial kernel:          {result['virial_kernel']}")
    print(f"Energy tag:             {result['energy_tag']}")
    print(f"Ground energy:          {result['ground_energy_h']: .12f} H")
    print(f"Triplet energy:         {result['triplet_energy_h']: .12f} H")
    print(f"Triplet excitation:     {result['triplet_excitation_ev']: .8f} eV")
    print(f"Virial K_if:            {result['virial_kif_ev']: .8f} eV")
    print(f"S1 excitation:          {result['s1_excitation_ev']: .8f} eV")
    if "transition_charge_sum" in result:
        print(
            "Transition-charge sum: "
            f"{result['transition_charge_sum']: .6e}"
        )
    if "open_shell_orbitals" in result:
        orbitals = result["open_shell_orbitals"]
        print(f"Open-shell orbitals:    {orbitals[0]} {orbitals[1]}")


def main() -> int:
    args = parse_args()
    hhubbards = unique_hhubbards(args.hhubbard)

    if args.template is not None:
        generate_inputs(
            args.template,
            args.ground_dir,
            args.triplet_dir,
            args.input_name,
            args.unpaired_electrons,
            hhubbards,
            args.allow_scc_gamma,
            args.overwrite_inputs,
        )

    ground_dir = validate_calculation_directory(
        args.ground_dir, args.input_name
    )
    triplet_dir = validate_calculation_directory(
        args.triplet_dir, args.input_name
    )
    if ground_dir == triplet_dir:
        raise DriverError("Ground and triplet directories must differ")

    kernel = verify_triplet_kernel(
        triplet_dir / args.input_name, args.allow_scc_gamma
    )

    if not args.reuse_results:
        executable = resolve_executable(args.executable)
        run_calculation(
            executable, ground_dir, args.log_name, args.results_name
        )
        run_calculation(
            executable, triplet_dir, args.log_name, args.results_name
        )

    ground_results = ground_dir / args.results_name
    triplet_results = triplet_dir / args.results_name
    if not ground_results.is_file():
        raise DriverError(f"Results file does not exist: {ground_results}")
    if not triplet_results.is_file():
        raise DriverError(f"Results file does not exist: {triplet_results}")

    result = combine_results(
        ground_results, triplet_results, kernel, hhubbards
    )
    print_report(result)

    output_path = args.output.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
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
