"""Molecular dynamics on a REAL structure of a REAL scanned gene.

This replaces the generic Trp-cage benchmark. It reads the structure manifest
produced by `python -m synatvis.fetch_structures_batch`, picks the structure for
the gene being scanned, and runs OpenMM on that.

Structure selection, in order:
  1. $MD_STRUCTURE  -- an explicit .pdb path, or a UniProt accession in the manifest
  2. the protein arriving on stdin is matched to the manifest by length
  3. otherwise the first manifest entry, clearly reported as a fallback

Whatever it uses is echoed back in the JSON as `structure_source`, so a result can
never silently be about a different protein than you think.

Setup:
    $env:MDSIM_CMD = '"<venv python>" "<path to this file>"'
    python -m synatvis scan <your.fasta> --plugins
"""
import gc
import json
import os
import sys
import tempfile

STRUCT_DIR = os.environ.get(
    "MD_STRUCTURE_DIR",
    r"C:\Users\Computer Arena\Downloads\synatvis_structures")
N_STEPS = 5000  # 10 ps at 2 fs -- a laptop-sized demonstration run, not a production run

raw = sys.stdin.read()
query_prot = "".join(l.strip() for l in raw.splitlines() if not l.startswith(">"))


def load_manifest():
    p = os.path.join(STRUCT_DIR, "structures_manifest.json")
    if not os.path.isfile(p):
        return []
    with open(p, encoding="utf-8") as fh:
        return json.load(fh).get("structures", [])


def choose(manifest):
    want = os.environ.get("MD_STRUCTURE", "").strip()
    if want:
        if os.path.isfile(want):
            return want, f"explicit path {os.path.basename(want)}"
        for m in manifest:
            if m["uniprot"].lower() == want.lower() or m["refseq"].lower() == want.lower():
                return os.path.join(STRUCT_DIR, m["file"]), \
                       f"{m['uniprot']} ({m['protein']}) selected by MD_STRUCTURE"
    if query_prot and manifest:
        best = min(manifest, key=lambda m: abs((m.get("length_aa") or 0) - len(query_prot)))
        gap = abs((best.get("length_aa") or 0) - len(query_prot))
        if gap <= max(15, 0.10 * len(query_prot)):
            return os.path.join(STRUCT_DIR, best["file"]), (
                f"{best['uniprot']} ({best['protein']}) matched to the scanned protein "
                f"by length {best['length_aa']} aa vs {len(query_prot)} aa")
    if manifest:
        m = manifest[0]
        return os.path.join(STRUCT_DIR, m["file"]), (
            f"FALLBACK -- no confident match to the scanned protein; used "
            f"{m['uniprot']} ({m['protein']}). Treat as a pipeline demonstration only.")
    return None, "no structures available"


manifest = load_manifest()
pdb_path, provenance = choose(manifest)
if not pdb_path or not os.path.isfile(pdb_path):
    print(json.dumps({"error": "no structure available",
                      "hint": "run: python -m synatvis.fetch_structures_batch --limit 100",
                      "structure_dir": STRUCT_DIR}))
    sys.exit(1)

from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation, NoCutoff, HBonds, DCDReporter
from openmm import LangevinMiddleIntegrator
from openmm.unit import kelvin, picosecond, picoseconds
import mdtraj as md

tmp = tempfile.mkdtemp()
fixed = os.path.join(tmp, "fixed.pdb")
traj_f = os.path.join(tmp, "traj.dcd")


def run():
    """Kept in a function so the DCD handle is released before it is read back --
    on Windows, reading it in the same scope raises 'premature end of file'."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater=False)   # crystallographic waters break the force field
    fixer.addMissingHydrogens(7.0)
    with open(fixed, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    pdb = PDBFile(fixed)
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    sim = Simulation(pdb.topology, system, integ)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy()
    sim.reporters.append(DCDReporter(traj_f, 10))
    sim.step(N_STEPS)


run()
gc.collect()

# mdtraj's DCD reader writes C-level chatter straight to fd 1, which would corrupt
# the single-JSON-line contract; redirect the file descriptor for just this call.
saved, devnull = os.dup(1), os.open(os.devnull, os.O_WRONLY)
os.dup2(devnull, 1)
try:
    traj = md.load(traj_f, top=fixed)
finally:
    os.dup2(saved, 1)
    os.close(devnull)
    os.close(saved)

traj.superpose(traj, 0)
print(json.dumps({
    "rmsd_nm": round(float(md.rmsd(traj, traj, 0)[-1]), 4),
    "radius_of_gyration_nm": round(float(md.compute_rg(traj)[-1]), 4),
    "sim_time_ns": round(N_STEPS * 0.002 / 1000, 4),
    "force_field": "amber14-all + implicit/gbn2 (OpenMM)",
    "temperature_k": 300,
    "structure_source": provenance,
}))
