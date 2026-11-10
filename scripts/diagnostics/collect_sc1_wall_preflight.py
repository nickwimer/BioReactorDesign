#!/usr/bin/env python3
"""Collect SC1 wall/MRF evidence without conversion, compilation or CFD.

Reads archived cases and existing scratch VTKs; writes only a new scratch
directory. This is an evidence report, never a physical validation PASS.
Requires numpy, pyvista and Foundation OF9 foamDictionary on PATH.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import numpy as np

from sc1_wall_mesh import read_file, zone_adjacency


ARCHIVE = Path('/projects/microenviro/final_runs/reacting_Julia_08282026')
VTK_ROOT = Path('/scratch/nwimer/sc1_surface_evolution_20261103')
HEIGHT_EDGES = np.arange(0.30, 1.551, 0.025)
HEIGHT_CENTERS = (HEIGHT_EDGES[:-1] + HEIGHT_EDGES[1:]) / 2
RADIAL_BANDS = ((0.0, 0.05), (0.35, 0.45), (0.45, 0.50), (0.50, 0.53))
THRESHOLDS = (0.1, 0.5, 0.9)
SECTORS = 12


def all_crossings(profile, threshold):
    """Keep upward/downward crossings; never interpolate across empty bins."""
    rows = []
    for i in range(len(profile) - 1):
        a, b = profile[i:i + 2]
        if not np.isfinite(a + b):
            continue
        direction = ('down' if a >= threshold and b < threshold else
                     'up' if a < threshold and b >= threshold else None)
        if direction:
            height = HEIGHT_CENTERS[i] + (threshold - a) / (b - a) * (
                HEIGHT_CENTERS[i + 1] - HEIGHT_CENTERS[i])
            rows.append(dict(direction=direction, height_m=float(height)))
    return rows


def surface_profiles(centers, volumes, alpha):
    """Cell-volume-weighted sector profiles, not sharp connected surfaces."""
    centers, volumes, alpha = (np.asarray(a, dtype=float) for a in (centers, volumes, alpha))
    if centers.shape != (len(alpha), 3) or volumes.shape != alpha.shape or alpha.ndim != 1:
        raise ValueError('Inconsistent cell arrays')
    if (not np.all(np.isfinite(centers)) or not np.all(np.isfinite(alpha))
            or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0)):
        raise ValueError('Invalid cell coordinates, liquid fractions or volumes; no cells discarded')
    y = centers[:, 1]
    radius = np.hypot(centers[:, 0], centers[:, 2])
    angles = np.mod(np.arctan2(centers[:, 2], centers[:, 0]), 2 * np.pi)
    sectors = np.minimum((angles * SECTORS / (2 * np.pi)).astype(int), SECTORS - 1)
    rows = []
    for lo, hi in RADIAL_BANDS:
        annulus = (radius >= lo) & (radius < hi)
        for sector in range(SECTORS):
            mask = annulus & (sectors == sector)
            inside = mask & (y >= HEIGHT_EDGES[0]) & (y < HEIGHT_EDGES[-1])
            bins = np.searchsorted(HEIGHT_EDGES, y[inside], side='right') - 1
            n = len(HEIGHT_CENTERS)
            v = np.bincount(bins, weights=volumes[inside], minlength=n)
            av = np.bincount(bins, weights=alpha[inside] * volumes[inside], minlength=n)
            count = np.bincount(bins, minlength=n)
            profile = np.full(n, np.nan)
            np.divide(av, v, out=profile, where=v > 0)
            above = mask & (y >= 1.02)
            rows.append(dict(radial_band_m=[lo, hi], sector=sector,
                             angle_degrees=[30 * sector, 30 * (sector + 1)],
                             cell_counts=count.tolist(), sampled_volume_m3=v.tolist(),
                             mean_alpha=[float(a) if np.isfinite(a) else None for a in profile],
                             occupied_height_bins=int(np.count_nonzero(v)),
                             liquid_above_1_02_L=float(1000 * np.sum(alpha[above] * volumes[above])),
                             crossings={str(t): all_crossings(profile, t) for t in THRESHOLDS}))
    return dict(cells=len(alpha), total_liquid_L=float(1000 * np.sum(alpha * volumes)),
                liquid_above_1_02_L=float(1000 * np.sum(alpha[y >= 1.02] * volumes[y >= 1.02])),
                alpha_min=float(alpha.min()), alpha_max=float(alpha.max()),
                out_of_bounds_cells=int(np.count_nonzero((alpha < 0) | (alpha > 1))),
                height_centers_m=HEIGHT_CENTERS.tolist(), sectors=rows)


def surface_summary(data):
    lines = [f"Liquid={data['total_liquid_L']:.6f} L; above 1.02 m={data['liquid_above_1_02_L']:.6f} L; "
             f"alpha range=[{data['alpha_min']:.6g}, {data['alpha_max']:.6g}], "
             f"out-of-bounds cells={data['out_of_bounds_cells']}"]
    for band in RADIAL_BANDS:
        rows = [r for r in data['sectors'] if r['radial_band_m'] == list(band)]
        downward = [[c['height_m'] for c in r['crossings']['0.5'] if c['direction'] == 'down']
                    for r in rows]
        heights = [max(h) for h in downward if h]
        envelope = f'{min(heights):.4f} to {max(heights):.4f} m' if heights else 'none'
        coverage = [r['occupied_height_bins'] for r in rows]
        lines.append(f'r={band[0]:.2f}-{band[1]:.2f}: highest downward alpha=0.5 crossings {envelope}; '
                     f'missing sectors={sum(not h for h in downward)}/12; '
                     f'multiple-downward sectors={sum(len(h) > 1 for h in downward)}/12; '
                     f'occupied height bins/sector={min(coverage)}-{max(coverage)}/{len(HEIGHT_CENTERS)}')
    return lines


def existing_vtk(root, name, time):
    folder = root / name / 'VTK'
    # Restrict to root volume meshes; never select a nested boundary export.
    candidates = sorted(folder.glob(f'*_{time}.vtk')) + sorted(folder.glob(f'*_{time}.vtu'))
    if len(candidates) != 1:
        raise ValueError(f'Expected one root volume VTK for {name}/{time}; found {candidates}')
    return candidates[0]


def evaluate_vtk(path):
    import pyvista as pv
    before = path.stat()
    mesh = pv.read(str(path))
    if 'alpha.liquid' not in mesh.cell_data:
        raise ValueError(f'Missing cell alpha.liquid in {path}')
    sizes = mesh.compute_cell_sizes(length=False, area=False, volume=True)
    data = surface_profiles(mesh.cell_centers().points, sizes.cell_data['Volume'],
                            mesh.cell_data['alpha.liquid'])
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'VTK changed during analysis: {path}')
    data['input'] = dict(path=str(path.resolve()), bytes=after.st_size, sha256=digest.hexdigest())
    return data


def prepare_output(out, inputs, scratch_roots=(Path('/scratch'), Path('/kfs3/scratch'))):
    resolved = out.resolve()
    if not any(resolved.is_relative_to(root.resolve()) for root in scratch_roots):
        raise ValueError('Output must be a new directory under /scratch or /kfs3/scratch')
    for src in inputs:
        src = src.resolve()
        if resolved == src or src in resolved.parents or resolved in src.parents:
            raise ValueError('Output and inputs must not overlap')
    if out.exists() or out.is_symlink():
        raise ValueError(f'Output already exists; preserve it: {out}')
    out.mkdir()  # parent must already exist; no existing directory is reused


def collect_configuration(case, name, out, foam_dictionary, issues):
    """Keep literal dictionaries and OF9-read boundary entries; no expansion/code execution."""
    records = []
    paths = set()
    for pattern in ('constant/*Properties*', 'constant/*Transport*', 'constant/global*',
                    'constant/g', 'constant/fvModels', 'constant/fvConstraints',
                    'constant/dynamicMeshDict', 'system/controlDict', 'system/fvSchemes',
                    'system/fvSolution', 'system/setFieldsDict*', 'system/topoSetDict*'):
        paths.update(p for p in case.glob(pattern) if p.is_file())
    for required in ('constant/MRFProperties', 'constant/phaseProperties',
                     'system/controlDict', 'system/fvSchemes', 'system/fvSolution'):
        if not (case / required).is_file() and not (case / (required + '.gz')).is_file():
            issues.append(f'{name}: missing required configuration {required}')
    with (out / f'{name}_configuration.txt').open('x') as report:
        for path in sorted(paths):
            try:
                data, meta = read_file(path)
                text = data.decode('utf-8')
                report.write(f'\nFILE {path}\nSHA256 {meta["sha256"]}\n{text}\n')
                records.append(meta)
                if '#include' in text or '#codeStream' in text:
                    issues.append(f'{name}: literal directive requires review in {path.relative_to(case)}')
            except (OSError, ValueError, UnicodeError) as exc:
                issues.append(f'{name}: configuration read failed: {path}: {exc}')
        available = [str(p.name) for p in case.iterdir() if p.is_dir() and
                     (p.name in ('0', '0.orig', '0.org', 'orig0') or numeric_time(p.name))]
        report.write('\nAVAILABLE ROOT/INITIAL DIRECTORIES: ' + ', '.join(sorted(available)) + '\n')
        # Initial templates plus the compared early checkpoint. Do not open 31.2 s runs.
        for time in ('0', '0.orig', '0.org', 'orig0', '9.9'):
            directory = case / time
            if not directory.is_dir():
                continue
            fields = sorted({p.with_suffix('') if p.suffix == '.gz' else p
                             for p in directory.iterdir() if p.is_file() and
                             (p.name.startswith(('U.', 'alpha.')) or p.name in ('p', 'p_rgh', 'p.gz', 'p_rgh.gz'))})
            for path in fields:
                report.write(f'\nBOUNDARY FIELD {path}\n')
                if not foam_dictionary:
                    issues.append(f'{name}: foamDictionary unavailable; cannot read {path}')
                    continue
                try:
                    _, meta = read_file(path)
                    # Foundation OF9 disables function entries by default without -expand.
                    # No mutation options, no libraries, and cwd is the new report directory.
                    result = subprocess.run([foam_dictionary, str(path), '-entry', 'boundaryField'],
                                            cwd=out, capture_output=True, text=True, timeout=60)
                    _, after = read_file(path)
                    if meta['sha256'] != after['sha256']:
                        raise ValueError(f'Input changed while querying boundaryField: {path}')
                    report.write(result.stdout + result.stderr + f'\nstatus={result.returncode}\n')
                    records.append(meta)
                    if result.returncode:
                        issues.append(f'{name}: boundaryField query failed: {path}')
                except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
                    issues.append(f'{name}: boundaryField query failed: {path}: {exc}')
        report.write('\nIncludes/substitutions are not expanded. No t=0 volume or integrated flux is inferred.\n')
    return records


def numeric_time(name):
    try:
        return np.isfinite(float(name))
    except ValueError:
        return False


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--emm', type=Path, default=ARCHIVE / 'emm/putida_sc1')
    parser.add_argument('--meml', type=Path, default=ARCHIVE / 'meml/putida_sc1_final')
    parser.add_argument('--vtk-root', type=Path, default=VTK_ROOT)
    parser.add_argument('--zone', default='rotor')
    args = parser.parse_args(argv)
    cases = dict(EMM=args.emm.resolve(), MEML_BASE=args.meml.resolve())
    prepare_output(args.out, [*cases.values(), args.vtk_root])
    out = args.out.resolve()
    issues, summary = [], ['SC1 wall/MRF evidence collection (no CFD launched)',
                          'Zone adjacency is not a verdict on stationary-wall configuration.',
                          'Sector crossings are volume averages; films/droplets are not classified as connected surfaces.']
    report = dict(cases={}, vtk_root=str(args.vtk_root.resolve()), zone=args.zone,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  mesh_reader_sha256=hashlib.sha256(Path(__file__).with_name('sc1_wall_mesh.py').read_bytes()).hexdigest())
    foam_dictionary = shutil.which('foamDictionary')
    if os.environ.get('WM_PROJECT_VERSION') != '9':
        issues.append('Foundation OF9 environment not identified; foamDictionary queries disabled')
        foam_dictionary = None
    for name, case in cases.items():
        print(f'Reading {name}: {case}', flush=True)
        entry = report['cases'][name] = dict(path=str(case), surfaces={})
        summary.append(f'\n{name}: {case}')
        try:
            entry['configuration_inputs'] = collect_configuration(case, name, out, foam_dictionary, issues)
        except (OSError, ValueError) as exc:
            issues.append(f'{name}: configuration collection incomplete: {exc}')
        try:
            mesh = zone_adjacency(case / 'constant/polyMesh', args.zone)
            entry['zone_adjacency'] = mesh
            summary.append(f'Zone {args.zone}: {mesh["zone_cells"]} cells')
            for patch in mesh['patches']:
                summary.append(f'  {patch["patch"]} ({patch["type"]}): '
                               f'{patch["zone_adjacent_faces"]}/{patch["faces"]} faces adjacent to zone')
        except (OSError, ValueError) as exc:
            issues.append(f'{name}: mesh adjacency unavailable: {exc}')
        for time in ('9.9', '20'):
            try:
                # Time-local mesh files would require their own addressing audit.
                if (case / time / 'polyMesh').exists():
                    issues.append(f'{name}/{time}: time-local mesh present; constant-mesh adjacency is insufficient')
                vtk = existing_vtk(args.vtk_root, name, time)
                print(f'Analyzing existing {vtk}', flush=True)
                data = evaluate_vtk(vtk)
                entry['surfaces'][time] = data
                summary.append(f'\n{name} t={time} s ({vtk})')
                summary.extend(surface_summary(data))
            except (OSError, ValueError, ImportError, KeyError, AttributeError) as exc:
                issues.append(f'{name}/{time}: surface analysis unavailable: {exc}')
    summary += ['\nREVIEW ITEMS (not a physical pass/fail):', *issues,
                't=0 inventory and integrated historical flux/source budgets remain unestablished.',
                'No existing run was altered or repeated. Review literal includes and all patch purposes.',
                f'Full report: {out / "report.json"}',
                f'Boundary/configuration evidence: {out / "EMM_configuration.txt"}, {out / "MEML_BASE_configuration.txt"}']
    report['review_items'] = issues
    with (out / 'report.json').open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (out / 'summary.txt').open('x') as stream:
        stream.write('\n'.join(summary) + '\n')
    with tarfile.open(out / 'evidence.tar.gz', 'x:gz') as archive:
        for name in ('summary.txt', 'report.json', 'EMM_configuration.txt', 'MEML_BASE_configuration.txt'):
            if (out / name).is_file():
                archive.add(out / name, arcname=name, recursive=False)
    # Print MRF dictionary evidence with the summary; the remaining configuration
    # is retained for upload rather than flooding the allocation terminal.
    for name, case in cases.items():
        try:
            data, _ = read_file(case / 'constant/MRFProperties')
            print(f'\n{name} literal MRFProperties:\n{data.decode("utf-8")}', flush=True)
        except (OSError, ValueError, UnicodeError):
            pass
    print('\n'.join(summary), flush=True)
    print(f'Evidence bundle: {out / "evidence.tar.gz"}', flush=True)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f'SC1 wall preflight stopped: {error}', file=sys.stderr)
        sys.exit(1)
