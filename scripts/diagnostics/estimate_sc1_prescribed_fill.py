#!/usr/bin/env python3
"""Estimate archived setFields liquid inventories using already exported meshes.

No OpenFOAM application is invoked. This is not a measured t=0 inventory.
VTK centers need not equal OpenFOAM finite-volume centroids, so also report
conservative box-selection bounds using each cell's bounding box. These bounds
assume valid cells whose finite-volume centroids lie within those bounds.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from collect_sc1_wall_preflight import prepare_output
from sc1_wall_mesh import Reader, read_file


def parse_fill(data):
    """Accept only scalar defaults and ordered boxToCell scalar overrides."""
    r = Reader(data)
    if r.binary:
        raise ValueError('setFieldsDict must be ASCII')

    def number():
        token = r.take()
        if token is None:
            raise ValueError('Missing scalar')
        value = float(token)
        if not np.isfinite(value):
            raise ValueError('Nonfinite scalar')
        return value

    def vector():
        r.expect(b'(')
        v = [number() for _ in range(3)]
        r.expect(b')')
        return v

    def fields():
        r.expect(b'(')
        values = {}
        while True:
            t = r.take()
            if t == b')':
                break
            if t != b'volScalarFieldValue':
                raise ValueError(f'Unsupported setFields entry: {t!r}')
            field_token = r.take()
            if field_token is None:
                raise ValueError('Missing scalar field name')
            field = field_token.decode('ascii')
            if field in values:
                raise ValueError('Duplicate scalar field')
            values[field] = number()
        r.expect(b';')
        return values

    r.expect(b'defaultFieldValues')
    default = fields()
    if 'alpha.liquid' not in default:
        raise ValueError('No explicit liquid default')
    r.expect(b'regions')
    r.expect(b'(')
    regions = []
    while True:
        t = r.take()
        if t == b')':
            break
        if t != b'boxToCell':
            raise ValueError(f'Unsupported selection: {t!r}')
        r.expect(b'{')
        r.expect(b'box')
        lower, upper = vector(), vector()
        if np.any(np.array(lower) > upper):
            raise ValueError('Inverted box')
        r.expect(b';')
        r.expect(b'fieldValues')
        values = fields()
        r.expect(b'}')
        regions.append(dict(lower=lower, upper=upper, values=values))
    r.expect(b';')
    r.finish()
    return dict(default=default, regions=regions)


def estimate(centers, bounds, volumes, prescription):
    centers, bounds, volumes = [np.asarray(a, dtype=float) for a in (centers, bounds, volumes)]
    n = len(volumes)
    if centers.shape != (n, 3) or bounds.shape != (n, 6) or volumes.shape != (n,) or not n:
        raise ValueError('Invalid mesh array shapes')
    if any(not np.all(np.isfinite(a)) for a in (centers, bounds, volumes)) or np.any(volumes <= 0):
        raise ValueError('Invalid mesh geometry/volumes')
    cell_lo, cell_hi = bounds[:, ::2], bounds[:, 1::2]
    if np.any(cell_lo > cell_hi) or np.any(centers < cell_lo) or np.any(centers > cell_hi):
        raise ValueError('Invalid cell bounds/centers')
    alpha = np.full(n, prescription['default']['alpha.liquid'])
    minimum, maximum = alpha.copy(), alpha.copy()
    region_counts = []
    for region in prescription['regions']:
        if 'alpha.liquid' not in region['values']:
            continue
        lo, hi = np.array(region['lower']), np.array(region['upper'])
        value = region['values']['alpha.liquid']
        chosen = np.all((centers >= lo) & (centers <= hi), axis=1)
        certain = np.all((cell_lo >= lo) & (cell_hi <= hi), axis=1)
        possible = np.all((cell_hi >= lo) & (cell_lo <= hi), axis=1)
        uncertain = possible & ~certain
        alpha[chosen] = value
        minimum[certain] = maximum[certain] = value
        minimum[uncertain] = np.minimum(minimum[uncertain], value)
        maximum[uncertain] = np.maximum(maximum[uncertain], value)
        region_counts.append(dict(center_selected_cells=int(chosen.sum()),
                                  box_selection_uncertain_cells=int(uncertain.sum())))
    if np.any(minimum < 0) or np.any(maximum > 1):
        raise ValueError('Prescribed liquid fraction is outside [0,1]')
    return dict(cells=n, mesh_bounds_m=[float(x) for pair in zip(cell_lo.min(axis=0), cell_hi.max(axis=0)) for x in pair],
                mesh_volume_L=float(1000 * volumes.sum()),
                prescribed_liquid_center_estimate_L=float(1000 * np.dot(alpha, volumes)),
                prescribed_liquid_box_selection_bounds_L=[float(1000 * np.dot(minimum, volumes)),
                                                         float(1000 * np.dot(maximum, volumes))],
                region_counts=region_counts)


def verify_vtk(path, expected):
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    if h.hexdigest() != expected['sha256'] or before.st_size != expected['bytes']:
        raise ValueError(f'VTK differs from completed preflight: {path}')
    return before.st_size, before.st_mtime_ns, before.st_ino


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    report_data = args.report.read_bytes()
    report = json.loads(report_data)
    inputs = [args.report.parent]
    for name in ('EMM', 'MEML_BASE'):
        inputs += [Path(report['cases'][name]['path']),
                   Path(report['cases'][name]['surfaces']['9.9']['input']['path']).parent]
    prepare_output(args.out, inputs)
    import pyvista as pv
    output = dict(report_sha256=hashlib.sha256(report_data).hexdigest(), cases={},
                  meaning='Prescribed-fill estimate, not measured t=0 or a closed conservation budget.',
                  uncertainty='Bounding-box selection bounds assume valid cell centroids inside each cell bounding box; VTK integration error is not included.')
    lines = [output['meaning'], output['uncertainty']]
    for name in ('EMM', 'MEML_BASE'):
        case = report['cases'][name]
        expected = [m for m in case['configuration_inputs'] if m['path'].endswith('/system/setFieldsDict')]
        if len(expected) != 1:
            raise ValueError(f'{name}: missing/ambiguous setFields provenance')
        data, meta = read_file(expected[0]['path'])
        if meta != expected[0]:
            raise ValueError(f'{name}: setFieldsDict changed since preflight')
        prescription = parse_fill(data)
        vtk_meta = case['surfaces']['9.9']['input']
        vtk = Path(vtk_meta['path'])
        before = verify_vtk(vtk, vtk_meta)
        print(f'Checking prescribed fill on existing {name} mesh: {vtk}', flush=True)
        mesh = pv.read(str(vtk))
        volumes = mesh.compute_cell_sizes(length=False, area=False, volume=True).cell_data['Volume']
        bounds = np.empty((mesh.n_cells, 6))
        for i in range(mesh.n_cells):
            mesh.GetCellBounds(i, bounds[i])
        result = estimate(mesh.cell_centers().points, bounds, volumes, prescription)
        if result['cells'] != case['surfaces']['9.9']['cells']:
            raise ValueError(f'{name}: VTK cell count changed')
        after = vtk.stat()
        if before != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise ValueError(f'{name}: VTK changed during reading')
        result.update(prescription=prescription, set_fields_input=meta, vtk_input=vtk_meta,
                      measured_9_9_liquid_L=case['surfaces']['9.9']['total_liquid_L'])
        output['cases'][name] = result
        lines.append(f"{name}: prescribed center estimate={result['prescribed_liquid_center_estimate_L']:.6f} L; "
                     f"box-selection bounds={result['prescribed_liquid_box_selection_bounds_L']}; "
                     f"measured 9.9 s={result['measured_9_9_liquid_L']:.6f} L; "
                     f"mesh bounds={result['mesh_bounds_m']}")
    lines += ['No fields or dictionaries were changed. No conversion, compilation or CFD was run.',
              f'REPORT={args.out.resolve()}']
    with (args.out / 'report.json').open('x') as f:
        json.dump(output, f, indent=2, allow_nan=False)
        f.write('\n')
    with (args.out / 'summary.txt').open('x') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines), flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, ImportError) as exc:
        print(f'SC1 prescribed-fill check stopped: {exc}', file=sys.stderr)
        sys.exit(1)
