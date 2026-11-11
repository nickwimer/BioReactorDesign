"""Read-only, restricted Foundation OpenFOAM mesh addressing reader.

Supports generated ASCII/binary owner, cellZones and boundary files, including
gzip and declared 32/64-bit labels. Rejects unsupported dictionary directives.
It reports zone adjacency, not whether MRF considers a patch stationary.
"""
import gzip
import hashlib
from pathlib import Path
import re

import numpy as np


TOKEN = re.compile(rb'\s+|/\*.*?\*/|//[^\r\n]*|"(?:\\.|[^"\\])*"|[{}();]|[^\s{}();"]+', re.S)


def read_file(path):
    path = Path(path)
    if not path.exists() and path.with_name(path.name + '.gz').is_file():
        path = path.with_name(path.name + '.gz')
    before = path.stat()
    raw = path.read_bytes()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (
            after.st_size, after.st_mtime_ns, after.st_ino):
        raise ValueError(f'Input changed during read: {path}')
    meta = dict(path=str(path.resolve()), bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest())
    return (gzip.decompress(raw) if path.suffix == '.gz' else raw), meta


class Reader:
    def __init__(self, data):
        self.data, self.pos = data, 0
        self.expect(b'FoamFile')
        self.header = self.dictionary()
        self.binary = self.header.get('format') == [b'binary']
        if self.header.get('format') not in ([b'ascii'], [b'binary']):
            raise ValueError('Missing/unsupported FoamFile format')
        self.dtype = None
        if self.binary:
            arch = self.header.get('arch', [])
            if len(arch) != 1:
                raise ValueError('Binary addressing requires explicit arch')
            match = re.fullmatch(rb'"(LSB|MSB);label=(32|64);scalar=(32|64)"', arch[0])
            if not match:
                raise ValueError('Unsupported binary arch')
            self.dtype = np.dtype(('<' if match[1] == b'LSB' else '>')
                                  + 'i' + str(int(match[2]) // 8))

    def take(self):
        while self.pos < len(self.data):
            m = TOKEN.match(self.data, self.pos)
            if m is None:
                raise ValueError(f'Unrecognized token at {self.pos}')
            self.pos = m.end()
            t = m[0]
            if t.isspace() or t.startswith((b'/*', b'//')):
                if t.startswith(b'/*') and not t.endswith(b'*/'):
                    raise ValueError('Unterminated comment')
                continue
            if t.startswith((b'#', b'$')):
                raise ValueError('Unexpanded directive/substitution in mesh addressing')
            return t
        return None

    def expect(self, value):
        actual = self.take()
        if actual != value:
            raise ValueError(f'Expected {value!r}, found {actual!r} at {self.pos}')

    def count(self, token=None):
        token = self.take() if token is None else token
        if token is None or not re.fullmatch(rb'\d+', token):
            raise ValueError(f'Expected nonnegative integer, found {token!r}')
        return int(token)

    def dictionary(self, labels=False):
        self.expect(b'{')
        result = {}
        while True:
            key = self.take()
            if key == b'}':
                return result
            if key is None:
                raise ValueError('Unterminated dictionary')
            key = key.decode('ascii').strip('"')
            if key in result:
                raise ValueError(f'Duplicate key: {key}')
            if labels and key == 'cellLabels':
                result[key] = self.label_list()
                self.expect(b';')
                continue
            values, depth = [], 0
            while True:
                t = self.take()
                if t is None or t in (b'{', b'}'):
                    raise ValueError(f'Unsupported nested/malformed entry: {key}')
                if t == b';' and depth == 0:
                    break
                if t == b'(':
                    depth += 1
                elif t == b')':
                    depth -= 1
                if depth < 0:
                    raise ValueError('Malformed list')
                values.append(t)
            result[key] = values

    def label_list(self):
        t = self.take()
        if t == b'List<label>':
            t = self.take()
        n = self.count(t)
        if self.binary:
            # Foundation binary UList emits no parentheses for an empty list.
            if n == 0:
                return np.empty(0, dtype=np.int64)
            self.expect(b'(')
            # Raw data starts immediately after '('; whitespace bytes can be labels.
            end = self.pos + n * self.dtype.itemsize
            if end > len(self.data):
                raise ValueError('Truncated binary label list')
            a = np.frombuffer(self.data, dtype=self.dtype, count=n, offset=self.pos).copy()
            self.pos = end
            self.expect(b')')
            return a
        opening = self.take()
        if opening == b'{':
            value = self.count()
            self.expect(b'}')
            return np.full(n, value, dtype=np.int64)
        if opening != b'(':
            raise ValueError('Expected ASCII label list')
        # Mesh labels can be large; parse the numeric payload in one operation.
        end = self.data.find(b')', self.pos)
        if end < 0:
            raise ValueError('Unterminated ASCII label list')
        payload = self.data[self.pos:end]
        payload = re.sub(rb'/\*.*?\*/|//[^\r\n]*', b' ', payload, flags=re.S)
        if not re.fullmatch(rb'[\s\d+-]*', payload):
            raise ValueError('Nonnumeric ASCII label list')
        tokens = payload.split()
        if any(not re.fullmatch(rb'\d+', t) for t in tokens):
            raise ValueError('Invalid/negative mesh labels')
        a = np.fromstring(payload.decode('ascii'), dtype=np.int64, sep=' ')
        if len(a) != n:
            raise ValueError('ASCII label count mismatch')
        self.pos = end + 1
        return a

    def named_list(self, labels=False):
        n = self.count()
        self.expect(b'(')
        result = {}
        for _ in range(n):
            name = self.take()
            if name is None:
                raise ValueError('Missing named-list entry')
            name = name.decode('ascii').strip('"')
            if name in result:
                raise ValueError('Duplicate named-list entry')
            result[name] = self.dictionary(labels=labels)
        self.expect(b')')
        return result

    def finish(self):
        if self.take() is not None:
            raise ValueError('Unexpected trailing mesh data')


def zone_adjacency(mesh_dir, zone_name='rotor'):
    """Count boundary faces whose owner cell belongs to a generated cell zone."""
    readers, metadata = {}, {}
    classes = dict(owner=b'labelList', neighbour=b'labelList',
                   cellZones=b'regIOobject', boundary=b'polyBoundaryMesh')
    for name in classes:
        data, metadata[name] = read_file(Path(mesh_dir) / name)
        readers[name] = Reader(data)
        if readers[name].header.get('object') != [name.encode()]:
            raise ValueError(f'Wrong mesh object: {name}')
        if readers[name].header.get('class') != [classes[name]]:
            raise ValueError(f'Unsupported mesh class: {name}')
    owner = readers['owner'].label_list()
    neighbour = readers['neighbour'].label_list()
    zones = readers['cellZones'].named_list(labels=True)
    boundary = readers['boundary'].named_list()
    for reader in readers.values():
        reader.finish()
    if zone_name not in zones:
        raise ValueError(f'Zone {zone_name!r} not present; available: {list(zones)}')
    z = zones[zone_name]
    if z.get('type') != [b'cellZone'] or 'cellLabels' not in z:
        raise ValueError('Unsupported/incomplete cell zone')
    cells = z['cellLabels']
    if not owner.size or np.any(owner < 0) or np.any(neighbour < 0) or np.any(cells < 0):
        raise ValueError('Invalid owner/zone addressing')
    max_cell = max(owner.max(), neighbour.max() if neighbour.size else -1)
    if len(np.unique(cells)) != len(cells) or (len(cells) and cells.max() > max_cell):
        raise ValueError('Duplicate/out-of-range zone cells')
    # Membership lookup without allocating max(label) entries.
    rows, previous_end = [], len(neighbour)
    for name, patch in boundary.items():
        def integer(key):
            value = patch.get(key, [])
            if len(value) != 1 or not re.fullmatch(rb'\d+', value[0]):
                raise ValueError(f'Missing/invalid {name}/{key}')
            return int(value[0])
        start, count = integer('startFace'), integer('nFaces')
        if start + count > len(owner) or start != previous_end:
            raise ValueError('Boundary addressing outside owner list or noncontiguous')
        previous_end = start + count
        touched = int(np.count_nonzero(np.isin(owner[start:start + count], cells)))
        rows.append(dict(patch=name, type=b' '.join(patch.get('type', [])).decode(),
                         faces=count, zone_adjacent_faces=touched))
    if previous_end != len(owner):
        raise ValueError('Boundary does not terminate at owner-list end')
    return dict(zone=zone_name, zone_cells=len(cells), patches=rows, inputs=metadata,
                interpretation='Adjacency only; compare actual MRF exclusions and patch purpose.')
