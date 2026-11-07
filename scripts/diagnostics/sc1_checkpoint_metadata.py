"""Prove a narrow outer-header path difference without changing any file.

This does not decode OpenFOAM fields or broadly normalize headers. Every byte
outside one validated location value must agree, including the full binary
tail and any embedded per-rank headers.
"""
import hashlib
from pathlib import Path
import re


COURANT_FILES=frozenset(('Co.liquid','Co.bubbleGas','Co.headspaceGas','Co.mixture'))
HEADER_LIMIT=8192
TOKEN=re.compile(rb'\s+|/\*.*?\*/|//[^\r\n]*|"(?:\\.|[^"\\])*"|[{};]|[^\s{};"]+',re.S)


def outer_header(stream, path):
    prefix=stream.read(HEADER_LIMIT)
    position=0

    def take():
        nonlocal position
        while position<len(prefix):
            match=TOKEN.match(prefix,position)
            if match is None:
                raise ValueError(f'Unrecognized outer header token: {path}')
            position=match.end()
            value=match[0]
            if value.startswith(b'/*') and not value.endswith(b'*/'):
                raise ValueError(f'Unterminated outer header comment: {path}')
            if value.isspace() or value.startswith((b'/*',b'//')):
                continue
            return value,match.start(),match.end()
        raise ValueError(f'Missing or oversized outer FoamFile header: {path}')

    if take()[0]!=b'FoamFile' or take()[0]!=b'{':
        raise ValueError(f'Expected outer FoamFile header: {path}')
    entries={}
    while True:
        key=take()[0]
        if key==b'}':
            break
        if not re.fullmatch(rb'[A-Za-z_][A-Za-z_0-9]*',key) or key in entries:
            raise ValueError(f'Invalid or duplicate header key {key!r}: {path}')
        values=[]
        while True:
            token=take()
            if token[0]==b';':
                break
            if token[0] in (b'{',b'}'):
                raise ValueError(f'Non-flat or unterminated outer header: {path}')
            values.append(token)
        if not values:
            raise ValueError(f'Empty outer header entry {key!r}: {path}')
        entries[key]=values

    for key,expected in ((b'class',b'decomposedBlockData'),(b'format',b'binary'),
                         (b'object',path.name.encode('ascii'))):
        values=entries.get(key,[])
        if len(values)!=1 or values[0][0]!=expected:
            raise ValueError(f'Unexpected outer header {key.decode()}: {path}')
    values=entries.get(b'location',[])
    if len(values)!=1 or not values[0][0].startswith(b'"') or b'\\' in values[0][0]:
        raise ValueError(f'Expected one unescaped quoted location: {path}')
    value,start,end=values[0]
    location=value[1:-1].decode('ascii')
    if not Path(location).is_absolute() or Path(location).resolve()!=path.resolve():
        raise ValueError(f'Outer location does not identify its source file: {path}')
    header=prefix[:position]
    normalized=header[:start+1]+b'<SC1_SOURCE_FILE_LOCATION>'+header[end-1:]
    return header,normalized,location


def verify_courant_location_metadata(control_file, normal_file):
    """Return proof only for an allowed name and complete remaining-byte match."""
    control_file,normal_file=Path(control_file),Path(normal_file)
    if control_file.name not in COURANT_FILES or normal_file.name!=control_file.name:
        raise ValueError('Metadata exception is restricted to the four observed Co files')
    before=[p.stat() for p in (control_file,normal_file)]
    with control_file.open('rb') as a,normal_file.open('rb') as b:
        header_a,normalized_a,location_a=outer_header(a,control_file)
        header_b,normalized_b,location_b=outer_header(b,normal_file)
        if normalized_a!=normalized_b:
            raise ValueError(f'Co header differs beyond its location value: {control_file.name}')
        digest_a=hashlib.sha256(header_a)
        digest_b=hashlib.sha256(header_b)
        payload=hashlib.sha256()
        size=0
        a.seek(len(header_a)); b.seek(len(header_b))
        while True:
            x,y=a.read(1024*1024),b.read(1024*1024)
            if x!=y:
                raise ValueError(f'Co bytes differ after the outer header: {control_file.name}')
            if not x:
                break
            digest_a.update(x); digest_b.update(y); payload.update(x)
            size+=len(x)
    after=[p.stat() for p in (control_file,normal_file)]
    identity=lambda s:(s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev)
    if any(identity(x)!=identity(y) for x,y in zip(before,after)):
        raise ValueError('Co input changed during metadata verification')
    return {'control_sha256':digest_a.hexdigest(),'normal_sha256':digest_b.hexdigest(),
            'control_location':location_a,'normal_location':location_b,
            'normalized_header_sha256':hashlib.sha256(normalized_a).hexdigest(),
            'payload_sha256':payload.hexdigest(),'payload_bytes':size}
