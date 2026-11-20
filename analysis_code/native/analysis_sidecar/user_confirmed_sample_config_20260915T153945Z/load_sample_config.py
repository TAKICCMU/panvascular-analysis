import csv, hashlib, json
from pathlib import Path

def load_sample_registry(config_path):
    cfg = json.loads(Path(config_path).read_text())
    path = Path(cfg['sample_registry'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != cfg['sample_registry_sha256']:
        raise ValueError('Sample registry checksum mismatch')
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len({r['sample_uid'] for r in rows}) != len(rows):
        raise ValueError('Duplicate sample_uid')
    return {r['sample_uid']: r for r in rows}
