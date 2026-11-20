from pathlib import Path
import json
MAP=json.loads((Path(__file__).resolve().parents[1]/'merge_records/PATH_MAP.json').read_text())
def resolve_legacy_path(value):
    value=str(value)
    for old,new in sorted(MAP.items(),key=lambda x:-len(x[0])):
        value=value.replace(old,new)
    return value
