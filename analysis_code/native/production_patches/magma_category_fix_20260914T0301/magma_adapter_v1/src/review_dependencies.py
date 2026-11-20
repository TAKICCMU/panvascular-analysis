"""Content-bound approval dependencies, separate from cache identity."""
from pathlib import Path
import hashlib,json,csv
from .core import Blocked,sha

SKIP_KEYS={'approval_receipt','root','release_support','output','out','outdir','output_dir','rscript','executable','python','checksum','input_sha256'}

def reviewed_dependencies(config, config_path):
    base=Path(config_path).resolve().parent
    files={}; seen=set()
    compat_manifest=None; metadata_path=None; bound_evidence={}
    if config.get('atlas',{}).get('release_compatibility_required'):
        # Historical provenance JSONs contain retired staging paths. Bind their
        # actual bytes, not every historical locator as a live runtime resource.
        compat_manifest=Path(config['atlas']['manifest']).resolve()
        manifest=json.loads(compat_manifest.read_text())
        support_path=Path(config['outputs']['release_support'])
        metadata_path=support_path/'metadata_review_sources.json'
        if sha(metadata_path)!=manifest.get('release_sidecar_sha256',{}).get(metadata_path.name):raise Blocked('Metadata source index binding differs')
        records=json.loads(metadata_path.read_text()).get('files',[])+[manifest.get('original_manifest',{}),manifest.get('publisher_script',{})]
        registry_path=support_path/'source_registry.tsv'
        if sha(registry_path)!=manifest.get('release_sidecar_sha256',{}).get(registry_path.name):raise Blocked('Source registry dependency binding differs')
        with registry_path.open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                for path_key,hash_key in [('version_source','version_source_sha256'),('mapping_path','mapping_sha256'),('label_usage_path','label_usage_sha256')]:
                    if row.get(path_key):records.append({'path':row[path_key],'sha256':row.get(hash_key)})
        for record in records:
            if not record.get('path') or not record.get('sha256'):raise Blocked('Missing explicit provenance dependency binding')
            key=str(Path(record['path']).resolve())
            if key in bound_evidence and bound_evidence[key]!=record['sha256']:raise Blocked('Conflicting provenance dependency hashes')
            bound_evidence[key]=record['sha256']
    def inspect(path,origin):
        path=Path(path).expanduser()
        if not path.is_absolute():path=base/path
        path=path.resolve();s=str(path)
        if s in seen:return
        seen.add(s)
        if not path.is_file():raise Blocked(f'Review dependency missing: {origin}: {path}')
        files[s]={'sha256':sha(path),'origin':origin}
        if s in bound_evidence:
            if files[s]['sha256']!=bound_evidence[s]:raise Blocked('Provenance dependency SHA differs: '+s)
            return
        if compat_manifest is not None and path==Path(config['atlas']['schema_contract']).resolve():
            # Schema includes example partial filenames and write-policy templates,
            # not live inputs. Its full bytes remain in the signed dependency map.
            return
        if metadata_path is not None and path==metadata_path.resolve():
            # Only explicit file records are live resources; URI/bundle-member
            # locators in the provenance narrative are not filesystem inputs.
            for record in json.loads(path.read_text())['files']:inspect(record['path'],s+'.files')
            return
        if compat_manifest is not None and path==compat_manifest:
            # Complete manifest bytes are bound. Live sidecars and their actual
            # source objects are traversed below; source/backup paths here are history.
            manifest=json.loads(path.read_text())
            for key in ['original_manifest','publisher_script']:walk(manifest[key],s+'.'+key)
            return
        if path.suffix.lower() in {'.json','.yaml','.yml'}:
            try:
                import yaml
                walk(yaml.safe_load(path.read_text()),s)
            except (UnicodeDecodeError,ValueError):raise Blocked(f'Invalid dependency manifest: {path}')
        elif path.suffix.lower() in {'.tsv','.csv'}:
            with path.open() as f:
                for row in csv.DictReader(f,delimiter='\t' if path.suffix=='.tsv' else ','):walk(row,s)
    def walk(obj,origin,key=''):
        if key in SKIP_KEYS:return
        if isinstance(obj,dict):
            for k,v in obj.items():walk(v,origin+'.'+str(k),str(k))
        elif isinstance(obj,list):
            for i,v in enumerate(obj):walk(v,origin+f'[{i}]',key)
        elif isinstance(obj,str):
            if obj.startswith(('/', './', '../')) or Path(obj).suffix.lower() in {'.tsv','.csv','.json','.yaml','.yml','.gz','.txt','.rds','.gs','.gmt','.raw','.h5ad','.npz'}:
                # Atlas content is already independently bound and verified by preflight.
                if obj==config.get('atlas',{}).get('path'):return
                candidate=Path(obj) if Path(obj).is_absolute() else base/obj
                if candidate.is_dir():return
                if candidate.suffix.lower() in {'.tsv','.csv','.json','.yaml','.yml','.gz','.txt','.rds','.gs','.gmt','.raw','.h5ad','.npz'} or candidate.is_file():inspect(candidate,origin)
    walk(config,'config')
    support=config.get('outputs',{}).get('release_support')
    if support:
        support=Path(support)
        if not support.is_dir():raise Blocked('Release support directory missing')
        for f in sorted(support.iterdir()):
            if f.is_file():inspect(f,'release_support')
    payload={'schema':'reviewed_dependencies_v1','files':dict(sorted(files.items()))}
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return digest,payload
