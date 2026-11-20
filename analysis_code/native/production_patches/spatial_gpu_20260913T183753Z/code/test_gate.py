import tempfile,shutil,json,os,sys,hashlib,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
import spatial_gpu as s
class GateTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.old=s.JOB;s.JOB=Path(self.tmp.name);(s.JOB/'configs').mkdir();self.original_env=os.environ.copy()
  for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
  os.environ['CUDA_VISIBLE_DEVICES']='GPU-synthetic-test'
  a=s.JOB/'synthetic_atlas';a.write_bytes(b'controlled tiny synthetic test');st=a.stat()
  c={'atlas':str(a),'atlas_sha256':'00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310','gpu_uuid':'GPU-synthetic-test','reference_mode':'descriptive','forced_genes':[],'excluded_targets':[]}
  p=s.JOB/'configs/config.json';p.write_text(json.dumps(c))
  r={'scope':['descriptive_reference','descriptive_mapping'],'cpu_threads':4,'bound_files':{str(p):s.sha(p)},'atlas_stat':[st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns],'atlas_sha256':c['atlas_sha256']};(s.JOB/'configs/PRODUCTION_RECEIPT.json').write_text(json.dumps(r))
 def tearDown(self):s.JOB=self.old;self.tmp.cleanup();os.environ.clear();os.environ.update(self.original_env)
 def test_valid_content_binding(self):self.assertEqual(s.guard()['reference_mode'],'descriptive')
 def test_config_drift(self):
  with open(s.JOB/'configs/config.json','a') as f:f.write(' ')
  with self.assertRaises(AssertionError):s.guard()
 def test_atlas_stat_drift(self):
  with open(s.JOB/'synthetic_atlas','ab') as f:f.write(b'x')
  with self.assertRaises(AssertionError):s.guard()
 def test_hidden_cuda(self):
  os.environ['CUDA_VISIBLE_DEVICES']=''
  with self.assertRaises(AssertionError):s.guard()
 def test_cpu_budget(self):
  os.environ['OMP_NUM_THREADS']='40'
  with self.assertRaises(AssertionError):s.guard()
if __name__=='__main__':unittest.main()
