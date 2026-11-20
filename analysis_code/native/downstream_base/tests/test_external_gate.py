import unittest,tempfile,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from external_gate import load_config,approved_preflight
class GateTests(unittest.TestCase):
 def test_current_config_cannot_execute_or_create_output(self):
  from src.core import Blocked
  with tempfile.TemporaryDirectory() as d:
   output=Path(d)/'no_scientific_output'
   with self.assertRaises(Blocked):approved_preflight(ROOT/'configs/analysis_config.yaml',output)
   self.assertFalse(output.exists())
 def test_cpu_budget_applied(self):
  c,p=load_config(ROOT/'configs/analysis_config.yaml')
  self.assertEqual(os.environ['CUDA_VISIBLE_DEVICES'],'')
  self.assertLessEqual(int(os.environ['OMP_NUM_THREADS']),20)
if __name__=='__main__':unittest.main()
