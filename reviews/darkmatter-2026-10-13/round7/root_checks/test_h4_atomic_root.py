from pathlib import Path
import sys,tempfile,unittest,json
from unittest.mock import patch
sys.path.insert(0,str(Path('reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs').resolve()))
import v7_crop_and_search_regions as m
class AtomicTests(unittest.TestCase):
 def test_zero_hit_signed_reuse_and_stale_rejection(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);q=p/'q';q.mkdir();(q/'one.pdb').write_text('query')
   out=p/'hits';tmp=p/'scratch'
   def run(cmd,check):Path(cmd[4]).write_text('')
   args=(str(q),'db',str(out),str(tmp),1,5000)
   with patch.object(m,'verified_database_identity',return_value={'component':'hash'}),patch.object(m.subprocess,'check_output',return_value='v10'),patch.object(m.subprocess,'run',side_effect=run) as runner:
    m.run_foldseek_search_atomic(*args,reuse=True)
    m.run_foldseek_search_atomic(*args,reuse=True)
    self.assertEqual(runner.call_count,1)
    (q/'one.pdb').write_text('changed')
    with self.assertRaises(AssertionError):m.run_foldseek_search_atomic(*args,reuse=True)
    (q/'one.pdb').write_text('query')
    with patch.object(m,'verified_database_identity',return_value={'component':'changed'}):
     with self.assertRaises(AssertionError):m.run_foldseek_search_atomic(*args,reuse=True)
    sig=Path(str(out)+'.signature.json');sig.unlink()
    with self.assertRaises(AssertionError):m.run_foldseek_search_atomic(*args,reuse=True)
if __name__=='__main__':unittest.main()
