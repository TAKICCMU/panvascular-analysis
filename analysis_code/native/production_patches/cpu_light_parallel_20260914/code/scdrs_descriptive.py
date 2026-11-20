#!/usr/bin/env python3
import os
for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"]:os.environ[k]="1"
os.environ["CUDA_VISIBLE_DEVICES"]=""
from pathlib import Path
import json,sys,time,hashlib,gzip
import numpy as np,pandas as pd
W=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(8*1024**2),b""):h.update(b)
 return h.hexdigest()
def validate(df,ids):
 cols=["raw_score","norm_score","mc_pval","pval","nlog10_pval","zscore"]
 if not df.index.is_unique or not ids.is_unique:raise ValueError("Duplicate cell IDs")
 if not df.index.equals(ids):raise ValueError("Cell IDs or order differ from certified atlas")
 if not set(cols)<=set(df):raise ValueError("Required native columns missing")
 if not np.isfinite(df[cols].to_numpy()).all():raise ValueError("Nonfinite native score")
 if not df[["mc_pval","pval"]].ge(0).all().all() or not df[["mc_pval","pval"]].le(1).all().all():raise ValueError("P outside [0,1]")
 return df
def read_col(group):
 if hasattr(group,"shape"):
  if group.dtype.kind in "OSU":return group.asstr()[:]
  return group[:]
 cats=read_col(group["categories"]);codes=group["codes"][:]
 return pd.Categorical.from_codes(codes,cats)
def summarize(df,labels):
 work=df.copy();work["label"]=pd.Series(labels,index=df.index).astype("string").fillna("MISSING_LABEL")
 work["native_p_lt_0_05"]=work["pval"]<0.05
 g=work.groupby("label",observed=True)
 out=g.agg(cell_number=("norm_score","size"),mean_norm_score=("norm_score","mean"),median_norm_score=("norm_score","median"),mean_zscore=("zscore","mean"),native_p_lt_0_05_count=("native_p_lt_0_05","sum"))
 out["native_p_lt_0_05_fraction"]=out["native_p_lt_0_05_count"]/out["cell_number"]
 out["fraction_of_all_cells"]=out["cell_number"]/len(df)
 return out
def test():
 ids=pd.Index(["a","b","c"]);df=pd.DataFrame({k:[0.1,0.2,0.3] for k in ["raw_score","norm_score","mc_pval","pval","nlog10_pval","zscore"]},index=ids)
 validate(df,ids)
 s=summarize(df,["L1a","L1b","L1a"]);assert s.cell_number.sum()==3 and s.loc["L1a","cell_number"]==2
 for kind,x,y in [("duplicate",df.set_axis(["a","a","c"]),ids),("order",df,ids[::-1]),("bad_p",df.assign(pval=2),ids),("nan",df.assign(norm_score=np.nan),ids)]:
  try:validate(x,y)
  except ValueError:pass
  else:raise AssertionError(kind)
 return ["valid_schema_and_all_cell_denominator","group_counts_preserved","reject_duplicate_ID","reject_reordered_ID","reject_bad_P","reject_nonfinite"]
def main():
 cfg=json.loads((W/"config.json").read_text());receipt=json.loads((W/"RECEIPT.json").read_text())
 for item in receipt["bound_files"]:
  if sha(item["path"])!=item["sha256"]:raise ValueError("Bound file drift")
 sys.path.insert(0,cfg["coordinator_code_dir"]);import coordinator as coord
 coord.verify()
 cpu=coord.cpu_budget(coord.C,{os.getpid(),os.getppid()})
 if cpu["actual_other_threads"]+1>40:raise RuntimeError("Project CPU budget unavailable")
 if coord.mount_budget(coord.C)["reasons"]:raise RuntimeError("Disk gate failed")
 out=Path(cfg["output"]);out.mkdir(exist_ok=False)
 def status(state,**kw):(W/"STATUS.json").write_text(json.dumps(dict(status=state,pid=os.getpid(),time=time.time(),**kw),indent=2))
 status("VALIDATING_COMPLETE_ORDINARY_SCORE")
 source=Path(cfg["score"])
 before=(source.stat().st_ino,source.stat().st_size,source.stat().st_mtime_ns);digest=sha(source)
 # Never read full_score: the native producer is still writing it.
 for fd in (Path("/proc")/str(cfg["native_pid"])/"fd").glob("*"):
  try:target=os.readlink(fd)
  except FileNotFoundError:continue
  if target==str(source):raise RuntimeError("Ordinary score is still open by producer")
 df=pd.read_csv(source,sep="\t",index_col=0) # Full gzip traversal verifies its end/CRC.
 import h5py
 with h5py.File(cfg["atlas"],"r") as h:
  g=h["obs"];key=g.attrs["_index"];ids=pd.Index(read_col(g[key]))
  validate(df,ids)
  labels={k:read_col(g[k]) for k in cfg["label_keys"]}
 if before!=(source.stat().st_ino,source.stat().st_size,source.stat().st_mtime_ns) or digest!=sha(source):raise ValueError("Score changed during read")
 status("DESCRIPTIVE_SUMMARY",n_cells=len(df))
 import matplotlib;matplotlib.use("Agg")
 import matplotlib.pyplot as plt
 for key,lab in labels.items():
  summary=summarize(df,lab);summary.to_csv(out/(key+"_native_score_summary.tsv"),sep="\t")
  x=summary.sort_values("mean_norm_score")
  fig,ax=plt.subplots(figsize=(10,max(4,len(x)*0.2)))
  ax.barh(x.index,x.mean_norm_score,color="#497ba8");ax.set_xlabel("Mean native scDRS normalized score")
  ax.set_title("Primary r0.1 | descriptive only; native group tests pending")
  fig.tight_layout();fig.savefig(out/(key+"_mean_score.png"),dpi=160);fig.savefig(out/(key+"_mean_score.pdf"));plt.close(fig)
 qc={"status":"PASS_ORDINARY_SCORE_ONLY","n_cells":len(df),"source":str(source),"source_sha256":digest,"native_group_tests_completed":False,"full_score_read":False,"atlas_sha256":cfg["atlas_sha256"],"quantiles":df[["norm_score","zscore","pval"]].quantile([0,.01,.25,.5,.75,.99,1]).to_dict()}
 (out/"ORDINARY_SCORE_QC.json").write_text(json.dumps(qc,indent=2))
 (out/"README_CN.md").write_text("首场景已完成普通score文件的描述性汇总；不读取仍在写入的full_score。全部细胞与原图谱顺序一一核对。各群均值/中位数与原生单细胞P<0.05比例仅用于描述，不能替代原生组检验、组层面FDR或供者层面证据。Unknown/Not_assessed不剔除；不使用临床数据定义组，不修改标签。三场景完整scDRS尚未完成。\n")
 manifest={p.name:sha(p) for p in out.iterdir() if p.is_file()}
 (out/"MANIFEST.json").write_text(json.dumps(manifest,indent=2))
 status("COMPLETE_DESCRIPTIVE_ONLY",n_cells=len(df),output=str(out))
if __name__=="__main__":
 if "--test" in sys.argv:
  print(json.dumps({"status":"PASS","tests":test()}))
 else:
  try:main()
  except Exception:
   import traceback
   (W/"FAILED.json").write_text(json.dumps({"error":traceback.format_exc(),"time":time.time()}));raise
