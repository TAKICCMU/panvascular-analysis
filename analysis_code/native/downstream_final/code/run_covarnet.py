from common import *
env={'COVARNET_WORKERS':'2','COVARNET_RANK_MIN':'2','COVARNET_RANK_MAX':'25','COVARNET_RANK_NRUN':'50','COVARNET_FINAL_NRUN':'200','COVARNET_STABILITY_NRUN':'30','COVARNET_N_BOOTSTRAP':'30','COVARNET_PARENT_NULL_PERMUTATIONS':'2000','COVARNET_SEED':'20260806','COVARNET_PARENT_NULL_SEED':'20260826'}
run_job('covarnet_main_and_refit_sensitivity',[RS,'--vanilla',W/'code/covarnet_nmf.R',W/'covarnet'],W/'logs/covarnet.log',reserve=6,env=env)
run_job('covarnet_parent_null',[RS,'--vanilla',W/'code/covarnet_parent_null.R',W/'covarnet'],W/'logs/covarnet_parent_null.log',reserve=2,env=env)
dump(W/'covarnet/COMPLETE.json',dict(status='COMPLETE_NATIVE_COVARNET_AND_INHERITED_SENSITIVITIES',files={str(p.relative_to(W/'covarnet')):sha(p) for p in (W/'covarnet').rglob('*') if p.is_file()}))
