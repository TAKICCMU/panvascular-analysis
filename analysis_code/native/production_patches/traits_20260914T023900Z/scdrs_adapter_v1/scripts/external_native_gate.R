# Native workers are reachable only after the common Python release preflight.
keys<-c('DOWNSTREAM_AUTH_CONFIG','DOWNSTREAM_AUTH_RECEIPT','DOWNSTREAM_AUTH_PYTHON','DOWNSTREAM_AUTH_GATE','DOWNSTREAM_AUTH_THREADS')
v<-Sys.getenv(keys,unset='')
if(any(!nzchar(v)))stop('Native external worker requires authenticated common release preflight')
threads<-as.integer(v[['DOWNSTREAM_AUTH_THREADS']]);if(!is.finite(threads)||threads<1||threads>20)stop('CPU worker cap20')
for(k in c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'))do.call(Sys.setenv,setNames(list(as.character(threads)),k))
status<-system2(v[['DOWNSTREAM_AUTH_PYTHON']],args=c(shQuote(v[['DOWNSTREAM_AUTH_GATE']]),'--config',shQuote(v[['DOWNSTREAM_AUTH_CONFIG']]),'--check-receipt',shQuote(v[['DOWNSTREAM_AUTH_RECEIPT']])))
if(status!=0L)stop('Common authorization/preflight receipt validation failed')
