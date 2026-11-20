from common import *
while not (W/'communication/ALL_VALID_RESULTS_INDEX.json').exists():
 dump(W/'status/communication_disease.json',dict(status='WAITING_FOR_VALIDATED_229_SAMPLE_L1_COMMUNICATION',updated=time.time()));time.sleep(45)
run_job('communication_disease',[PY,W/'code/communication_disease_calculate.py'],W/'logs/communication_disease.log',reserve=3)
