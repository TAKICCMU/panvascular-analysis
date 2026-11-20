"""Only COMPLETE unlocks downstream dependencies; process errors always dominate."""
from enum import Enum
class Status(str,Enum):
    PENDING='PENDING'
    RUNNING='RUNNING'
    COMPLETE='COMPLETE'
    PREPARED_NOT_EXECUTED='PREPARED_NOT_EXECUTED'
    PARTIAL='PARTIAL'
    BLOCKED='BLOCKED'
    BLOCKED_DEPENDENCY='BLOCKED_DEPENDENCY'
    FAILED='FAILED'

def receipt_state(returncode, receipt):
    if returncode!=0:return 'BLOCKED' if returncode==3 else 'FAILED'
    if not receipt or receipt.get('status') not in {s.value for s in Status}:return 'FAILED'
    state=receipt['status']
    return state if state in {'COMPLETE','PREPARED_NOT_EXECUTED','PARTIAL','BLOCKED','FAILED'} else 'FAILED'
