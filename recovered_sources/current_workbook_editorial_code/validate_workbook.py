from pathlib import Path
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter as L
import csv,math,json
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921';B=O/'轻量备份';S=O/'待替换';src=P/'Supplementary_Tables.xlsx'
conversions=list(csv.reader((O/'numeric_text_conversions.tsv').open(),delimiter='\t'))[1:];text_changes=list(csv.reader((O/'label_text_changes.tsv').open(),delimiter='\t'))[1:]
s=(O/'code/format_workbook.py').read_text();exec(s[s.index('a=load_workbook(B/src.name,read_only=True)'):])
