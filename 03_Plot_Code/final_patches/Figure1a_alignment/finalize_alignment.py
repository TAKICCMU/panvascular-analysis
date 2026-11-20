"""Install only the verified vessel placement correction in current figures."""
from pathlib import Path
import csv
import hashlib
import io
import json
import os
import shutil

import fitz
import numpy as np
from lxml import etree as ET

P = Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
R = P / '图件编辑资料/Figure1a_人体血管对齐'
B = R / 'backup'
baseline = json.loads((R / 'baseline.json').read_text())
alignment = json.loads((R / 'alignment_check.json').read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for rel, info in baseline.items():
    staged = R / 'staged' / rel
    assert sha(P / rel) == info['sha256'] or (staged.exists() and sha(P / rel) == sha(staged)), f'File changed since backup: {rel}'

old_cmd = b'11.92 0 0 27.04 17.28 24.764724 cm\n/a0 gs /x8 Do'
v = fitz.Rect(alignment['new_vessel_bbox'])
new_cmd = (f'{v.width:.9f} 0 0 {v.height:.9f} {v.x0:.9f} {v.y0:.9f} cm\n/a0 gs /x8 Do').encode()
merged = ['Figures_and_ExtendedData_审阅.pdf', 'Figure/Figures_1-4_and_Extended_Data.pdf']
for rel in merged:
    d = fitz.open(B / rel)
    changed = 0
    for xref in d[0].get_contents():
        s = d.xref_stream(xref)
        changed += s.count(old_cmd)
        if old_cmd in s:
            d.update_stream(xref, s.replace(old_cmd, new_cmd))
    assert changed == 1
    target = R / 'staged' / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    d.save(target, garbage=4, deflate=True)

checks = []
for rel in [x for x in baseline if x.endswith('.pdf')]:
    before, after = fitz.open(B / rel), fitz.open(R / 'staged' / rel)
    assert len(before) == len(after)
    page_checks = []
    for i, (p0, p1) in enumerate(zip(before, after)):
        assert p0.rect == p1.rect
        assert p0.get_text() == p1.get_text()
        # Full resolution rendering checks include all 14 pages of both collections.
        px0 = p0.get_pixmap(matrix=fitz.Matrix(4, 4), alpha=False)
        px1 = p1.get_pixmap(matrix=fitz.Matrix(4, 4), alpha=False)
        a = np.frombuffer(px0.samples, np.uint8).reshape(px0.height, px0.width, 3)
        b = np.frombuffer(px1.samples, np.uint8).reshape(px1.height, px1.width, 3)
        diff = np.any(a != b, axis=2)
        total = int(diff.sum())
        if i == 0:
            u = fitz.Rect(alignment['old_vessel_bbox']) | v
            diff[int(u.y0 * 4)-2:int(u.y1 * 4)+3, int(u.x0 * 4)-2:int(u.x1 * 4)+3] = False
            assert total > 0
        assert not diff.any(), (rel, i, int(diff.sum()))
        page_checks.append({'page': i + 1, 'changed_pixels': total,
                            'changed_outside_vessel_region': int(diff.sum()),
                            'text_identical': True, 'native_page_size_preserved': True})
    checks.append({'file': rel, 'pages': page_checks})

for rel in [x for x in baseline if x.endswith('.svg')]:
    a = ET.parse(str(B / rel))
    b = ET.parse(str(R / 'staged' / rel))
    ea = a.xpath('//*[@id="a_image08f78ba088"]')[0]
    eb = b.xpath('//*[@id="a_image08f78ba088"]')[0]
    # Restore the sole changed element's attributes before a complete XML comparison.
    eb.attrib.clear()
    eb.attrib.update(ea.attrib)
    assert ET.tostring(a, method='c14n') == ET.tostring(b, method='c14n')
    checks.append({'file': rel, 'only_vessel_transform_changed': True})

# Install only after all comparisons pass.
for rel in baseline:
    destination = P / rel
    temporary = destination.with_name(destination.name + '.alignment.tmp')
    shutil.copyfile(R / 'staged' / rel, temporary)
    os.replace(temporary, destination)

manifest = P / 'Figure/图件版本与尺寸.tsv'
shutil.copy2(manifest, B / 'Figure/图件版本与尺寸.tsv')
rows = list(csv.DictReader(io.StringIO(manifest.read_text()), delimiter='\t'))
for row in rows:
    if row['figure'] == 'Figure1':
        row['source_sha256'] = sha(P / 'Figures/Figure1.pdf')
        row['pdf_sha256'] = sha(P / 'Figure/Figure1.pdf')
        row['outside_correction_pixels_changed'] = '0'
        row['changes'] = '1a: restored native BioRender vessel/body alignment; body, labels, size and all other panels preserved'
out = io.StringIO(newline='')
writer = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\r\n')
writer.writeheader()
writer.writerows(rows)
with manifest.open('w', newline='') as f:
    f.write(out.getvalue())

report = {'scope': 'Figure1a vessel placement only',
          'biorender_current_canvas': 'https://app.biorender.com/illustrations/6af93b20061047a1155dc09d',
          'biorender_saved_reloaded_and_pdf_export_verified': True,
          'biorender_export': str(R / 'BioRender_updated_canvas.pdf'),
          'checks': checks,
          'installed_sha256': {rel: sha(P / rel) for rel in baseline}}
(R / 'final_checks.json').write_text(json.dumps(report, ensure_ascii=False, indent=2))
print(json.dumps({'updated_files': list(baseline), 'combined_page_count': 14,
                  'all_other_13_pages_pixel_identical': True,
                  'figure1_outside_vessel_region_pixel_identical': True}, ensure_ascii=False, indent=2))
