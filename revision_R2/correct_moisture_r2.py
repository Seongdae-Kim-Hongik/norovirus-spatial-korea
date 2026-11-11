"""R2 source correction. No modification of R1 input or outcomes.
Facility-row arithmetic mean of positive moisture values for positive sludge production.
Subtotal/total rows and transferred-only zero-production rows do not enter the mean.
A generation-weighted mean and positive-value-only mean are provided as sensitivities.
Missing physical measurements remain NA; no missing/structural value becomes 0%.
"""
from pathlib import Path
import pandas as pd,numpy as np,openpyxl,json,hashlib
R=Path(__file__).resolve().parent;S=R/'input';O=R/'output';O.mkdir(exist_ok=True)
f=pd.read_csv(S/'R1_고정패널_검증용.csv');f=f[f.region.isin(pd.read_csv(S/'R1_분석지역223.csv').region)].copy();assert len(f)==1115 and f.cases.sum()==4801
prov={'서울특별시':'서울시','부산광역시':'부산시','대구광역시':'대구시','인천광역시':'인천시','광주광역시':'광주시','대전광역시':'대전시','울산광역시':'울산시','세종특별자치시':'세종시','강원특별자치도':'강원도','전북특별자치도':'전라북도','제주특별자치도':'제주도'}
allowed=set(prov)|{'경기도','강원도','충청북도','충청남도','전라북도','전라남도','경상북도','경상남도'};rows=[];excluded=[];sources=[]
for y in range(2020,2025):
 path=S/f'공식_{y}_하수도통계.xlsx';sources.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()));wb=openpyxl.load_workbook(path,read_only=True,data_only=True);ws=wb[next(k for k in wb.sheetnames if '하수찌꺼기 발생' in k)]
 for i,r in enumerate(ws.iter_rows(values_only=True),1):
  if len(r)<10 or r[0] not in allowed:continue
  if r[2] in ['합계','소계'] or not r[3] or str(r[3]).isdigit() or str(r[3]).endswith('개소'):
   excluded.append(dict(year=y,row=i,province=r[0],district=r[1],label=r[2],plant=r[3]));continue
  if not r[1] and r[0]!='세종특별자치시':continue
  reg=(prov.get(r[0],r[0])+str(r[1] or '세종시')).replace('인천시미추홀구','인천시남구').replace('경상북도군위군','대구시군위군')
  rows.append(dict(year=y,region=reg,source_row=i,plant=r[3],generated=r[4],own_generated=r[5],external_received=r[6],transfer_to=r[7],transfer_volume=r[8],moisture=r[9]))
 wb.close()
p=pd.DataFrame(rows)
for c in ['generated','own_generated','external_received','transfer_volume','moisture']:p[c]=pd.to_numeric(p[c],errors='coerce')
p['valid_moisture']=p.moisture.gt(0)&p.moisture.le(100);p['include_primary']=p.valid_moisture&p.generated.gt(0)
p['duplicate_plant_name']=p.duplicated(['region','year','plant'],False)
records=[]
for (reg,y),g in p.groupby(['region','year']):
 v=g[g.include_primary]; records.append(dict(region=reg,year=y,plant_n=len(g),positive_generation_n=int(g.generated.gt(0).sum()),valid_n=len(v),generated=g.generated.sum(min_count=1),moisture_mean=v.moisture.mean(),moisture_weighted=np.average(v.moisture,weights=v.generated) if len(v) else np.nan,positive_only_mean=g.loc[g.valid_moisture,'moisture'].mean(),zero_n=int(g.moisture.eq(0).sum()),missing_n=int(g.moisture.isna().sum()),missing_positive_generation_n=int((g.generated.gt(0)&~g.valid_moisture).sum()),coverage_mass=v.generated.sum()/g.generated.sum() if g.generated.sum()>0 else np.nan))
a=pd.DataFrame(records);panel=f[['region','year','cases','population','함수율(%,탈수기준)']].merge(a,how='left',on=['region','year'],validate='one_to_one');panel=panel.rename(columns={'함수율(%,탈수기준)':'legacy_moisture'})
panel['measurement_available']=panel.moisture_mean.notna().astype(int)
panel['status']=np.select([panel.plant_n.isna(),panel.positive_generation_n.eq(0),panel.moisture_mean.isna(),panel.coverage_mass.lt(.999999)],['no_facility_row','no_positive_sludge_generation','positive_generation_moisture_unavailable','partial_measured_generation'],default='measured_generation_complete')
assert len(panel)==1115 and panel.cases.sum()==4801
panel.to_csv(O/'corrected_moisture_panel.csv',index=False);p.to_csv(O/'facility_rows_with_inclusion.csv',index=False);pd.DataFrame(excluded).to_csv(O/'excluded_total_subtotal_rows.csv',index=False)
summary={'district_years':len(panel),'cases':int(panel.cases.sum()),'positive_measurement_n':int(panel.measurement_available.sum()),'status':panel.status.value_counts().to_dict(),'physical_min':panel.moisture_mean.min(),'physical_max':panel.moisture_mean.max(),'physical_mean':panel.moisture_mean.mean(),'physical_sd':panel.moisture_mean.std(),'positive_below30':int(panel.moisture_mean.lt(30).sum()),'duplicate_name_rows':int(p.duplicate_plant_name.sum()),'input_sources':sources}
(O/'source_correction_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2))
