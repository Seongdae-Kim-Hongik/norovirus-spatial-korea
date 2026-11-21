from pathlib import Path
import pandas as pd,numpy as np,json
import os
p=Path(os.getenv("R2_CODE_ROOT","."))/"output/cpo_diagnostics"
a=pd.concat([pd.read_csv(p/'loo_predictive.csv'),pd.read_csv(p/'loo_extra_predictive.csv')],ignore_index=True);assert len(a)==29 and not a.duplicated(['group','index']).any();a['draws']=20000
q=pd.read_csv(p/'precision_Namwon.csv').iloc[0];i=(a.group=='Total')&(a['index']==925);a.loc[i,'new']=q.cpo;a.loc[i,'mcse']=q.mcse;a.loc[i,'draws']=q.n;a['relative_mcse']=a.mcse/a.new
b=pd.read_csv(p/'gcpo_all_warning_locations.csv');a=a.merge(b[['group','index','gcpo','old_failure']],on=['group','index']);a['gcpo_relative_difference']=(a.gcpo-a.new)/a.new;assert np.isfinite(a.new).all() and (a.new>0).all() and(a.relative_mcse<.05).all();a.to_csv(p/'final_predictive_comparison.csv',index=False)
e=pd.concat([pd.read_csv(p/'loo_effects.csv'),pd.read_csv(p/'loo_extra_effects.csv')],ignore_index=True);assert e[['group','heldout']].drop_duplicates().shape[0]==32
for prefix in ['original','new']:
 e[prefix+'_IRR']=np.exp(e[prefix+'_mean']);e[prefix+'_lower']=np.exp(e[prefix+'_lo']);e[prefix+'_upper']=np.exp(e[prefix+'_hi']);e[prefix+'_credible']=(e[prefix+'_lo']>0)|(e[prefix+'_hi']<0)
e.to_csv(p/'final_effects.csv',index=False)
rows=[]
for (g,tm),v in e[e.heldout>0].groupby(['group','term']):
 rows.append(dict(group=g,term=tm,original_IRR=v.original_IRR.iloc[0],original_lower=v.original_lower.iloc[0],original_upper=v.original_upper.iloc[0],minimum_IRR=v.new_IRR.min(),maximum_IRR=v.new_IRR.max(),credible_original=bool(v.original_credible.iloc[0]),credible_in_deletions=int(v.new_credible.sum()),deletions=len(v)))
s=pd.DataFrame(rows);s.to_csv(p/'single_deletion_summary.csv',index=False)
flips=e[(e.original_credible!=e.new_credible)&(e.term!='(Intercept)')].copy();flips.to_csv(p/'credible_classification_changes.csv',index=False)
print('29 cells, unique district-years',a[['region','year']].drop_duplicates().shape[0]);print('CPO range',a.new.min(),a.new.max(),'maxRSE',a.relative_mcse.max(),'maxGCPOreldiff',a.gcpo_relative_difference.abs().max());print(s[s.term.isin(['livestock_total_z','sludge_moisture_z','school_gw_use_z'])].to_string(index=False));print('COMBINED',e[(e.heldout==0)&e.term.isin(['livestock_total_z','sludge_moisture_z','school_gw_use_z'])][['group','term','new_IRR','new_lower','new_upper']].to_string(index=False));print('FLIPS',flips[['group','heldout','term','new_IRR','new_lower','new_upper']].to_string(index=False))
