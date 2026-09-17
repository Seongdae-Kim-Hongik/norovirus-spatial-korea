# Evidence grading for every covariate entered in the principal models.
# Grading rules were fixed before the evidence fits were inspected:
#   Robust       : credible in the principal model AND credible in the same direction in >=80% of the sensitivity specifications
#                  in which the covariate was entered AND same direction in the maximum-likelihood fit
#   Supported    : credible in the principal model AND same direction (not necessarily credible) in >=80% of sensitivity specifications
#   Suggestive   : not credible in the principal model AND posterior direction probability >=0.95 AND same direction in >=80% of specifications
#   Inconclusive : all other covariates
import pandas as pd, numpy as np, math, os
import sys
# usage: python3 assemble_evidence.py <principal_evidence.csv> <ridge31_evidence.csv> <sensitivity_all.csv> <mle_crosscheck_glmmTMB.csv> <out.csv>
P=pd.read_csv(sys.argv[1])
RG=pd.read_csv(sys.argv[2])
SX=pd.read_csv(sys.argv[3]); SX=SX[SX.spec!="principal"]
X=pd.read_csv(sys.argv[4]).drop_duplicates(["model","term"])
SPECS=sorted(SX.spec.unique())
rows=[]
for r in P.itertuples():
    d = 1 if r.IRR >= 1 else -1
    cred = (r.lo > 1) or (r.hi < 1)
    pdir = r.P_pos if d == 1 else r.P_neg
    s = SX[(SX.model == r.model) & (SX.term == r.term)]
    n_av = len(s); n_same = int(((np.sign(np.log(s.IRR)) == d)).sum()); n_cred = int((((s.lo > 1) & (d == 1)) | ((s.hi < 1) & (d == -1))).sum())
    thr = math.ceil(0.8 * n_av) if n_av else 999
    m = X[(X.model == r.model) & (X.term == r.term)]
    mle_same = (not m.empty) and (np.sign(np.log(m.MLE_IRR.iloc[0])) == d)
    mle_cred = (not m.empty) and ((m.MLE_lo.iloc[0] > 1 and d == 1) or (m.MLE_hi.iloc[0] < 1 and d == -1))
    g = RG[(RG.model == r.model) & (RG.term == r.term)]
    rg_txt = "—" if g.empty else f"{g.IRR.iloc[0]:.2f} ({g.lo.iloc[0]:.2f}–{g.hi.iloc[0]:.2f}){'*' if (g.lo.iloc[0]>1 or g.hi.iloc[0]<1) else ''}"
    rg_pdir = np.nan if g.empty else (g.P_pos.iloc[0] if d == 1 else g.P_neg.iloc[0])
    if cred and n_cred >= thr and mle_same: grade = "Robust"
    elif cred and n_same >= thr: grade = "Supported"
    elif (not cred) and pdir >= 0.95 and n_same >= thr: grade = "Suggestive"
    else: grade = "Inconclusive"
    rows.append(dict(model=r.model, term=r.term, IRR=r.IRR, lo=r.lo, hi=r.hi, credible=cred, direction="+" if d==1 else "−", P_dir=pdir,
                     n_spec=n_av, n_same=n_same, n_cred=n_cred, mle_same=mle_same, mle_cred=mle_cred, ridge31=rg_txt, ridge31_Pdir=rg_pdir, grade=grade))
E=pd.DataFrame(rows)
order={"Robust":0,"Supported":1,"Suggestive":2,"Inconclusive":3}
E["o"]=E.grade.map(order); E=E.sort_values(["model","o","P_dir"],ascending=[True,True,False]).drop(columns="o")
E.to_csv(sys.argv[5],index=False)
pd.set_option("display.width",250)
print(E[E.grade!="Inconclusive"].round(3).to_string(index=False))
print(E.groupby(["model","grade"]).size().unstack(fill_value=0))
print("ridge31 credible:"); print(RG[(RG.lo>1)|(RG.hi<1)][["model","term","IRR","lo","hi"]].round(2).to_string(index=False))
