# Pool run 1 (offset_no_plus1_M4.R) and run 2 (public sensitivity_offset_no_plus1.R test) with the saved principal fits.
suppressPackageStartupMessages(library(INLA))
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4)
O0<-normalizePath(args[1]);O1<-normalizePath(args[2]);O2<-normalizePath(args[3]);dir.create(args[4],recursive=TRUE,showWarnings=FALSE);setwd(args[4])
fx<-function(p){f<-readRDS(p);s<-f$summary.fixed;data.frame(term=rownames(s),IRR=exp(s$mean),lo=exp(s[['0.025quant']]),hi=exp(s[['0.975quant']]),DIC=f$dic$dic,WAIC=f$waic$waic)}
cred<-function(t)paste(sort(t$term[t$term!='(Intercept)'&(t$lo>1|t$hi<1)]),collapse=';')
ld<-function(a,b){m<-merge(a,b,by='term');m<-m[m$term!='(Intercept)',];max(abs(log(c(m$IRR.x,m$lo.x,m$hi.x))-log(c(m$IRR.y,m$lo.y,m$hi.y))))}
ad<-function(a,b){m<-merge(a,b,by='term');m<-m[m$term!='(Intercept)',];max(abs(c(m$IRR.x-m$IRR.y,m$lo.x-m$lo.y,m$hi.x-m$hi.y)))}
rows<-list();sm<-list()
for(g in c('Total','Urban','Rural')){
 P<-c(p1_saved1=file.path(O0,paste0(g,'_corrected__M4__cases__1.rds')),p1_saved2=file.path(O0,paste0(g,'_corrected__M4__cases__2.rds')),
      p1_run1=file.path(O1,paste0(g,'_corrected_p1fresh__M4__cases__1.rds')),p1_run2=file.path(O2,paste0(g,'_corrected_p1fresh__M4__cases__1.rds')),
      nop1_run1a=file.path(O1,paste0(g,'_corrected_nop1__M4__cases__1.rds')),nop1_run1b=file.path(O1,paste0(g,'_corrected_nop1__M4__cases__2.rds')),
      nop1_run2a=file.path(O2,paste0(g,'_corrected_nop1__M4__cases__1.rds')),nop1_run2b=file.path(O2,paste0(g,'_corrected_nop1__M4__cases__2.rds')))
 T<-lapply(P,fx);p1<-grep('^p1',names(T),value=TRUE);np<-grep('^nop1',names(T),value=TRUE)
 for(k in names(T)){t<-T[[k]];rows[[paste(g,k)]]<-data.frame(group=g,fit=k,credible=cred(t),DIC=t$DIC[1],WAIC=t$WAIC[1])}
 pr<-function(A,B)unlist(lapply(A,function(a)lapply(B,function(b)if(a!=b)c(ld(T[[a]],T[[b]]),ad(T[[a]],T[[b]])) else NULL)),recursive=FALSE)
 off<-do.call(rbind,pr(p1,np));rp<-do.call(rbind,c(pr(p1,p1),pr(np,np)))
 sm[[g]]<-data.frame(group=g,max_log_offset=max(off[,1]),max_log_repeat=max(rp[,1]),max_abs_offset=max(off[,2]),max_abs_repeat=max(rp[,2]),
   DIC_p1=paste(range(round(sapply(T[p1],function(t)t$DIC[1]),1)),collapse='–'),DIC_nop1=paste(range(round(sapply(T[np],function(t)t$DIC[1]),1)),collapse='–'),
   WAIC_p1=paste(range(round(sapply(T[p1],function(t)t$WAIC[1]),1)),collapse='–'),WAIC_nop1=paste(range(round(sapply(T[np],function(t)t$WAIC[1]),1)),collapse='–'))
 if(g=='Rural')print(do.call(rbind,lapply(names(T),function(k){t<-T[[k]];r<-t[t$term=='school_gw_use_z',];data.frame(fit=k,IRR=r$IRR,lo=r$lo,hi=r$hi)})),digits=6)
 if(g=='Total')print(do.call(rbind,lapply(names(T),function(k){t<-T[[k]];r<-t[t$term%in%c('livestock_total_z','sludge_moisture_z'),];data.frame(fit=k,term=r$term,IRR=round(r$IRR,3),lo=round(r$lo,3),hi=round(r$hi,3))})))
}
A<-do.call(rbind,rows);S<-do.call(rbind,sm);write.csv(A,'pooled_fits_credible.csv',row.names=FALSE);write.csv(S,'pooled_summary.csv',row.names=FALSE);print(A);print(S,digits=4)
