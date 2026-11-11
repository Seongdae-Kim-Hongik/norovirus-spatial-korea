suppressPackageStartupMessages({library(INLA);library(dplyr);library(spdep);library(glmmTMB);library(car)})
INLA::inla.setOption(num.threads='4:1')
B<-normalizePath(Sys.getenv("R2_CODE_ROOT","."));R<-B;O<-file.path(B,"output");dir.create(O,showWarnings=FALSE);dir.create(file.path(O,"fits"),showWarnings=FALSE)
source(file.path(R,'r1_covariate_dictionary.R'))
x<-readRDS(file.path(B,'input/analysis_input.rds'));d<-as.data.frame(x$cor_merged);d<-d[d$region %in% x$shp_main$region & d$population>0,];d<-d[order(d$region,d$year),];g_main<-x$g_main
stopifnot(nrow(d)==1115,sum(d$cases)==4801)
a<-read.csv(file.path(O,'corrected_moisture_panel.csv'),check.names=FALSE);d<-left_join(d,a[,setdiff(names(a),c('cases','population'))],by=c('region','year'))
phase<-'cases'
ob<-read.csv(file.path(B,'input/outbreaks_R1_verified.csv'));d$outbreaks<-ob$outbreaks[match(paste(d$region,d$year),paste(ob$region,ob$year))];stopifnot(sum(d$outbreaks)==220)
stopifnot(nrow(d)==1115,sum(d$cases)==4801)
write.csv(d[,c('region','year','cases','outbreaks','population')],file.path(O,paste0('frozen_outcome_panel_',phase,'.csv')),row.names=FALSE)
orig<-read.csv(file.path(B,'input/principal_determinants_M4.csv'))
dict<-rbind(TV_v6,TV_RURAL);dict<-dict[!duplicated(dict$eng),];dict$safe<-paste0(dict$eng,'_z')
terms_by_group<-split(orig$term,orig$model);pc_bym<-list(prec.unstruct=list(prior='pc.prec',param=c(.5,.01)),prec.spatial=list(prior='pc.prec',param=c(.5,.01)));pc_prec<-list(prec=list(prior='pc.prec',param=c(.5,.01)))
randoms<-c(M1='',M2="+ f(idarea,model='besag',graph=g_main,scale.model=TRUE,hyper=pc_prec)",M3="+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)",M4="+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+f(idarea_time,model='iid',hyper=pc_prec)",M5="+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+f(idtime,model='rw1',hyper=pc_prec)",M6="+ f(idarea,model='bym',graph=g_main,scale.model=TRUE,hyper=pc_bym)+f(idtime,model='rw1',hyper=pc_prec)+f(idarea_time,model='iid',hyper=pc_prec)")
scale_records<-list()
make_data<-function(group='Total',variant='corrected'){
 z<-if(group=='Total')d else d[d$area_2==ifelse(group=='Urban','도시','농촌'),];ts<-terms_by_group[[group]]
 for(tm in ts){cd<-dict$code[match(tm,dict$safe)];v<-as.numeric(z[[cd]]);form<-x$FIX_FORMS[[cd]];if(grepl('사육두수',cd))v<-v/10000 else if(!is.null(form)&&form=='log1p')v<-log1p(pmax(v,0));z[[tm]]<-(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)}
 if(variant!='legacy'){
  val<-if(variant=='weighted')z$moisture_weighted else if(variant=='positive_only')z$positive_only_mean else z$moisture_mean
  if(variant=='exclude_below30')val[val<30]<-NA_real_
  av<-is.finite(val);mu<-mean(val[av]);sg<-sd(val[av]);z$sludge_moisture_z<-ifelse(av,(val-mu)/sg,0);z$moisture_observed<-as.integer(av);ts<-c(ts,'moisture_observed')
  scale_records[[paste(group,variant,sep='_')]]<<-data.frame(group=group,variant=variant,N=nrow(z),observed=sum(av),mean=mu,SD=sg)
  if(variant=='complete_case')z<-z[av,]
 }
 z$idarea<-match(z$region,x$shp_main$region);z$idtime<-match(z$year,2020:2024);z$idarea_time<-seq_len(nrow(z));z$year_f<-factor(z$year)
 if(variant=='complete_case')ts<-setdiff(ts,'moisture_observed')
 stopifnot(all(complete.cases(z[,c('cases','population',ts)])))
 list(data=z,terms=ts)
}
safe_inla<-function(...){
 for(attempt in 1:5){
  fit<-tryCatch(INLA::inla(...,num.threads='4:1'),error=function(e){cat('INLA_ERROR attempt',attempt,conditionMessage(e),'\n');NULL})
  if(!is.null(fit) && !is.null(fit$dic$dic) && is.finite(fit$dic$dic) && fit$dic$dic>0 && all(is.finite(fit$summary.fixed$mean)) && all(fit$summary.fixed$sd>0))return(fit)
  cat('REJECT_INVALID_NUMERICAL_FIT attempt',attempt,'DIC',if(is.null(fit))NA else fit$dic$dic,'\n')
 }
 stop('INLA failed all five attempts; do not report this fit')
}
rows<-list();diagnostics<-list()
runfit<-function(dd,tag,model='M4',outcome='cases',repeat_id=1){
 dat<-dd$data;ts<-dd$terms;key<-paste(tag,model,outcome,repeat_id,sep='__');cache<-file.path(O,'fits',paste0(key,'.rds'))
 fm<-as.formula(paste(outcome,'~',paste(ts,collapse='+'),'+offset(log(population+1))',randoms[[model]]))
 cat(format(Sys.time()),'START',key,'N',nrow(dat),'\n');flush.console()
 if(file.exists(cache)){cf<-readRDS(cache);if(!is.finite(cf$dic$dic)||cf$dic$dic<=0){file.rename(cache,paste0(cache,'.rejected'));cat('REJECT_CACHED_INVALID',key,'\n')}}
 if(file.exists(cache))fit<-readRDS(cache) else {
  fit<-safe_inla(fm,family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE)
  if(model!='M1'){
   ref<-safe_inla(as.formula(paste(outcome,'~',paste(ts,collapse='+'),'+offset(log(population+1))')),family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE),verbose=FALSE)
   tries<-0;while(fit$dic$dic>ref$dic$dic+50 && tries<5){tries<-tries+1;cat('DEGENERACY REFIT',key,tries,'\n');fit<-safe_inla(fm,family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE)}
   stopifnot(fit$dic$dic<=ref$dic$dic+50)
  }
  saveRDS(fit,cache)
 }
 s<-fit$summary.fixed;rr<-data.frame(tag=tag,model=model,outcome=outcome,repeat_id=repeat_id,term=rownames(s),N=nrow(dat),cases=sum(dat[[outcome]]),IRR=exp(s$mean),lo=exp(s[['0.025quant']]),hi=exp(s[['0.975quant']]),DIC=fit$dic$dic,WAIC=fit$waic$waic)
 rows[[key]]<<-rr;write.csv(do.call(rbind,rows),file.path(O,'all_model_estimates.csv'),row.names=FALSE)
 vifmax<-tryCatch(max(car::vif(lm(as.formula(paste(outcome,'~',paste(ts,collapse='+'))),data=dat))),error=function(e)NA_real_)
 diagnostics[[key]]<<-data.frame(key=key,N=nrow(dat),DIC=fit$dic$dic,WAIC=fit$waic$waic,max_VIF=vifmax,CPO_failure=sum(fit$cpo$failure>0,na.rm=TRUE));write.csv(do.call(rbind,diagnostics),file.path(O,'model_diagnostics.csv'),row.names=FALSE)
 cat('DONE',key,'DIC',fit$dic$dic,'\n');flush.console();fit
}
# Preserve R1 covariate sets and transformations; only the audited moisture exposure changes.
# Availability is an administrative indicator, not an indicator of treatment coverage.
# Its centered zero design entry denotes an unavailable measurement, never 0% physical moisture.
runfit(make_data('Total','legacy'),'Total_legacy')
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group,'corrected');write.csv(dd$data[,c('region','year','cases','outbreaks','population',dd$terms)],file.path(O,paste0('design_',group,'.csv')),row.names=FALSE)
 for(rep in 1:2)runfit(dd,paste0(group,'_corrected'),'M4','cases',rep)
}
for(m in c('M1','M2','M3','M5','M6'))runfit(make_data(),'Total_corrected',m)
runfit(make_data('Total','legacy'),'Total_legacy')
for(v in c('weighted','positive_only','exclude_below30','complete_case'))runfit(make_data('Total',v),paste0('Total_',v))
if(phase=='outbreaks')for(group in c('Total','Urban','Rural'))runfit(make_data(group),paste0(group,'_corrected'),'M4','outbreaks')
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group);dd$data<-dd$data[!(dd$data$region=='전라북도남원시'&dd$data$year==2024),];runfit(dd,paste0(group,'_drop_Namwon2024'))
}
# Paediatric battery retains R1 anchors for comparability; does not re-select on R2 significance.
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group);z<-dd$data;ped<-c(child_0_4_z='영유아비율_0_4',child_5_9_z='아동비율_5_9',ped_gastro_adm_z='예방가능입원율_소아위장관염',daycare_per1k_z='유아천명당보육시설수')
 for(tm in names(ped)){v<-as.numeric(z[[ped[[tm]]]]);form<-x$FIX_FORMS[[ped[[tm]]]];if(!is.null(form)&&form=='log1p')v<-log1p(pmax(v,0));z[[tm]]<-as.numeric(scale(v))}
 # All R1 anchors, reconstructed in each stratum even when absent from that stratum's principal set.
 for(tm in c('livestock_total_z','fiscal_auto_z'))if(!tm%in%names(z)){cd<-dict$code[match(tm,dict$safe)];v<-as.numeric(z[[cd]]);if(grepl('사육두수',cd))v<-v/10000;z[[tm]]<-as.numeric(scale(v))}
 ts<-c('livestock_total_z','sludge_moisture_z','moisture_observed','fiscal_auto_z',names(ped));z<-z[complete.cases(z[,ts]),]
 # VIF pruning is based only on predictors, with deterministic removal of the largest VIF.
 repeat{vv<-car::vif(lm(as.formula(paste('cases~',paste(ts,collapse='+'))),data=z));if(max(vv)<5||length(ts)<3)break;drop<-names(which.max(vv));cat('PAED VIF REMOVE',group,drop,max(vv),'\n');ts<-setdiff(ts,drop)}
 runfit(list(data=z,terms=ts),paste0(group,'_paediatric'))
}
write.csv(do.call(rbind,scale_records),file.path(O,'moisture_scaling.csv'),row.names=FALSE)
# Independent likelihood check with the same fixed effects, no spatial random effects (R1 S11 convention).
mles<-list()
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group);fm<-as.formula(paste('cases~',paste(dd$terms,collapse='+'),'+offset(log(population+1))'));fit<-glmmTMB(fm,family=nbinom2,data=dd$data);s<-summary(fit)$coefficients$cond
 mles[[group]]<-data.frame(group=group,term=rownames(s),IRR=exp(s[,1]),lo=exp(s[,1]-1.96*s[,2]),hi=exp(s[,1]+1.96*s[,2]),pdHess=fit$sdr$pdHess);write.csv(do.call(rbind,mles),file.path(O,'mle_crosscheck.csv'),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'));cat('ALL REQUESTED FITS COMPLETED\n')
