# Offset sensitivity for the principal M4 (added in v2.4.1, 2026-11-21).
# The published models use offset(log(population + 1)), a zero-population guard inherited from R1.
# All analysed district-year populations are positive, so this script refits the principal M4 in each
# stratum with offset(log(population)) under otherwise identical settings (two runs), refits the
# published +1 specification once more in the same session, and compares them with the saved
# principal fits written by reanalyse_corrected_moisture.R (run that script first).
# Differences between saved fits include run-to-run numerical variation and do not isolate the effect of +1.
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
 fm<-as.formula(paste(outcome,'~',paste(ts,collapse='+'),OFF,randoms[[model]]))
 cat(format(Sys.time()),'START',key,'N',nrow(dat),'\n');flush.console()
 if(file.exists(cache)){cf<-readRDS(cache);if(!is.finite(cf$dic$dic)||cf$dic$dic<=0){file.rename(cache,paste0(cache,'.rejected'));cat('REJECT_CACHED_INVALID',key,'\n')}}
 if(file.exists(cache))fit<-readRDS(cache) else {
  fit<-safe_inla(fm,family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE)
  if(model!='M1'){
   ref<-safe_inla(as.formula(paste(outcome,'~',paste(ts,collapse='+'),OFF)),family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE),verbose=FALSE)
   tries<-0;while(fit$dic$dic>ref$dic$dic+50 && tries<5){tries<-tries+1;cat('DEGENERACY REFIT',key,tries,'\n');fit<-safe_inla(fm,family='nbinomial',data=dat,control.inla=list(control.vb=list(enable=FALSE)),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),control.predictor=list(link=1),verbose=FALSE)}
   stopifnot(fit$dic$dic<=ref$dic$dic+50)
  }
  saveRDS(fit,cache)
 }
 s<-fit$summary.fixed;rr<-data.frame(tag=tag,model=model,outcome=outcome,repeat_id=repeat_id,term=rownames(s),N=nrow(dat),cases=sum(dat[[outcome]]),IRR=exp(s$mean),lo=exp(s[['0.025quant']]),hi=exp(s[['0.975quant']]),DIC=fit$dic$dic,WAIC=fit$waic$waic)
 rows[[key]]<<-rr;write.csv(do.call(rbind,rows),file.path(O,'offset_sensitivity_estimates.csv'),row.names=FALSE)
 vifmax<-tryCatch(max(car::vif(lm(as.formula(paste(outcome,'~',paste(ts,collapse='+'))),data=dat))),error=function(e)NA_real_)
 diagnostics[[key]]<<-data.frame(key=key,N=nrow(dat),DIC=fit$dic$dic,WAIC=fit$waic$waic,max_VIF=vifmax,CPO_failure=sum(fit$cpo$failure>0,na.rm=TRUE));write.csv(do.call(rbind,diagnostics),file.path(O,'offset_sensitivity_diagnostics.csv'),row.names=FALSE)
 cat('DONE',key,'DIC',fit$dic$dic,'\n');flush.console();fit
}

cat('MIN_POPULATION_ALL',min(d$population),'\n')
OFF<-'+offset(log(population))'
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group,'corrected');cat('MIN_POPULATION',group,min(dd$data$population),'N',nrow(dd$data),'\n')
 for(rep in 1:2)runfit(dd,paste0(group,'_corrected_nop1'),'M4','cases',rep)
}
OFF<-'+offset(log(population+1))'
for(group in c('Total','Urban','Rural'))runfit(make_data(group,'corrected'),paste0(group,'_corrected_p1fresh'),'M4','cases',1)
# Comparison with the saved principal fits (repeat_id 1 and 2).
fx<-function(k){f<-readRDS(file.path(O,'fits',paste0(k,'.rds')));s<-f$summary.fixed;data.frame(term=rownames(s),IRR=exp(s$mean),lo=exp(s[['0.025quant']]),hi=exp(s[['0.975quant']]),DIC=f$dic$dic,WAIC=f$waic$waic)}
cred<-function(t)paste(sort(t$term[t$term!='(Intercept)'&(t$lo>1|t$hi<1)]),collapse=';')
ld<-function(a,b){m<-merge(a,b,by='term');m<-m[m$term!='(Intercept)',];max(abs(log(c(m$IRR.x,m$lo.x,m$hi.x))-log(c(m$IRR.y,m$lo.y,m$hi.y))))}
summ<-list()
for(g in c('Total','Urban','Rural')){
 T<-list(p1_rep1=fx(paste0(g,'_corrected__M4__cases__1')),p1_rep2=fx(paste0(g,'_corrected__M4__cases__2')),p1_fresh=fx(paste0(g,'_corrected_p1fresh__M4__cases__1')),
         nop1_rep1=fx(paste0(g,'_corrected_nop1__M4__cases__1')),nop1_rep2=fx(paste0(g,'_corrected_nop1__M4__cases__2')))
 summ[[g]]<-data.frame(group=g,min_population=min(make_data(g,'corrected')$data$population),
  credible_p1_rep1=cred(T$p1_rep1),credible_p1_rep2=cred(T$p1_rep2),credible_p1_fresh=cred(T$p1_fresh),credible_nop1_rep1=cred(T$nop1_rep1),credible_nop1_rep2=cred(T$nop1_rep2),
  max_logdiff_offset=ld(T$p1_rep1,T$nop1_rep1),max_logdiff_repeat_p1=max(ld(T$p1_rep1,T$p1_rep2),ld(T$p1_rep1,T$p1_fresh),ld(T$p1_rep2,T$p1_fresh)),max_logdiff_repeat_nop1=ld(T$nop1_rep1,T$nop1_rep2),
  DIC_p1=T$p1_rep1$DIC[1],DIC_nop1=T$nop1_rep1$DIC[1],WAIC_p1=T$p1_rep1$WAIC[1],WAIC_nop1=T$nop1_rep1$WAIC[1])
}
S<-do.call(rbind,summ);write.csv(S,file.path(O,'offset_sensitivity_summary.csv'),row.names=FALSE);print(S)
cat('OFFSET SENSITIVITY COMPLETED\n')
