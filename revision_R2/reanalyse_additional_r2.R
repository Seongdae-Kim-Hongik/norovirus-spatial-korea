# Additional R2 checks; the outcome panel is byte-for-byte traceable to R1 aggregates.
p<-file.path(Sys.getenv("R2_CODE_ROOT","."),"reanalyse_corrected_moisture.R");s<-paste(readLines(p,warn=FALSE),collapse='\n');eval(parse(text=strsplit(s,'# Preserve R1 covariate sets',fixed=TRUE)[[1]][1]))
old<-read.csv(file.path(O,'all_model_estimates.csv'));rows<-split(old,paste(old$tag,old$model,old$outcome,old$repeat_id,sep='__'));dg<-read.csv(file.path(O,'model_diagnostics.csv'));diagnostics<-split(dg,dg$key)
ob<-read.csv(file.path(B,'input/outbreaks_R1_verified.csv'));d$outbreaks<-ob$outbreaks[match(paste(d$region,d$year),paste(ob$region,ob$year))];stopifnot(sum(d$outbreaks)==220,sum(d$cases)==4801)
# Same posterior model and priors. Start hyperparameter optimisation at a finite
# overdispersion/precision point to avoid numerical failure near default starting values.
safe_inla<-function(formula,...){
 fs<-paste(deparse(formula),collapse=' ');nlatent<-length(gregexpr('f\\(',fs)[[1]]);if(!grepl('f\\(',fs))nlatent<-0L
 nh<-1+nlatent+as.integer(grepl("model = \"bym\"|model='bym'|model = \"bym2\"|model='bym2'",fs))
 for(attempt in 1:5){
  fit<-tryCatch(if(grepl('bym2',fs)) INLA::inla(formula,...,num.threads='4:1') else INLA::inla(formula,...,num.threads='4:1',control.mode=list(theta=c(-3,rep(4,nh-1)),restart=TRUE)),error=function(e){cat('INLA_ERROR',attempt,conditionMessage(e),'\n');NULL})
  if(!is.null(fit)&&is.finite(fit$dic$dic)&&fit$dic$dic>0&&all(is.finite(fit$summary.fixed$mean)))return(fit)
  cat('REJECT_INVALID_NUMERICAL_FIT',attempt,'\n')
 };stop('Failed five numerical attempts')
}
for(group in c('Total','Urban','Rural'))runfit(make_data(group),paste0(group,'_corrected'),'M4','outbreaks')
runfit(make_data('Total'),'Total_corrected','M5',repeat_id=2)
for(v in c('exclude_below30','complete_case'))runfit(make_data('Rural',v),paste0('Rural_',v))
# Source-documented filled years, retaining all other R1 measurements.
filled<-list('자체처리량(톤/년)_건조후처리(2차).3'=2020:2021,'폐수방류량'=2024,'민방위용_개소수'=2023:2024,'민방위용_이용량'=2023:2024,'화장실다녀온후손씻기실천율_조율'=2021:2022,'독거노인비율'=2021,'기초생활수급자수율'=2021:2024,'usual_handwash_rate_adj'=2021:2024,'학교용_이용량'=2023:2024,'사육두수(두)_합계'=2024,'농가수(호)_한육우'=2024,'농가수(호)_합계'=2024,'사육두수(두)_가금'=2024,'사육두수(두)_한육우'=2024)
base_d<-d
for(sens in c('timeinv','linear','omit_filled_covariates','observed2020_2023',paste0('drop_year_',2020:2024),'drop_Anseong2024','drop_both_districts')){
 d<-base_d
 if(sens%in%c('timeinv','linear'))for(v in intersect(names(filled),names(d))){
  values<-as.numeric(d[[v]]);values[d$year%in%filled[[v]]]<-NA_real_;out<-values
  for(rg in unique(d$region)){k<-which(d$region==rg);ok<-is.finite(values[k]);if(!any(ok))next
   if(sens=='timeinv')out[k]<-mean(values[k][ok]) else if(sum(ok)>=2)out[k]<-approx(d$year[k][ok],values[k][ok],xout=d$year[k],rule=2)$y else out[k]<-values[k][ok][1]
  };d[[v]]<-out
 }
 dd<-make_data('Total')
 if(sens=='omit_filled_covariates')dd$terms<-setdiff(dd$terms,dict$safe[dict$code%in%names(filled)])
 if(sens=='observed2020_2023'){dd$data<-dd$data[dd$data$year<=2023,];dd$terms<-setdiff(dd$terms,dict$safe[dict$code%in%names(filled)[sapply(filled,function(y)any(y<=2023))]])}
 if(grepl('drop_year_',sens))dd$data<-dd$data[dd$data$year!=as.integer(sub('drop_year_','',sens)),]
 if(sens=='drop_Anseong2024')dd$data<-dd$data[!(dd$data$region=='경기도안성시'&dd$data$year==2024),]
 if(sens=='drop_both_districts')dd$data<-dd$data[!dd$data$region%in%c('전라북도남원시','경기도안성시'),]
 runfit(dd,paste0('Total_',sens))
}
d<-base_d
base_bym<-pc_bym;base_prec<-pc_prec
for(setting in c('pc_sd03','pc_sd10','pc_tail05')){
 u<-switch(setting,pc_sd03=.3,pc_sd10=1,pc_tail05=.5);alpha<-if(setting=='pc_tail05').05 else .01
 pc_bym<-list(prec.unstruct=list(prior='pc.prec',param=c(u,alpha)),prec.spatial=list(prior='pc.prec',param=c(u,alpha)));pc_prec<-list(prec=list(prior='pc.prec',param=c(u,alpha)));runfit(make_data(),paste0('Total_',setting))
}
pc_bym<-base_bym;pc_prec<-base_prec
randoms[['BYM2']]<-"+f(idarea,model='bym2',graph=g_main,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(.5,.01)),phi=list(prior='pc',param=c(.5,.5))))+f(idarea_time,model='iid',hyper=pc_prec)"
fitb<-runfit(make_data(),'Total_BYM2','BYM2');write.csv(fitb$summary.hyperpar,file.path(O,'bym2_hyperparameters.csv'))
# Alternative adjacency graphs use projected district centroids; same 223 districts.
suppressPackageStartupMessages(library(sf));sh<-x$shp_main;if(is.na(st_crs(sh)))sh<-st_set_crs(sh,5179);coords<-st_coordinates(st_centroid(st_transform(sh,5179)));graphs<-list(Rook=poly2nb(sh,queen=FALSE,snap=.01));for(k in 2:7)graphs[[paste0('kNN',k)]]<-knn2nb(knearneigh(coords,k=k),sym=TRUE)
g_orig<-g_main
for(nm in names(graphs)){gp<-file.path(O,paste0('graph_',nm,'.graph'));nb2INLA(file=gp,nb=graphs[[nm]]);g_main<-inla.read.graph(gp);runfit(make_data(),paste0('Total_graph_',nm))};g_main<-g_orig
# Updated residual spatial summaries and directional probabilities.
fit<-readRDS(file.path(O,'fits/Total_corrected__M4__cases__1.rds'));dd<-make_data()$data;res<-dd$cases-fit$summary.fitted.values$mean[seq_len(nrow(dd))];agg<-aggregate(res,list(region=dd$region),mean);rv<-agg$x[match(sh$region,agg$region)];rates<-aggregate(cbind(cases,population)~region,dd,sum);cr<-rates$cases/rates$population*1e5;cr<-cr[match(sh$region,rates$region)];lw<-nb2listw(x$nb_obj,style='W');mp<-moran.test(cr,lw);mr<-moran.test(rv,lw)
write.csv(data.frame(stage=c('crude','residual'),I=c(mp$estimate[1],mr$estimate[1]),p=c(mp$p.value,mr$p.value)),file.path(O,'moran_updated.csv'),row.names=FALSE)
write.csv(data.frame(region=sh$region,crude=cr,residual=rv),file.path(O,'moran_plot_data.csv'),row.names=FALSE)
sp<-fit$summary.random$idarea[seq_len(nrow(sh)),];write.csv(data.frame(region=sh$region,sp),file.path(O,'posterior_district_effects.csv'),row.names=FALSE)
probs<-list();for(group in c('Total','Urban','Rural')){fit<-readRDS(file.path(O,'fits',paste0(group,'_corrected__M4__cases__1.rds')));for(tm in names(fit$marginals.fixed)){pr<-1-inla.pmarginal(0,fit$marginals.fixed[[tm]]);probs[[paste(group,tm)]]<-data.frame(group=group,term=tm,Pr_positive=pr,Pr_direction=max(pr,1-pr))}}
write.csv(do.call(rbind,probs),file.path(O,'posterior_direction_probabilities.csv'),row.names=FALSE);write.csv(do.call(rbind,scale_records),file.path(O,'moisture_scaling_additional.csv'),row.names=FALSE)
cat('ADDITIONAL CHECKS COMPLETE\n')
