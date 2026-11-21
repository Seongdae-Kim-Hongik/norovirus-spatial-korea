suppressPackageStartupMessages(library(INLA));B<-normalizePath(Sys.getenv("R2_CODE_ROOT","."));O<-file.path(B,'output/cpo_diagnostics');x<-readRDS(file.path(B,'input/analysis_input.rds'));g_main<-x$g_main
pc_bym<-list(prec.unstruct=list(prior='pc.prec',param=c(.5,.01)),prec.spatial=list(prior='pc.prec',param=c(.5,.01)));pc_prec<-list(prec=list(prior='pc.prec',param=c(.5,.01)))
dir.create(file.path(O,'loo_fits'),showWarnings=FALSE);rows<-list();effects<-list();if(file.exists(file.path(O,'loo_extra_predictive.csv'))){old<-read.csv(file.path(O,'loo_extra_predictive.csv'));rows<-split(old,paste(old$group,old$index,sep='_'))};if(file.exists(file.path(O,'loo_extra_effects.csv'))){old<-read.csv(file.path(O,'loo_extra_effects.csv'));effects<-split(old,paste(old$group,old$heldout,sep='_'))}
for(group in c('Total')){
 f<-readRDS(file.path(B,'output/fits',paste0(group,'_corrected__M4__cases__1.rds')));idxs<-c(455L,578L)
 for(idx in idxs){
  key<-paste(group,idx,sep='_');if(key%in%names(rows)||(idx==0&&key%in%names(effects)))next;cat(format(Sys.time()),'START',key,'\n');flush.console();a<-f$.args;exclude<-if(idx==0)idxs else idx;y<-a$data$cases[exclude];a$data$cases[exclude]<-NA;a$control.compute<-list(config=TRUE,return.marginals.predictor=TRUE,dic=TRUE,waic=TRUE);a$control.mode<-list(theta=f$mode$theta,restart=TRUE);a$num.threads<-'4:1';a$verbose<-FALSE
  dest<-file.path(O,'loo_fits',paste0(key,'.rds'));z<-if(file.exists(dest))readRDS(dest) else do.call(INLA::inla,a);saveRDS(z,dest)
  stopifnot(is.finite(z$dic$dic),z$dic$dic>0,all(is.finite(z$summary.fixed$mean)))
  ef<-z$summary.fixed;orig<-f$summary.fixed;effects[[key]]<-data.frame(group=group,heldout=idx,term=rownames(ef),original_mean=orig$mean,new_mean=ef$mean,original_lo=orig[['0.025quant']],original_hi=orig[['0.975quant']],new_lo=ef[['0.025quant']],new_hi=ef[['0.975quant']]);write.csv(do.call(rbind,effects),file.path(O,'loo_extra_effects.csv'),row.names=FALSE)
  if(idx>0){s<-inla.posterior.sample(20000,z,selection=list(Predictor=idx),num.threads='1:1',seed=261121L);eta<-sapply(s,function(v)v$latent[1]);size<-sapply(s,function(v)v$hyperpar[1]);stopifnot(grepl('size for the nbinomial',names(s[[1]]$hyperpar)[1]),all(size>0),abs(mean(eta)-z$summary.linear.predictor$mean[idx])<(5*sd(eta)/sqrt(length(eta))+.05));v<-exp(dnbinom(y,size=size,mu=exp(eta),log=TRUE));rows[[key]]<-data.frame(group=group,index=idx,region=f$.args$data$region[idx],year=f$.args$data$year[idx],cases=y,old=f$cpo$cpo[idx],new=mean(v),mcse=sd(v)/sqrt(length(v)),eta_mean=mean(eta),eta_summary=z$summary.linear.predictor$mean[idx],offset=z$offset.linear.predictor[idx]);write.csv(do.call(rbind,rows),file.path(O,'loo_extra_predictive.csv'),row.names=FALSE)}
  cat('DONE',key,'DIC',z$dic$dic,'\n');flush.console()
 }
}
cat('COMPLETE\n')
