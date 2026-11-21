suppressPackageStartupMessages(library(INLA));B<-normalizePath(Sys.getenv("R2_CODE_ROOT","."));O<-file.path(B,'output/cpo_diagnostics');x<-readRDS(file.path(B,'input/analysis_input.rds'));g_main<-x$g_main
pc_bym<-list(prec.unstruct=list(prior='pc.prec',param=c(.5,.01)),prec.spatial=list(prior='pc.prec',param=c(.5,.01)));pc_prec<-list(prec=list(prior='pc.prec',param=c(.5,.01)))

dir.create(O,recursive=TRUE,showWarnings=FALSE)
for(group in c('Total','Urban','Rural')){
 f<-readRDS(file.path(B,'output/fits',paste0(group,'_corrected__M4__cases__1.rds')));a<-f$.args;a$control.compute<-list(config=TRUE,dic=TRUE,waic=TRUE,control.gcpo=list(enable=TRUE,num.level.sets=-1,correct.hyperpar=TRUE));a$control.mode<-list(theta=f$mode$theta,restart=TRUE);a$num.threads<-'2:1';a$verbose<-FALSE;cat('START',group,'\n');flush.console();z<-do.call(INLA::inla,a);saveRDS(z,file.path(O,paste0('gcpo_',group,'.rds')));print(names(z$gcpo));print(head(z$gcpo));cat('DONE',group,'\n');flush.console()
}
