# Run after gcpo.R and the direct-refit scripts.
suppressPackageStartupMessages(library(INLA))
B<-normalizePath(Sys.getenv('R2_CODE_ROOT','.')); O<-file.path(B,'output/cpo_diagnostics')
rows<-list()
for(group in c('Total','Urban','Rural')) {
 f<-readRDS(file.path(B,'output/fits',paste0(group,'_corrected__M4__cases__1.rds')))
 z<-readRDS(file.path(O,paste0('gcpo_',group,'.rds')))
 idx<-which(f$cpo$failure>0)
 if(group=='Total') idx<-sort(unique(c(idx,455L,578L)))
 stopifnot(length(idx)==switch(group,Total=17L,Urban=11L,Rural=1L))
 rows[[group]]<-data.frame(group=group,index=idx,gcpo=z$gcpo$gcpo[idx],old_failure=f$cpo$failure[idx])
}
write.csv(do.call(rbind,rows),file.path(O,'gcpo_all_warning_locations.csv'),row.names=FALSE)
