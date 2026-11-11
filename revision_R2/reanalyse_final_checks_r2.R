p<-file.path(Sys.getenv("R2_CODE_ROOT","."),"reanalyse_additional_r2.R");s<-paste(readLines(p,warn=FALSE),collapse='\n');eval(parse(text=strsplit(s,"for(group in c('Total','Urban','Rural'))runfit",fixed=TRUE)[[1]][1]))
# Repeat the paediatric battery with the identical covariate specification.
for(group in c('Total','Urban','Rural')){
 dd<-make_data(group);z<-dd$data;ped<-c(child_0_4_z='영유아비율_0_4',child_5_9_z='아동비율_5_9',ped_gastro_adm_z='예방가능입원율_소아위장관염',daycare_per1k_z='유아천명당보육시설수')
 for(tm in names(ped)){v<-as.numeric(z[[ped[[tm]]]]);form<-x$FIX_FORMS[[ped[[tm]]]];if(!is.null(form)&&form=='log1p')v<-log1p(pmax(v,0));z[[tm]]<-as.numeric(scale(v))}
 for(tm in c('livestock_total_z','fiscal_auto_z'))if(!tm%in%names(z)){cd<-dict$code[match(tm,dict$safe)];v<-as.numeric(z[[cd]]);if(grepl('사육두수',cd))v<-v/10000;z[[tm]]<-as.numeric(scale(v))}
 ts<-c('livestock_total_z','sludge_moisture_z','moisture_observed','fiscal_auto_z',names(ped));z<-z[complete.cases(z[,ts]),]
 repeat{vv<-car::vif(lm(as.formula(paste('cases~',paste(ts,collapse='+'))),data=z));if(max(vv)<5||length(ts)<3)break;ts<-setdiff(ts,names(which.max(vv)))}
 write.csv(z[,c('region','year','cases','population',ts)],file.path(O,paste0('design_paediatric_',group,'.csv')),row.names=FALSE)
 runfit(list(data=z,terms=ts),paste0(group,'_paediatric'),repeat_id=2)
}
for(sens in c('full','drop_Namwon2024','drop_both_districts')){
 dd<-make_data();dd$data$urban<-as.integer(dd$data$area_2=='도시');dd$terms<-c(dd$terms,'urban')
 for(v in c('livestock_total_z','sludge_moisture_z','ww_effluent_z')){nm<-paste0(v,'_urban');dd$data[[nm]]<-dd$data[[v]]*dd$data$urban;dd$terms<-c(dd$terms,nm)}
 if(sens=='drop_Namwon2024')dd$data<-dd$data[!(dd$data$region=='전라북도남원시'&dd$data$year==2024),]
 if(sens=='drop_both_districts')dd$data<-dd$data[!dd$data$region%in%c('전라북도남원시','경기도안성시'),]
 runfit(dd,paste0('Total_interaction_',sens))
}
# Regularisation check of the corrected principal fixed-effect set (no reselection).
base_safe<-safe_inla
safe_inla<-function(formula,...)base_safe(formula,...,control.fixed=list(prec=1))
for(group in c('Total','Urban','Rural'))runfit(make_data(group),paste0(group,'_ridge'))
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo_final.txt'));cat('FINAL CHECKS COMPLETE\n')
