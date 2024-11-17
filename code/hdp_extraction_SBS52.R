library(hdp)
library(ggplot2)
library(RColorBrewer)
library(lsa)
library(lattice)
chlist_file<-commandArgs(T)[1]
hdp_input_path<-commandArgs(T)[2]
output_path<-commandArgs(T)[3]
prefix<-commandArgs(T)[4]
#------------extract HDP results---------
chlist <- vector("list", 10)
for (i in 1:10){
	  chlist[[i]] <- readRDS(paste0(chlist_file,i,".rds"))
}

mut_example_multi <- hdp_multi_chain(chlist)


mut_example_multi <- hdp_extract_components(mut_example_multi)

hdp_exposures=mut_example_multi@comp_dp_distn[["mean"]][2:dim(mut_example_multi@comp_dp_distn[["mean"]])[1],]
input_for_hdp = read.table(hdp_input_path,check.names = F,header=T)
input_for_hdp_sum <- apply(input_for_hdp,1,sum)
input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>0,]
rownames(hdp_exposures)[(nrow(hdp_exposures)-nrow(input_for_hdp)+1):nrow(hdp_exposures)]=rownames(input_for_hdp)
write.csv(hdp_exposures,paste0(output_path, "/",prefix,"_HDP_exposure.csv"),quote = F)


hdp_sigs=data.frame(t(mut_example_multi@comp_categ_distn[["mean"]][1:dim(mut_example_multi@comp_categ_distn[["mean"]])[1],]))
categories = paste0(substr(colnames(input_for_hdp),5,5),'[',substr(colnames(input_for_hdp),1,3),']',substr(colnames(input_for_hdp),7,7))
rownames(hdp_sigs)<-colnames(categories)
write.csv(hdp_sigs,paste0(output_path, "/",prefix,"_HDP_sigs.csv"),quote = F)

