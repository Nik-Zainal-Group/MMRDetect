
#' Substitution signature of MMR gene knockouts
#'
#' @format A data frame with nine variables:
#' \describe{
#' \item{\code{MutationType}}{96 substitution types}
#' \item{\code{MLH1}}{MLH1 KO}
#' \item{\code{MSH2}}{MSH2 KO}
#' \item{\code{MSH6}}{MSH6 KO}
#' \item{\code{PMS2}}{PM2 KO}
#' }
#'
#'
"MMRKO_subsig"

#' Indel signature of MMR gene knockouts
#'
#' @format A data frame with nine variables:
#' \describe{
#' \item{\code{indelsubtype}}{Indel types}
#' \item{\code{MLH1}}{MLH1 KO}
#' \item{\code{MSH2}}{MSH2 KO}
#' \item{\code{MSH6}}{MSH6 KO}
#' \item{\code{PMS2}}{PM2 KO}
#' }
#'
#'
"MMRKO_indelsig"



#' Compute substitution signature exposure for a given catalogue
#' 
#' @param sub_cat substitution catalogue
#' @param tissue_type tissue type
#' @param MMR_subsig96 mismatch repair gene knockout substitution signatures
#' @param tissue_subsig96 tissue-specific substitution signatures
#' 
#' @return Substitution signature exposure
#' 
#' @export
MMRDetect.compute.subcatalogue.exposure <- function(sub_cat, tissue_type,MMR_subsig96=MMRKO_subsig,tissue_subsig96){
  
  
  # Generate catalouge for subs_highburden
  sub_catalouge <- sub_cat
  
  
  # Step 2:  Signature fitting for subs using Andrea's tissue-specific signatures 
  
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  MMR1Sig <- c("Breast_A", "Colorectal_F", "Liver_E","Stomach_H", "Uterus_C", "Uterus_J")
  MMR2Sig <- c("Biliary_E", "Breast_D", "Colorectal_E", "Ovary_F", "Pancreas_H", "Stomach_B", "Uterus_D", "Kidney_D")
  
  tissue_MMR1Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR1Sig)]
  tissue_MMR2Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR2Sig)]
  
  
  if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)>0){
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)==0){
    # add PMS2
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)>0){
    # add MLH1
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR1")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)==0){
    
    # add MLH1 and PMS2 
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    names(selected_tissueSig)[dim(selected_tissueSig)[2]-1] <- c("MMR1")
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  
  # Choose the samples that have exposure in MMR signatures
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% c(tissue_MMR1Sig, tissue_MMR2Sig,"MMR2","MMR1")),]
  
  MMRsig_sample$AndreaSig <- rownames(MMRsig_sample)
  MMRsig_sample_melt <- reshape2::melt(MMRsig_sample,c("AndreaSig"))
  names(MMRsig_sample_melt) <- c("AndreaSig","Sample","exposure")
  MMRsig_sample_melt_dcast <- reshape2::dcast(MMRsig_sample_melt,Sample~AndreaSig,value.var="exposure")
  
  if(dim(MMRsig_sample_melt_dcast)[2]==2){
    MMRsig_sample_melt_dcast$MMR_sum <- MMRsig_sample_melt_dcast[,-1]
  }else{
    MMRsig_sample_melt_dcast$MMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
    
  }
  write.table(MMRsig_sample_melt_dcast,paste0("subsig_exposure_",tissue_type,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  return(MMRsig_sample_melt_dcast) 
}

#' Compute cosine similarity between given indel profile and MMR gene knockout indel profilesct
#' 
#' @param indels list of indels
#' @param tissue_type tissue type
#' @param MMR_sig_indel mismatch repair gene knockout indel signatures
#' 
#' @return cosine similarity between given indel list and MMR gene knockout indel profiles
#' 
#' @export
MMRDetect.compute.Repindel.similarity <- function(indels, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  indel_classied <- indel_classifier(indels)
  
  indel_classied_rep <- indel2610_classied[!indel2610_classied$indeltype_short%in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),]
  Sample_MMR <- data.frame(table(indel_classied_rep$Sample))
  names(Sample_MMR) <- c("Sample","RepIndel_num")
  
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  # Sample_MMR <- data.frame(table(indels$Sample))
  #  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- plyr::ddply(cossim_allsample,c("Sample"),plyr::summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- plyr::ddply(cossim_allsample,c("Sample"),plyr::summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  MMRsig_2 <- merge(Sample_MMR,MMRsig_2, by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}

#' Compute cosine similarity between given indel profile and MMR gene knockout indel profiles
#' 
#' @param indel_cat  indels catalogue
#' @param tissue_type tissue type
#' @param MMR_sig_indel mismatch repair gene knockout indel signatures
#' 
#' @return cosine similarity between given indel profile and MMR gene knockout indel profiles
#' 
#' @export
MMRDetect.compute.Repindelcatalogue.similarity <- function(indel_cat, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  # Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  
  Sample_MMR <- data.frame("Sample"=names(indel_cat)[3:dim(indel_cat)[2]], "RepIndel_num"=colSums(indel_cat[!indel_cat$indelsubtype %in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),3:dim(indel_cat)[2]]))
  names(Sample_MMR) <- c("Sample","RepIndel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_cat[,c("indelsubtype", "type",as.character(Sample_MMR[j,"Sample"]))]
      
      cossim <- Calculae_Cossim_catalogue_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
      cossim_allsample <- rbind(cossim_allsample,cossim)
      cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
      cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
      
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- plyr::ddply(cossim_allsample,c("Sample"),plyr::summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- plyr::ddply(cossim_allsample,c("Sample"),plyr::summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  MMRsig_2 <- merge(Sample_MMR,MMRsig_2, by="Sample")
  
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  
  return(MMRsig_2) 
}

#' Compute cosine similarity between given sub profiles and MMR gene knockout sub profiles
#' 
#' @param mcat  substitution catalogue
#' @param scat mismatch repair gene knockout sub signatures
#' 
#' @return cosine similarity between given sub profiles and MMR gene knockout sub profiles
#' 
#' @export
MMRDetect.compute.subcatalogue.similarity <- function(mcat,scat=MMRKO_subsig){
  mut_sig <- merge(scat,mcat,by="MutationType")
  sig_cat <- mut_sig[,2:dim(scat)[2]]
  mut_cat <- mut_sig[,(dim(scat)[2]+1):dim(mut_sig)[2]]
  mut_cat <- mut_cat/colSums(mut_cat)[col(mut_cat)]
  a <- apply(sig_cat, 2, function(x) SigCossim(mut_cat,x))
  subsim <- as.data.frame(a)
  subsim$Sample <- rownames(subsim)
  subsim$Sample <- sub(".DNA","-DNA",subsim$Sample)
  subsim$maxcossim <- apply(subsim[,-which(names(subsim)=="Sample")],1,max)
  
  return(subsim)
  
}

#' Compute cosine similarity between one sample and MMR gene knockouts
#' 
#' @param SingleSample_classified_indels  classified indel list
#' @param Sig mismatch repair gene knockout indel signatures
#' 
#' @return cosine similarity between given sub profiles and MMR gene knockout sub profiles
#' 
#' @export
Generate_CossimVector_SingleSample_RepIndel <- function(SingleSample_classified_indels, Sig=MMRKO.indelsig){
  mut_catalogue <-  gen_indelmuttype_MMRD(SingleSample_classified_indels,"Sample","indeltype_short")
  return(Calculae_Cossim_catalogue_RepIndel(mut_catalogue, Sig))
}

#' Compute cosine similarity between one sample indel profile and MMR gene knockouts
#' 
#' @param singlesample_indel_catalogue   indel catalogue of one sample
#' @param Sig mismatch repair gene knockout indel signatures
#' 
#' @return cosine similarity between given sub profiles and MMR gene knockout sub profiles
#' 
#' @export
Calculae_Cossim_catalogue_RepIndel <- function(singlesample_indel_catalogue, Sig){
  total_subs_sig <- merge(singlesample_indel_catalogue,Sig, by="indelsubtype")
  
  cossim_del <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Del",3],total_subs_sig[total_subs_sig$type=="Del",4]))
  cossim_ins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Ins",3],total_subs_sig[total_subs_sig$type=="Ins",4]))
  cossim <- abs(cos_similarity(total_subs_sig[,3],total_subs_sig[,4]))
  
  return(c(cossim_del, cossim_ins,cossim))
}

#' Compute cosine similarity between two vectors
#' 
#' @param v1   vector 1
#' @param v2 vector 2
#' 
#' @return cosine similarity
#' 
#' @export
cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}

#' Compute cosine similarity between a set of catalogue and a given signatures
#' 
#' @param mcat   mutation catalogue
#' @param sigx a given signature
#' 
#' @return a list of cosine similarity
#' 
#' @export
SigCossim <- function(mcat, sigx){
  
  return(apply(mcat, 2, function(x) cos_similarity(x,sigx)))
}

#' Compute variables for MMRDetect
#' 
#' @param sub_cat substitution catalogue
#' @param indel_cat indel catalogue
#' @param tissue_type tissue type
#' @param MMR_subsig96 mismatch repair gene knockout substitution signatures
#' @param MMR_sig_indel mismatch repair gene knockout indel signatures
#' @param tissue_subsig96 tissue-specific substitution signatures
#' 
#' @return variables
#' 
#' @export
MMRDetect.compute.variables <- function(sub_cat, indel_cat, tissue_type,MMR_subsig96=MMRKO_subsig,MMR_sig_indel=MMRKO_indelsig, tissue_subsig96=PancanSig){
  
  sub_MMR_expo <- MMRDetect.compute.subcatalogue.exposure(sub_cat, tissue_type,MMR_subsig96,tissue_subsig96)
  indel_similarity <- MMRDetect.compute.Repindelcatalogue.similarity(indel_cat, tissue_type,MMR_sig_indel)
  sub_similiarity <- MMRDetect.compute.subcatalogue.similarity(sub_cat, MMR_subsig96)
  
  MMRDetect_variables <- merge(sub_MMR_expo[,c("Sample","MMR_sum")],indel_similarity,by="Sample")
  MMRDetect_variables <- merge(MMRDetect_variables,sub_similiarity[,c("Sample","maxcossim")],by="Sample")
  write.table(MMRDetect_variables,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  return(MMRDetect_variables) 
}
