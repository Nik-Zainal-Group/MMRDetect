#library(ggplot2)
#library(plyr) 
#library(reshape2)
#library(gridExtra)
#library(scales)
#library(factoextra)
#library("VariantAnnotation")
#library('GenomicFeatures')
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(signature.tools.lib)
#library("glmnet")

#load("../data/MMRKO_indelsigCT.rda")
#load("../data/indelsig_template.rda")
#load("../data/MMRKO_indelsig.rda")
#load("../data/MMRKO_subsig.rda")
#load("../data/PancanSig.rda")
#MMRDclassifier <- readRDS("../data/MMRDetect.rds")

#' Generate 96 channel catalogue for substitutions
#'
#' @param CTsubs A list of substitutions: Chrom, Pos, Ref, Alt, Sample column
#' @param SampleCol Sample column
#' @return 96 channel catalogue for substitutions
#' @export
GenCatalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-1, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+1))
  
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  base <- c("A","C","G","T")
  mutationtype <- NULL
  for(i in 1:length(mutation)){
    for(j in 1:length(base)){
      for(k in 1:length(base)){
        
        mutationtype <- c(mutationtype,paste0(base[j],"[",mutation[i],"]",base[k]))
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("MutationType"=mutationtype)
 # muttype_freq_template$Mutation <- substr(muttype_freq_template$MutationType,3,5)
  #read.table("./MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(Biostrings::complement(Biostrings::DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(Biostrings::complement(Biostrings::DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(Biostrings::complement(Biostrings::DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(Biostrings::complement(Biostrings::DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-reshape2::dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}


##########################
#
#  Indel classification
#
##########################
#' Generate indel catalogue for substitutions
#'
#' @param CTsubs A list of substitutions: Chrom, Pos, Ref, Alt, Sample column
#' @param SampleCol Sample column
#' @return 96 channel catalogue for substitutions

#' @export
indel_classifier <- function(indels){
  
  indel.data <- indels
  indel.data[indel.data$Chrom=="23","Chrom"]="X"
  indel.data[indel.data$Chrom=="24","Chrom"]="Y"
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df.tab(indel.data)
  indel.df.max100 <- indel.df[indel.df$indel.length<=100 & indel.df$indel.length>0,]
  # indel classification
  indel.classified.df <- mh_indel_v2(indel.df.max100)
  
  # remove the insertions with repcount =10
  indel.classified.df <- indel.classified.df[indel.classified.df$repcount<10,]
  
  indel.classified_details <- subtype_indel_mmrd(indel.classified.df)
  
  indel.classified_details[indel.classified_details$Chrom=="X","Chrom"]=23
  indel.classified_details[indel.classified_details$Chrom=="Y","Chrom"]=24
  indel_template2 <- indelsig_template
  names(indel_template2)[1] <- c("Subtype")
  indel.classified_details <- merge(indel.classified_details,indel_template2,by="Subtype")
  
  return(indel.classified_details)
  
}

prepare.indel.df.tab <- function(indel.data) {
  
  if (nrow(indel.data)>0) {
    
    
    indel.data$ref.length <- nchar(indel.data$Ref)
    indel.data$alt.length <- nchar(indel.data$Alt)
    indel.data$indel.length <- abs(indel.data$ref.length - indel.data$alt.length)
    
    #
    indel.data$Type <- NULL
    indel.data[indel.data$ref.length==1,"Type"] <- "Ins"
    indel.data[indel.data$alt.length==1,"Type"] <- "Del"
    indel.data[indel.data$alt.length!=1 & indel.data$ref.length!=1,"Type"] <- "Complex"
    
    # sequence of change
    indel.data$change <- NULL
    indel.data[indel.data$Type=="Complex","change"] <- substr( as.character(indel.data[indel.data$Type=='Complex',"Ref"]),1,1e5)
    indel.data[indel.data$Type=="Ins","change"] <- substr( as.character(indel.data[indel.data$Type=='Ins',"Alt"]),2,1e5)
    indel.data[indel.data$Type=="Del","change"] <- substr( as.character(indel.data[indel.data$Type=='Del',"Ref"]),2,1e5)
    
    # 5'-neighbor base
    indel.data$pre_context <- substr( as.character(indel.data[,"Ref"]),1,1)
    
    # 27bp before and after change                     
    indel.data$extend5 = indel.data$Pos-indel.data$indel.length-25;
    indel.data$extend3 = indel.data$Pos+indel.data$indel.length + indel.data$indel.length+25;
    
    
    indel.data$slice5 <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data$Chrom), indel.data$extend5, indel.data$Pos))
    indel.data$slice3 <- NULL
    indel.data[indel.data$Type=="Del","slice3"] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","extend3"]))
    indel.data[indel.data$Type=="Ins","slice3"] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","extend3"]))
    
    # 1bp before and after change                     
    indel.data$slice5_1bp <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos, indel.data$Pos))
    indel.data$slice3_1bp <- NULL
    indel.data[indel.data$Type=="Del","slice3_1bp"] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1))
    indel.data[indel.data$Type=="Ins","slice3_1bp"] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+1))
    
   
    # Pyrimidine represnetation for 1bp indels
    indel.data$change.pyr <- indel.data$change
    indel.data$slice3_1bp_pyr <- indel.data$slice3_1bp
    indel.data$slice5_1bp_pyr <- indel.data$slice5_1bp
    
    #    indel.data$slice3_2bp_pyr <- indel.data$slice3_2bp
    #    indel.data$slice5_2bp_pyr <- indel.data$slice5_2bp
    
    indel.data[indel.data$change=="A","change.pyr"] <- "T"
    indel.data[indel.data$change=="A","slice5_1bp_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.data[indel.data$change=="A","slice3_1bp"])))
    indel.data[indel.data$change=="A","slice3_1bp_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.data[indel.data$change=="A","slice5_1bp"])))
    
  
    
    indel.data[indel.data$change=="G","change.pyr"] <- "C"
    indel.data[indel.data$change=="G","slice5_1bp_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.data[indel.data$change=="G","slice3_1bp"])))
    indel.data[indel.data$change=="G","slice3_1bp_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.data[indel.data$change=="G","slice5_1bp"])))
    
   
    return(indel.data)
    
    
  }
}

#' find mh indels
#'
#' @param indel.df A list of indels with flanking sequence information: Chrom, Pos, Ref, Alt, Sample column
#' @return annotated indel list

#' @export
mh_indel_v2 <- function(indel.df) {
  
  # indel.df needs following columns:
  # indel.type
  # change
  # slice3
  # slice5
  # indel.length
  
  
  classification <- rep (NA, nrow(indel.df))
  repcount_all <- rep (0, nrow(indel.df))
  mhcount_all <- rep (0, nrow(indel.df))
  slice3_nonrep_all <- rep (NA, nrow(indel.df))
  if (nrow(indel.df)>0) {
    for (i in 1:nrow(indel.df)) {
      print(i)
      if (indel.df$Type[i]=='Del') { # the classification is only for deletions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        #Look for Microhomology first and then for repeats - tandem/normal
        r_ref = repcaller_del(as,bs,cs,indel.df$indel.length[i])
        repcount_ref <- r_ref$countrep # number of times there is a repeat
        
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        classification[i] <- finalcall
        repcount_all[i] <- repcount_ref
        mhcount_all[i] <- mhcount
        slice3_nonrep_all[i] <- substr(bs,repcount*nchar(rept)+1,repcount*nchar(rept)+1)
        #print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall
      } # end if
      
      if (indel.df$Type[i]=='Ins') { # the classification is only for insertions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        
        classification[i] <- finalcall
        repcount_all[i] <- repcount
        mhcount_all[i] <- mhcount
        slice3_nonrep_all[i] <- substr(bs,repcount*nchar(rept)+1,repcount*nchar(rept)+1)
        
      }
    } # end for
    
    classification[is.na(classification)] <- '-'
    slice3_nonrep_all[is.na(slice3_nonrep_all)] <- '-'
    #  mhcount_all[is.na(mhcount_all)] <- '-'
    indel.df$repcount <- repcount_all
    indel.df$mhcount <- mhcount_all
    indel.df$slice3_nonrep <- slice3_nonrep_all
    indel.df$slice3_nonrep_pyr <- indel.df$slice3_nonrep
    indel.df$slice5_nonrep_pyr <- indel.df$slice5_1bp_pyr
    indel.df[indel.df$change=="A","slice5_nonrep_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.df[indel.df$change=="A","slice3_nonrep"])))
    indel.df[indel.df$change=="A","slice3_nonrep_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.df[indel.df$change=="A","slice5_1bp"])))
    indel.df[indel.df$change=="G","slice5_nonrep_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.df[indel.df$change=="G","slice3_nonrep"])))
    indel.df[indel.df$change=="G","slice3_nonrep_pyr"] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(indel.df[indel.df$change=="G","slice5_1bp"])))
    
    
    indel.df$classification <- classification
  } else {
  }
  
  indel.df
}  # end the main mh function  


#' classify indels 
#'
#' @param indel.df A list of indels with flanking sequence information: Chrom, Pos, Ref, Alt, Sample column
#' @return 96 channel catalogue for substitutions

#' @export
subtype_indel_mmrd <- function(indel.df) {
  
  indel.df$Subtype <- NULL
  indel.df$TypeS <- ""
  indel.df[indel.df$Type=="Ins","TypeS"] <- "+"
  indel.df[indel.df$Type=="Del","TypeS"] <- "-"
  indel.df[indel.df$Type=="Complex","TypeS"] <- "Complex"
  
  indel.new <- NULL
  #[+C], [+T], [-C], [-T]
  
  #  1 rep
  indel.c <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1)
  if(dim(indel.c)[1]>0){
    indel.c$Subtype <- paste0(indel.c$slice5_nonrep_pyr,"|[",indel.c$TypeS,indel.c$change.pyr,"]Rep=",indel.c$repcount,"|",indel.c$slice3_nonrep_pyr)
    indel.c$classification <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]")
    indel.new <- rbind(indel.new,indel.c)
    
  }

  # [+>1]Mh, [->1]Mh
  indel.e1 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated")
  if(dim(indel.e1)[1]>0){
    indel.e1$Subtype <-paste0("[",indel.e1$TypeS,"]Mh")
    indel.e1$classification <- paste0("[",indel.e1$TypeS,"]Mh")
    indel.new <- rbind(indel.new,indel.e1)
    
  }
   
  
  # [+>1]Rep, [->1]Rep
  indel.h <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount<5)
  if(dim(indel.h)[1]>0){
    indel.h$Subtype <-paste0("[",indel.h$TypeS,">1]Rep<=4")
    indel.h$classification <- paste0("[",indel.h$TypeS,">1]Rep")
    indel.new <- rbind(indel.new,indel.h)
    
  }
 
  indel.i <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount>=5)
  if(dim(indel.i)[1]>0){
    indel.i$Subtype <-paste0("[",indel.i$TypeS,">1]Rep>=5")
    indel.i$classification <- paste0("[",indel.i$TypeS,">1]Rep")
    indel.new <- rbind(indel.new,indel.i)
    
  }
 
  # [+>1]Other, [->1]Other
  indel.j <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="None")
  if(dim(indel.j)[1]>0){
    indel.j$Subtype <-paste0("[",indel.j$TypeS,">1]Other")
    indel.j$classification <- paste0("[",indel.j$TypeS,">1]Other")
    indel.new <- rbind(indel.new,indel.j)
  }
  
  
  # Complex
  indel.k <- subset(indel.df,Type=="Complex")
  if(dim(indel.k)[1]>0){
    indel.k$Subtype <-"Complex"
    indel.k$classification <- "Complex"
    indel.new <- rbind(indel.new,indel.k)
    
  }
  
  
  #a=data.frame(table(indel.new$VariantID))
  #a=a[order(a$Freq,decreasing=T),]
  return(indel.new)
}

#' micro-homology between the deleted sequence and 3' sequence (number of bases)
#' 
#' @param d A list of indels with flanking sequence information: Chrom, Pos, Ref, Alt, Sample column
#' @return annotated indels

#' @export
mhcaller <- function(d, prime3) {
  # d is the change of a deletion
  # prime 3 is the sequence 3' of the deleted segment
  countmh = 0
  seq = '-'
  # Dealing with microhomology or lack of microhmomology in the first position
  
  
  for (i in 1:(nchar(d))) { # for all substrings
    
    if (substr(d, 1, i)==substr(prime3,1,i)) {
      countmh <- countmh + 1
      seq <- substr(d, 1, i)
    }
  }
  result <- list()
  result$countmh <- countmh # number of basepairs matching between the deleted fragment and the 3' segment
  result$seq <- seq
  result
}

repcaller <- function(d, prime3, prime5, l) {
  # d : deletion
  # prime3 : 3' context
  # prime5: 5' context
  # l : length of change
  
  result <- list() 
  result$countrep = 0
  result$rep<-''
  # This is for counting single base repeats
  # if the length of change is 1
  if (l==1) {
    countrep <- 0
    i <- 1
    while (substr(d,1,1)==substr(prime3,i,i)) {
      countrep <- countrep +1
      i <- i + 1
    }
    result$countrep <- countrep
    result$rep <- d
    return(result)
  } else if (d==substr(prime3, 1, nchar(d))) { # This is for counting whole deletion/DI repeats that are present in the flanking region  
    
    countrep = tandemcount_v2(d,prime3)
    rep = findsmallestrep(d)
    countrep = max(countrep,tandemcount_v2(rep,prime3))
    result$countrep <- countrep
    result$rep <- rep
    return(result)
  } else {   # This is for counting anything in between 1bp and whole del repetition                                     
    rep = '-'
    
    for (t in seq(from=(nchar(d)-1), to=2)) {  # Look for repeats of 2bp to n-1 bp of the length of the indel
      if (grepl(substring(d,1,t), prime3)) {
        countrep = tandemcount_v2(substr(d, 1, t),prime3)
        rep = findsmallestrep(substr(d, 1, t))
        unit = tandemcount_v2(rep,d)*nchar(rep)
        # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
        # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated. 
        # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false
        if (substr(d,1,unit) == substr(d, unit+1,unit*2)) {
          #print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
          countrep = max(countrep, tandemcount_v2(rep,prime3))         
          result$countrep <- countrep
          result$rep <- rep
          return(result)
        } else {
          result$countrep <- 0
          result$rep <- '-'
          return(result)
        }
      } # endif
    } # end for
    result$countrep <- 0 # in case no repeat 2bp or longer is fount
    result$rep <- '-'
    return(result)
  } # end else
} # end repcaller
repcaller_del <- function(d, prime3, prime5, l){
  # d : deletion
  # prime3 : 3' context
  # prime5: 5' context
  # l : length of change
  
  prime3 <- paste0(d,prime3)
  result <- list() 
  result$countrep = 0
  result$rep<-''
  # This is for counting single base repeats
  # if the length of change is 1
  if (l==1) {
    countrep <- 0
    i <- 1
    while (substr(d,1,1)==substr(prime3,i,i)) {
      countrep <- countrep +1
      i <- i + 1
    }
    result$countrep <- countrep
    result$rep <- d
    return(result)
  } else if (d==substr(prime3, 1, nchar(d))) { # This is for counting whole deletion/DI repeats that are present in the flanking region  
    
    countrep = tandemcount_v2(d,prime3)
    rep = findsmallestrep(d)
    countrep = max(countrep,tandemcount_v2(rep,prime3))
    result$countrep <- countrep
    result$rep <- rep
    return(result)
  } else {   # This is for counting anything in between 1bp and whole del repetition                                     
    rep = '-'
    
    for (t in seq(from=(nchar(d)-1), to=2)) {  # Look for repeats of 2bp to n-1 bp of the length of the indel
      if (grepl(substring(d,1,t), prime3)) {
        countrep = tandemcount_v2(substr(d, 1, t),prime3)
        rep = findsmallestrep(substr(d, 1, t))
        unit = tandemcount_v2(rep,d)*nchar(rep)
        # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
        # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated. 
        # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false
        if (substr(d,1,unit) == substr(d, unit+1,unit*2)) {
          #print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
          countrep = max(countrep, tandemcount_v2(rep,prime3))         
          result$countrep <- countrep
          result$rep <- rep
          return(result)
        } else {
          result$countrep <- 0
          result$rep <- '-'
          return(result)
        }
      } # endif
    } # end for
    result$countrep <- 0 # in case no repeat 2bp or longer is fount
    result$rep <- '-'
    return(result)
  } # end else
} # end repcaller
# tandemcount_v2, by Xueqing, only counts continuous repeated patterns
tandemcount_v2 <- function(pat,string) {
  sum = 0
  
  start.pos <- seq(from=1, to=nchar(string)-nchar(pat)+1, by=nchar(pat))
  s = 1
  
  while (substr(string, s, s+nchar(pat)-1) == pat) {
    sum <- sum + 1
    s = s + nchar(pat)
  }
  
  sum
  
}

#' This finds the smallest repeating subunit of the deletion
findsmallestrep <- function(d) {
  d.factors <- as.numeric(sort(unlist(get_all_factors(nchar(d))), decreasing=TRUE))
  rep.unit <- ''
  
  for (f in d.factors) {
    no.repeats <- nchar(d)/f
    rep.string <- paste0(rep(substring(d,1,f), no.repeats), collapse='')
    if (d==rep.string) {
      rep.unit <- substring(d,1,f)
    }
  }
  rep.unit
}
#' get all factors of an integer
get_all_factors <- function(n)
{
  prime_factor_tables <- lapply(
    setNames(n, n), 
    function(i)
    {
      if(i == 1) return(data.frame(x = 1L, freq = 1L))
      plyr::count(as.integer(gmp::factorize(i)))
    }
  )
  lapply(
    prime_factor_tables, 
    function(pft)
    {
      powers <- plyr::alply(pft, 1, function(row) row$x ^ seq.int(0L, row$freq))
      power_grid <- do.call(expand.grid, powers)
      sort(unique(apply(power_grid, 1, prod)))
    }
  )
}

finalcaller <- function(mhcount, replength, rept) {
  # mhcount
  # replength
  # rept
  
  
  if (replength >= mhcount) {
    if ((replength/nchar(rept)) >= 1) {
      return("Repeat-mediated")
    } else if (mhcount > 0) {
      return("Microhomology-mediated")
    } else {
      return("None")
    }
  } else {
    if (mhcount > 0) {
      return("Microhomology-mediated")
    } else { 
      return("None")
    }
  }
  
}

#' Generate 45 channel indel catalogue
#' 
#' @param muts_list A list of indels with flanking sequence information: Chrom, Pos, Ref, Alt, Sample column
#' @param Sample_col Sample column name
#' @param muttype_col indel classification column name
#' @return Indel catalogue

#' @export
gen_indelmuttype_MMRD <- function(muts_list, Sample_col, muttype_col){
  indel_template <- indelsig_template
  indel_template_uniq <- unique(indel_template[,c(muttype_col,"type")])
  names(indel_template_uniq) <- c("indelsubtype","type")
  indel_catalogue <- data.frame(table(muts_list[,Sample_col],muts_list[,muttype_col]))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- reshape2::dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template_uniq,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}
