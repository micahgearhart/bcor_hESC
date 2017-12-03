#this was from vince buffalo's book
#swl<-function (x) {
#  sapply(strsplit(sub("([\\d+MYX]+):(\\d+)-(\\d+)", "\\1;;\\2;;\\3",x, perl=TRUE),";;"),"[[",2)
#}

swl<- function(x) {as.numeric(sapply(strsplit(as.character(x),":|-"),function(x) x[2]))}
wls<- function(x) { paste0(seqnames(x),":",start(x),"-",end(x))}


strstrip<- function(s,i) sapply(strsplit(s,"_"),function(x) x[i])

label_wrap <- function(variable, value) {
  lapply(strwrap(as.character(value), width=25, simplify=FALSE),
         paste, collapse="\n")
}

#http://stackoverflow.com/questions/17319487/median-and-quartile-on-violin-plots-in-ggplot2
median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out)
}


narrowPeakToGRanges<-function(file) {
  x <- read.table(file,stringsAsFactors=F)
  gr <-GRanges(seqnames=x$V1, ranges = IRanges(start=x$V2, end=x$V3),
               strand="*", score=x$V5, e=x$V7,summit=x$V10)
  return(gr)
}

broadPeakToGRanges<-function(file,name="peak") {
  x <- read.table(file,stringsAsFactors=F)
  gr <-GRanges(seqnames=x$V1, ranges = IRanges(start=x$V2, end=x$V3),
               strand="*", score=x$V5, e=x$V7,p=x$V8,q=x$V9)
  names(gr)<-paste0(name,"_",formatC(1:length(gr),width=5,format="d",flag="0"))
  return(gr)
}

center<-function(gr) {
  GRanges(seqnames = seqnames(gr),
          IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                         (width(gr)+1)/2),width=1))
}

setClass("fileset", representation( filename="character",count="numeric",labels="character"))

countFileset<-function (x) {
  for (i in 1:length(x@filename)) {
    x@count[i]<-countBam(x@filename[i],param=ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE),))$records
  }
  return(x)
}

reds <- colorRampPalette(brewer.pal(n = 9, "Reds"))

p_param <- PileupParam(distinguish_nucleotides=FALSE,
                       distinguish_strands=FALSE,
                       min_nucleotide_depth=0)

#plotCounts function for DESeq DDS objects
gg_plotCounts<-function(x="ENSG00000163508",d=dds) {
  if (substr(x,1,4)=="ENSG") {
    title<-unique(hgnc[grep(x,hgnc$ensembl_gene_id),"hgnc_symbol"])
  } else {
    title<-x
    x<-unique(hgnc[grep(paste0("^",title,"$"),hgnc$hgnc_symbol),"ensembl_gene_id"])
  }
  
  plotCounts(d,x,intgroup=c("sample","time"),returnData=T) %>%
    mutate(l2=log2(count + 0.5)) %>% 
    ggplot(aes(x=time, y=l2)) + facet_wrap(~sample)+
    geom_point(position=position_jitter(w=0.1,h=0)) + ggtitle(paste0(x,":  ",title)) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, 
                 geom = "crossbar", width = 0.35,color="blue",size=0.2) +
    expand_limits(x=0, y = 0) + ylab("log2(normalized counts)")+ xlab("")+
    theme_bw()
  
}


pileup_by_pos <- function(i,env = parent.frame()) {
  style<-seqlevelsStyle(seqinfo(BamFile(dataset@filename[which(files==i)])))[1]
  
  pileup(dataset@filename[which(files==i)], 
         scanBamParam=sb_param[[style]], 
         pileupParam=p_param) %>% 
    mutate(which_label=gsub("chr","",which_label)) %>%
    mutate(seqnames=gsub("chr","",seqnames)) %>%
    #mutate(count=ceiling(count*1e6/dataset@count[which(files==i)])) %>% 
    mutate(genotype=dataset@labels[which(files==i)]) %>% 
    mutate(startpos=0.001*plyr::round_any(pos-swl(which_label)-pad - 0.5*width(gr[1]),window)) %>% 
    mutate(startpos=ifelse(strandHash[which_label]=="-",-1*startpos,startpos)) %>% 
    group_by(genotype,which_label,startpos) %>%
    #summarize(count=sum(count)) %>% 
    summarize(count=sum(count)*1e6/(dataset@count[which(files==i)]*window)) %>% 
#    mutate(count=log2(count)) %>% 
    ungroup() 
}






####################
# TORNADO FUNCTION #
####################

tornado <- function(gr,dataset,pad=250,ord=0,window=5,color="red2") {
  environment(pileup_by_pos) <- environment()
  
  #if gr is not uniform width, find the center of the gr and set width to 1
  if (length(gr) > 1 & var(width(gr))!=0) {
    gr<-GRanges(seqnames = seqnames(gr),
                IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                               (width(gr)+1)/2),width=1))
  }
  
  #add which_label to gr
  gr2<-gr+pad
  gr2$which_label<-paste0(seqnames(gr2),":",start(gr2),"-",end(gr2))
  #remove duplicated regions
  gr2<-gr2[!duplicated(gr2$which_label)]
  
  #set the color to whatever the user wants.
  reds<-colorRampPalette(c('snow', color))
  
  #this tells Rsamtools what and where we need data out of the file
  #sb_param = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
  #                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  sb_param<-list()
  seqlevelsStyle(gr2)<-"NCBI"
  sb_param[["NCBI"]] = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                                    flag = scanBamFlag(isUnmappedQuery = FALSE))
  seqlevelsStyle(gr2)<-"UCSC"
  sb_param[["UCSC"]] = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                                    flag = scanBamFlag(isUnmappedQuery = FALSE))
  
  #scan the bam files
  
  files <- gsub("\\.bam","",sapply(strsplit(dataset@filename,"\\/"), function (x) x[length(x)]))
  names(files)<-dataset@labels
  
  #plot labeller
  #plot_labeller <- function(variable,value){ return(names(which(files==value))) }
  

  
  #print(paste0("Bam file read in:  ",t1))
  
  #swl<- function(x) {as.numeric(strsplit(as.character(x),":|-")[[1]][2])}
  
  #if not ordering use the first file to do the ordering
  filter_file <- ifelse(ord==0,files[1],files[ord])
  
  
  which_factor_levels<-rev(as.character(gr2$which_label))
  
  if (ord==0) {
    spHash<-swl(which_factor_levels)
    names(spHash)<-which_factor_levels
    
  } else {
    #print("ORDER IS NOT 0")
    temp2 <- resL %>%
      bind_rows(.id="genotype") %>%
      dplyr::filter(genotype==filter_file) %>%
      group_by(which_label) %>%
      summarize(s = sum(count)) %>%
      arrange(s) %>%
      mutate(which_label = as.character(which_label)) %>%
      dplyr::select(which_label)
    
    temp1<-which_factor_levels[!which_factor_levels %in% temp2$which_label]
    which_factor_levels<-c(temp1,temp2$which_label)
    
    spHash<-swl(which_factor_levels)
    names(spHash)<-which_factor_levels
    
  }
  
  #strandHash
  strandHash<-as.character(strand(gr2[match(which_factor_levels,gr2$which_label)]))
  names(strandHash)<-which_factor_levels
  
  
  suppressWarnings(  
                     
                     bplapply(files,pileup_by_pos) %>% 
                       bind_rows() %>% 
                       mutate(genotype=factor(genotype,levels=names(files))) %>%
                       mutate(which_label = factor(which_label,levels=names(spHash))) %>% 
                       mutate(count=log2(count)) %>% 
                       group_by(genotype) %>%   #scale across genotype
                        dplyr::mutate(count=(count-min(count))/(max(count)-min(count))) %>%
                       ungroup() %>%
                       ggplot(aes(x=startpos,y=which_label,fill=count)) +
                       geom_raster() +
                       scale_fill_gradient2(
                         #low="white",mid=reds(256)[16],high=color,midpoint=0.0,
                         #low = reds(256)[1],
                                          low="white",
                                          mid = reds(256)[64],
                                           high = reds(256)[256],
                                            midpoint = 0.15, space = "Lab",
                         name = "Normalized Log2(cpm)",limits=c(0, 1)) +
                       
                       ylab("Genomic Regions") + xlab("Position Relative to Peak Center (kb)") +
                       facet_grid(. ~ genotype,scales="free_y",labeller=label_wrap) +
                       theme_bw() + theme(strip.background = element_blank(),
                                          axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank(),
                                          panel.grid.major=element_blank(),
                                          panel.grid.minor=element_blank())
  )
  
  
}





######################
# Twister Function #
######################

twister <- function(gr,dataset,pad=250,ord=1,window=5,color="red2",ya=c(0,20),sample_name="ChIP") {
  environment(pileup_by_pos) <- environment()
  
  splitbygr<-FALSE
  
  #if gr is a GRangeslist, list names kept in names()
  if( class(gr)=="GRangesList") {
    splitbygr<-TRUE
    normalizerHash=unlist(lapply(gr,length))
    gr<-unlist(gr)
  } else if (class(gr)!="GRanges") { 
    stop("The Range must be a GRanges object or a GRangesList")
  } 
  
  
  #if gr is not uniform width, find the center of the gr and set width to 1
  if (length(gr) > 1 & var(width(gr))!=0) {
    gr<-GRanges(seqnames = seqnames(gr),
                IRanges(start=start(gr)+ifelse(width(gr)%%2==0,width(gr)/2,
                                               (width(gr)+1)/2),width=1))
  }
  
  #add which_label to gr
  gr2<-gr+pad
  gr2$which_label<-paste0(seqnames(gr2),":",start(gr2),"-",end(gr2))
  
  
  #set the color to whatever the user wants.
  #reds<-colorRampPalette(c('snow', color))
  
  #this tells Rsamtools what and where we need data out of the file
  sb_param<-list()
  seqlevelsStyle(gr2)<-"NCBI"
  sb_param[["NCBI"]] = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                          flag = scanBamFlag(isUnmappedQuery = FALSE))
  seqlevelsStyle(gr2)<-"UCSC"
  sb_param[["UCSC"]] = ScanBamParam(what = c("pos", "rname", "strand","qwidth"), which = gr2,
                                    flag = scanBamFlag(isUnmappedQuery = FALSE))
  
  #scan the bam files
  
  files <- gsub("\\.bam","",sapply(strsplit(dataset@filename,"\\/"), function (x) x[length(x)]))
  names(files)<-dataset@labels
  
  
  #if not ordering use the first file to do the ordering
  filter_file <- ifelse(ord==0,files[1],files[ord])
  
  
  which_factor_levels<-rev(as.character(gr2$which_label))
  spHash<-swl(which_factor_levels)
  names(spHash)<-which_factor_levels
  
  #strandHash
  strandHash<-as.character(strand(gr2[match(which_factor_levels,gr2$which_label)]))
  names(strandHash)<-which_factor_levels
  
  #grHash
  if(splitbygr) {
    grHash<-sapply(strsplit(names(gr2),"\\."), function(x) x[1])
    names(grHash)<-gr2$which_label
  } else {
    grHash<-rep("test1",length(gr2))
    names(grHash)<-gr2$which_label
    normalizerHash=length(gr2)
    names(normalizerHash)="test1"
  }
  

  suppressWarnings( 
    
    bplapply(files,pileup_by_pos) %>% 
      bind_rows() %>% 
      mutate(genotype=factor(genotype,levels=names(files))) %>%
      #mutate(which_label = factor(which_label,levels=names(spHash))) %>% 
      mutate(whichgr=factor(grHash[which_label],levels=names(normalizerHash))) %>% 
      dplyr::select(genotype,startpos,count,whichgr) %>% 
      group_by(genotype,startpos,whichgr) %>%   #sum across which_labels
      summarize(count=sum(count)) %>%
      ungroup() %>%
      dplyr::filter((startpos != max(startpos)) & (startpos != min(startpos))) %>% #edge effects
      dplyr::mutate(normalizer = normalizerHash[as.character(whichgr)]) %>% 
      dplyr::mutate(count = count/normalizer) %>% 
     # dplyr::mutate(count = log2(count)) %>% 
     # tidyr::spread(genotype,count) %>% 
    #  dplyr::mutate(count = CHIP - INPUT) %>% 
     # dplyr::mutate(genotype= sample_name) %>% 
      group_by(genotype) %>%   #scale across genotype
      #dplyr::mutate(count = scale(count,center=TRUE,scale=FALSE)) %>% 
      #dplyr::mutate(count = count - min(count)) %>% 
      dplyr::mutate(count=(count-min(count))/(max(count)-min(count))) %>%
      ungroup() %>%
      ggplot(aes(x=startpos,y=count)) +
      geom_line(aes(color=whichgr)) +
      #stat_smooth(aes(color=whichgr),method="loess") +
      ylab("Normalized Counts per gene") + xlab("Position Relative to TSS (kb)") +
      facet_grid(.~ genotype) +
     # scale_color_manual(values=c("#EF3A38","#BFBFBF","#A2AED9"),name="")+
      scale_color_manual(values=color,name="")+
      #ylim(ya) +
      #  facet_grid(. ~ genotype,scales="free_y",labeller=label_wrap) +
      theme_bw() + theme(panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank())
      
  )
  

  }



