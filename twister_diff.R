################################################
# TORNADO DIFF FUNCTION FOR SUBTRACTING INPUT #
################################################

tornado_diff <- function(gr,dataset,pad=250,ord=0,window=5,color="red2") {
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
      tidyr::spread(genotype,count,fill=0) %>% 
      dplyr::mutate(count=H3K36me2-Input) %>% 
      dplyr::mutate(genotype=factor("H3K36me2")) %>% 
      dplyr::select(genotype,which_label,startpos,count) %>% 
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



################################################
# Twister Diff Function for Subtracting Input #
################################################ 

twister_diff <- function(gr,dataset,pad=250,ord=1,window=5,color="red2",ya=c(0,20),sample_name="ChIP") {
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
      tidyr::spread(genotype,count,fill=0) %>% 
      dplyr::mutate(count=H3K36me2-Input) %>% 
      dplyr::mutate(genotype=factor("H3K36me2")) %>% 
      dplyr::select(genotype,startpos,whichgr,count) %>% 
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