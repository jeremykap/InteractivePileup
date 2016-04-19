library(Rsamtools)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(plyr)
library(plotly)

prepareBedForPileup <- function(bed) {
  bedTable = read.table(bed)
  names(bedTable) = c("chr","start","stop","amp_name","gc","numBase","unkn","gene")
  return (bedTable)
}

piPlot <- function(pPlot) {
  toPlot <- pPlot[['toPlot']]
  bedSubset <- pPlot[['bedSubset']]
  location <- pPlot[['location']]
  refFrame <- pPlot[['refFrame']]
  cols = c("A"="green", "C"="blue", "G"="brown","T"="red","-"="magenta","+"="purple","Reference"="grey")
  bedCols = c("black","blue","green")
  piPlot <- plot_ly(data=toPlot[toPlot$strand=="-"& toPlot$panel=="pileup",],x=pos,y=y,marker=list(color=cols[color]),type="bar",hoverinfo="all",showlegend=F,fixedrange=TRUE,legendgroup="group",name="-",xaxis="x1",yaxis="y2",text=text)
  piPlot <- add_trace(p=piPlot,data=toPlot[toPlot$strand=="+"& toPlot$panel=="pileup",],type="bar",x=pos,y=y,marker=list(color=cols[color]),hoverinfo="all",fixedrange=TRUE,legendgroup="group", showlegend=FALSE, name="+",xaxis="x1",yaxis="y1",text=text) 
  piPlot <- add_trace(p=piPlot,data=toPlot[toPlot$panel=="errors",],x=pos,y=y,marker=list(color=cols[color]),hoverinfo="all",type="bar",showlegend=FALSE,legendgroup="group",name="Errors",xaxis="x1",yaxis="y3",text=text)
  piPlot <- add_trace(p=piPlot,data=refFrame,x=pos,y=y,text=base,mode="text", textfont=list(color=cols[refFrame$base]),type="scatter",xaxis="x1",yaxis="y2",hoverinfo="none",showlegend=FALSE,legendgroup="group",name="REF")
  if(!is.null(bedSubset)) {
    piPlot <- add_trace(p=pPlot,data=bedSubset,x=start,y=y,showlegend=F,type="scatter",mode="lines",hoverinfo="text",xaxis="x1", yaxis="y4",legendgroup="group",text=text,color=bedCols[count],colors="Accent")
  }
  layout(piPlot, barmode="stack",
         shapes = list(list(type="rect",fillcolor="orange",opacity=.2,x0=location-.5,x1=location+.5,xref="x1",y0=min(toPlot$y)-50,y1=max(toPlot$y)+50,yref="y2"),
                       list(type="rect",fillcolor="orange",opacity=.2,x0=location-.5,x1=location+.5,xref="x1",y0=min(toPlot$y)-50,y1=max(toPlot$y)+50,yref="y1")),
         yaxis = list(anchor = 'x',title="Forward",range=c(0,max(toPlot$y)+100),domain = c(.50, .90),fixedrange=TRUE), 
         yaxis2 = list(anchor = 'x',title="Reverse",range=c(min(toPlot$y)-100,0),domain=c(.1,.50),fixedrange=TRUE), 
         yaxis3=list(anchor="x", range=c(0,1), domain=c(0,0.080),fixedrange=TRUE,showticklabels=FALSE),
         yaxis4=list(anchor="x", range=c(0,1), domain=c(.95,1),fixedrange=TRUE,showticklabels=FALSE),
         xaxis = list(title="",exponentformat="none",hoverformat="g",type="linear",tickformat="g"),title="Pileup",hovermode="closest")
  
}
preparePileupData <- function(bamFile,chrom,location,bedTable=NULL,errorDat=NULL) {
  
  start = location-50
  end = location+50
  ref = getSeq(Hsapiens,names=chrom,start,end)
  chr = tolower(chrom)
  if (chr=="chrx") {
    chr <- "chrX"
  } else if (chr=="chry") {
    chr <- "chrY"
  }
  chrom <- chr
  scan_params = ScanBamParam(which=GRanges(seqnames = chrom,IRanges(start = start, end = end,names = chrom)))
  pileup_params = PileupParam(min_base_quality =0, min_map=0,max_depth=10000,min_nucleotide_depth = 0,include_insertions=T)
  p = pileup(file=bamFile,scanBamParam = scan_params,pileupParam = pileup_params)
  p$mismatch = logical(nrow(p))
  p$plotX = numeric(nrow(p))
  p$newPos = p$pos - start +1
  p$strandMult = ifelse(p$strand=="+",1,-1)
  p$y = p$count * p$strandMult
  p$color = factor(character(nrow(p)),levels=c("Reference","-","+","A","C","G","T")) 
  for (i in 1:nrow(p)) {
    if ( p[i,]$nucleotide != as.character(ref[ p[i,]$newPos])){
      p[i,]$mismatch =T
      p[i,]$color = as.character(p[i,]$nucleotide)
    } else{
      # p[i,]$y = 0 
      p[i,]$color = "Reference"
    }
  }
  
  mismatches <- p
  errorSubset = data.frame()
  if (!is.null(errorDat)) {
    errorSubset = errorDat[errorDat$chr==chrom & errorDat$pos>=start & errorDat$pos<=end,]
  }
  bedSubset=NULL
  if (!is.null(bedTable)){
    bedSubset = bedTable[(bedTable$chr==chrom & bedTable$start>=start & bedTable$stop <= end),]

    bedSubset = ddply(.data = bedSubset,~start,summarise,text=paste(amp_name,collapse=", "),count=length(amp_name))
    bedSubset$y = .5
    
  } 
  
  mismatches = mismatches[order(mismatches$color,decreasing = T),]
  toPlot = data.frame(pos=mismatches$pos,y=mismatches$y,color=mismatches$color,panel=factor("pileup",levels=c("bed","pileup","errors")),strand=factor(mismatches$strand,levels=c("+","-","*")))
  if (nrow(errorSubset)>0) {
    for (l in errorSubset$pos) {
      current = errorSubset[errorSubset$pos==l,]
      errorSubset[nrow(errorSubset)+1,] = c(chrom,l,current[1,]$gene,ref=current[1,]$ref,num_samps=current[1,]$num_samps,type2=current[1,]$ref,noise=as.numeric((1-sum(as.numeric(current$noise)))))
    }
    errorSubset$noise=as.numeric(errorSubset$noise)
    errorSubset$pos = as.numeric(errorSubset$pos)
    errorSubset$type2 = ifelse(errorSubset$type2=="Ins","+",ifelse(errorSubset$type2=="Del","-",errorSubset$type2))
    
    errorSubset$panel = "errors"
    errorSubset$noise = errorSubset$noise
    errorSubset$strand = "+"
    for (row in 1:nrow(errorSubset)) {
      toPlot[nrow(toPlot)+1,] = c(pos=as.numeric(errorSubset[row,]$pos),y=as.numeric(as.numeric(errorSubset[row,]$noise)),color=errorSubset[row,]$type2,panel=errorSubset[row,]$panel,strand=errorSubset[row,]$strand)
    }
  }
  toPlot$pos= as.numeric(toPlot$pos)
  toPlot$y = as.numeric(toPlot$y)
  toPlot$color = as.character(toPlot$color)
  refFrame = data.frame(pos=seq(start,end,1),base=ref,panel=factor("pileup",levels=c("bed","pileup","errors")))
  row.names(refFrame) = seq(start,end,1)
  refFrame$y = min(toPlot$y)-25
  toPlot$text = ifelse(toPlot$color=="Reference",
                       paste("Ref:",refFrame[as.character(toPlot$pos),]$base),
                       paste("Ref:",refFrame[as.character(toPlot$pos),]$base,"\n ALT:",toPlot$color))
  cols = c("T"="red","G"="brown","C"="blue","A"="green","-"="magenta","+"="purple","Reference"="grey")
  bedCols = c("black","blue","green")
  cols = c("T"="red","G"="brown","C"="blue","A"="green","-"="magenta","+"="purple","Reference"="grey")
  bedCols = c("black","blue","green")
  
  
  return(list(location=location,toPlot=toPlot,bedSubset=bedSubset,refFrame=refFrame))
}