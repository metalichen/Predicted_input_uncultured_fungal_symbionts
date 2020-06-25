#Script for generating GC vs coverage plots. Takes a fasta file, a binning file provided by CONCOCT, and a coverage file (see coverage.sh for how coverage files are generated)
library(dplyr)
library(tidyr)
library(plotly)
library(Biostrings)

###generating GC to coverage plot
#function to get a contig average coverage
getCov<-function(coverage_file){
  colnames(coverage_file) <- c("contig", "position", "coverage") #name columns
  cov.contig <- as.data.frame(tapply(coverage_file$coverage, coverage_file$contig, mean,na.rm=T)) #calculate avg coverage per contig
  colnames(cov.contig)<-c('coverage') #rename column
  return(cov.contig)
}

getGCcov<-function(fasta,coverage_file,...){
  #process fasta
  contig_names<-names(fasta) #get contig names
  nucl_freq<-alphabetFrequency(fasta)
  nucl_freq<-as.data.frame(nucl_freq) %>% mutate(gc=100*(G+C+S)/(A+T+W+G+C+S)) #get gc content
  dataset1<-data.frame(contig_names,nucl_freq$gc)
  row.names(dataset1)<-dataset1$contig_names
  #process coverage file
  dataset2<-getCov(coverage_file)
  #combine
  final_dataset<-merge(dataset1,dataset2,by = "row.names")
  return(final_dataset)
}

#read files  
cov<-read.delim('path/to/coverage/file',header=F)
fasta<-readDNAStringSet('path/to/fasta/file')
binning<-read.csv2('path/to/concoct/file',header=F)
colnames(binning)<-c('contig_names',"bin")
binning$bin<-paste('bin',binning$bin,sep='')

#get the data
gc_dataset1<-getGCcov(fasta,cov) #get GC and coverage data
gc_dataset2<-left_join(gc_dataset1[2:4],binning) #add info on the bin assignment
gc_dataset2$bin<-factor(gc_dataset2$bin,levels=c("bin0",'bin1',"bin2",'bin3',"bin4","bin5","bin6","bin7","bin8",
                                                 "bin10","bin11","bin12","bin13","bin14","bin15","bin16","bin17", "bin18","bin19","bin20","bin21","bin22",
                                                 "bin23","bin24","bin25","bin26","bin27",'bin9',"unbinned"))

write.table(gc_dataset2,"gc_plot.txt",sep='\t',quote = F,row.names = F)

#make a graph
gc_dataset2<-read.delim("gc_plot.txt",as.is =T)
colnames(gc_dataset2)<-c("contig","nucl_freq.gc","coverage","bin")
gc_dataset2[is.na(gc_dataset2$bin),4]<-'unbinned'
gc_dataset2$bin<-as.factor(gc_dataset2$bin)

palette<-c('#47306a','#6cd445','#a94de4','#cddc32','#592ea8','#64d57c','#c656c2','#aec04f','#6574dc','#dd9c2b','#908acf',
           '#a7b696','#d84472','#5fbb8c','#c74694','#779c4f','#652a42','#cab176','#6c9abe','#b25038','#a1a7b2','#b27c36','#c789b3','#534a25',
           '#b35055','#373e3c','#bb8a82','#db482b','#61c2bc')
gc_plot<-plot_ly(data=gc_dataset2,x=~nucl_freq.gc,y=~coverage,color=~bin,colors = palette,marker=list(size=4)) %>%add_markers()%>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))

htmlwidgets::saveWidget(as.widget(gc_plot), "gc_plot.html", selfcontained = FALSE) 

#calculate median coverage per bin
coverage_table<-gc_dataset2 %>% group_by(bin) %>% summarise(median=median(coverage),medianGC=median(nucl_freq.gc))
write.table(coverage_table,"bin_coverage.txt",sep='\t',quote = F,row.names = F,col.names=T)

