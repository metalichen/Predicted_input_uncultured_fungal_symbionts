library(dplyr)
library (stringr)
library (seqRFLP)
library(Biostrings)

#load the fun annotate table
genome_ann<-read.delim('/file/with/funannotate/annotation/tables',na.strings='')

#filter gene models labelled as secreted by SignalP
genome_signalp<-genome_ann[!is.na(genome_ann$Secreted),c(1,20)]

#write them as fasta
genome_signalp_faa<-dataframe2fas(genome_signalp,file='genome_signalp.faa')



#take the fasta file, analyze it with the TNHMM web server. The lines below process the output file produced by TMHMM
genome_tmhmm<-read.delim('genome_tmhmm.txt',sep='\t',header=F)
genome_tmhmm_readable<-genome_tmhmm %>% separate(V2, into=c('t1','len'),sep = "=",extra="merge") %>% separate(V3, into=c('t2','ExpAA'),sep = "=",extra="merge") %>%
  separate(V4, into=c('t3','First60'),sep = "=",extra="merge") %>%   separate(V5, into=c('t4','PredHel'),sep = "=",extra="merge") %>%
  separate(V6, into=c('t5','Topology'),sep = "=",extra="merge") %>% select(-c(t1,t2,t3,t4,t5))
write.table(genome_tmhmm_readable,"genome_tmhmm_readable.txt",sep='\t',quote = F,row.names = F)

genome_signalp2<-data.frame(genome_signalp$GeneID)
colnames(genome_signalp2)<-'V1'
genome_secreted<-left_join(genome_signalp2,genome_tmhmm_readable)
genome_secreted$PredHel<-as.numeric(genome_secreted$PredHel)
genome_secreted$len<-as.numeric(genome_secreted$len)
genome_secreted$ExpAA<-as.numeric(genome_secreted$ExpAA)
genome_secreted$First60<-as.numeric(genome_secreted$First60)

genome_secreted <- mutate(genome_secreted,delta=ExpAA-First60) 
genome_secreted <-genome_secreted[genome_secreted$delta<18 & genome_secreted$PredHel<2,]
write.table(genome_secreted,"genome_secreted_signalp_tmhmm.txt",sep='\t',quote = F,row.names = F)

genome_secreted_signalp_tmhmm_fa<-genome_ann[genome_ann$GeneID %in% genome_secreted$V1,c(1,20)]
genome_secreted_signalp_tmhmm_fa<-dataframe2fas(genome_secreted_signalp_tmhmm_fa,file='genome_secreted_signalp_tmhmm.faa')
genome_secreted_signalp_tmhmm_fa<-readAAStringSet('genome_secreted_signalp_tmhmm.faa')
writeXStringSet(genome_secreted_signalp_tmhmm_fa,'genome_secreted_signalp_tmhmm_batch1.faa')



#take the resulting file and analyze it with the WolfPSORT web server. The lines below process the output file produced by WolfPSORT
genome_wolf<-as.character(read.delim('genome_wolf_output.txt',header = F,as.is=T))
genome_wolf<-as.data.frame(unlist(strsplit(genome_wolf,'clgrpredictedgene_'))[-1])
colnames(genome_wolf)<-'V1'
genome_wolf2<-separate(genome_wolf,'V1',into=c('ID','pred'),extra='merge')
genome_wolf2$ID<-paste('clgrpredictedgene_',genome_wolf2$ID,sep='')

predictions<-c('plas', 'extr', 'cyto_nucl', 'nucl', 'cyto', 'pero', 'mito', 'E.R.', 'vacu', 'golg','cyto_mito','cysk')

for (p in predictions){
  regex<-paste(".*?",p,':([0-9]+).*$',sep='')
  vector<-str_extract(genome_wolf2$pred,regex)
  vector[!is.na(vector)]<-sub(regex, "\\1", vector[!is.na(vector)])
  genome_wolf2<-cbind(genome_wolf2,as.numeric(as.character(vector)))
  colnames(genome_wolf2)[ncol(genome_wolf2)]<-p}

genome_wolf2<- genome_wolf2 %>% mutate(X = extr/rowSums(.[3:14],na.rm = T))
genome_wolf2_filter<-genome_wolf2[genome_wolf2$X>0.6 & !is.na(genome_wolf2$X),1]
write.table(genome_wolf2_filter,"genome_secreted_signalp_tmhmm_wolf.txt",sep='\t',quote = F,row.names = F,col.names=F)



#write FINAL secretome fastas
genome_wolf2_fa<-genome_ann[genome_ann$GeneID %in% genome_wolf2_filter,c(1,20),]
genome_wolf2_fa<-dataframe2fas(genome_wolf2_fa,file='genome_secreted_signalp_tmhmm_wolf.faa')

#make a summary for the content of secretome
genome_wolf2_filter<-read.delim("genome_secreted_signalp_tmhmm_wolf.txt",header=F)[,1]

genome_secretome_summary<-data.frame()
genome_secretome_summary[1,1]<-length(genome_wolf2_filter)
genome_secretome_summary[1,2]<-length(genome_wolf2_filter)/nrow(genome_ann)
genome_secretome_summary[1,3]<-nrow(genome_ann[genome_ann$GeneID %in% genome_wolf2_filter & !is.na(genome_ann$CAZyme),])
genome_secretome_summary[1,4]<-nrow(genome_ann[genome_ann$GeneID %in% genome_wolf2_filter & !is.na(genome_ann$CAZyme),])/length(genome_wolf2_filter)
genome_secretome_summary[1,5]<-nrow(genome_ann[genome_ann$GeneID %in% genome_wolf2_filter & !is.na(genome_ann$Protease),])
genome_secretome_summary[1,6]<-nrow(genome_ann[genome_ann$GeneID %in% genome_wolf2_filter & !is.na(genome_ann$Protease),])/length(genome_wolf2_filter)
genome_effector<-read.delim('genome_effectorp.txt',sep='\t',header=T)
genome_secretome_summary[1,7]<- nrow(genome_effector[genome_effector$Prediction=='Effector',]) 
genome_secretome_summary[1,8]<- nrow(genome_ann[genome_ann$GeneID %in% genome_wolf2_filter & (genome_ann$Product!='hypothetical protein' | !is.na(genome_ann$PFAM)) & is.na(genome_ann$CAZyme) & is.na(genome_ann$Protease),])

colnames(genome_secretome_summary)<-c('secretome','genome_in_secretome','cazy','secretome_in_cazy','merops','secretome_in_merops','effectors','annotated_rest')
write.table(genome_secretome_summary,"genome_secretome_summary.txt",sep='\t',quote = F,row.names = F)
