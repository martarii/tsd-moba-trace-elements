library(data.table)
library(foreign)
library(stringr)
library(dplyr)

#Output file path/name
outname <- ""

#Input files
PCfile <- ""
INFOfile <- ""
KeyFile <- ""
TraceElementFile <- ""

#Read and clean input

PCs <- fread(PCfile) #Principal components

INFO <- read.spss(INFOfile.to.data.frame=T) #Connects pregnancy IDs and unique mother IDs, so that we can identify mothers who participated or have been genotyped multiple times
INFO$M_ID_2944 <- str_replace_all(string=INFO$M_ID_2944,pattern="\t",repl="") #Remove whitespace
INFO$M_ID_2944 <- str_replace_all(string=INFO$M_ID_2944,pattern=" ",repl="") #Remove whitespace

TE <- read.spss(TraceElementFile,to.data.frame=T)
trace_elements = trace_elements[!(trace_elements$PREG_ID_2944==41735|trace_elements$PREG_ID_2944==87746)] #Remove preg ids without trace element measurement

#Order by age and remove duplicates, keeping the measurements from the latest pregnancies.
trace_elements=trace_elements[order(trace_elements[,"M_ID_2944"],-trace_elements[,"AgeSample"],]
trace_elements = trace_elements[!duplicated(trace_elements$M_ID_2944),]

#MoBa keyfile
MoBa_mother <- read.spss(KeyFile,to.data.frame=T)

MoBa_mother$SENTRIX_ID <- str_replace_all(string=MoBa_mother$SENTRIX_ID,pattern=" ",repl="")
MoBa_mother$SENTRIX_ID <- str_replace_all(string=MoBa_mother$SENTRIX_ID,pattern=\t",repl="")
MoBa_mother$M_ID_2944 <- str_replace_all(string=MoBa_mother$M_ID_2944,pattern=" ",repl="")
MoBa_mother$M_ID_2944 <- str_replace_all(string=MoBa_mother$M_ID_2944,pattern=\t",repl="")
MoBa_mother$BATCH <- str_replace_all(string=MoBa_mother$BATCH,pattern=" ",repl="")
MoBa_mother$BATCH <- str_replace_all(string=MoBa_mother$BATCH,pattern=\t",repl="")

#Select HARVEST batch (7 individuals from other batches are excluded)
MoBa <- MoBa_mother[MoBa_mother$BATCH=="HARVEST",]

#Combine data
tmp1 <- merge(PCs,MoBa, by.x="FID",by.y="SENTRIX_ID",all=F)
tmp1 <- tmp1[!duplicated(tmp1$M_ID_2944),]
tmp2 <- merge(trace_elements,tmp1,by="M_ID_2944",all=F)

tmp2$MATID <- 0
tmp2$PATID <- 0

#Final output
pheno <- tmp2[,c("FID","IID","MATID","PATID","Mn_Mangan","Co_Kobolt","Cu_Kopper","Zn_Sink","As_Arsen","Se_Selen","Mo_Molybdeno","Cd_Kadmium","Tl_Tallium","Pb_Bly","Hg_total","AgeSample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

write.table(pheno,outname,quote=F,col.names=T,row.names=F)
