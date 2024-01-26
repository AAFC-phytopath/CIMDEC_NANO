summ <- read.csv("/mnt/DATA_2/Biovigilance/20231108_1208_MC-115710_FAX57101_dd2339fd/sequencing_summary_FAX57101_dd2339fd_920c1ccf.txt", sep="\t" )
setwd("/mnt/DATA_2/Biovigilance/Tech_compare/summary/")

head(summ)
###################################
# compare barcodes 25 and 27 

bc25 <- subset(summ, barcode_arrangement=="barcode25")
head(bc25)

bc27 <- subset(summ, barcode_arrangement=="barcode27")
head(bc27)

str(bc25)

library(ggplot2)

plot_bc25<-ggplot(bc25, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc25")

plot_bc27<-ggplot(bc27, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc27")

###############################
# compare barcodes 26 and 28 

bc26 <- subset(summ, barcode_arrangement=="barcode26")
head(bc26)

bc28 <- subset(summ, barcode_arrangement=="barcode28")
head(bc28)

library(ggplot2)

plot_bc26<-ggplot(bc26, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc26")

plot_bc28<-ggplot(bc28, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc28")

library(gridExtra)

grid.arrange(plot_bc25, plot_bc27, plot_bc26, plot_bc28,
             nrow = 2, ncol = 2)

###################################
# compare barcodes 25 and 27 

plot_bc25b<-ggplot(bc25, aes(x=mean_qscore_template)) + 
  geom_boxplot()+
  xlim(0,25)+
  ggtitle("bc25")

plot_bc27b<-ggplot(bc27, aes(x=mean_qscore_template)) + 
  geom_boxplot()+
  xlim(0,25)+
  ggtitle("bc27")

###############################
# compare barcodes 26 and 28 

plot_bc26b<-ggplot(bc26, aes(x=mean_qscore_template)) + 
  geom_boxplot()+
  xlim(0,25)+
  ggtitle("bc26")

plot_bc28b<-ggplot(bc28, aes(x=mean_qscore_template)) + 
  geom_boxplot()+
  xlim(0,25)+
  ggtitle("bc28")

library(gridExtra)

grid.arrange(plot_bc25b, plot_bc27b, plot_bc26b, plot_bc28b,
             nrow = 2, ncol = 2)




###################################
# compare barcodes 31 and 32 

bc31 <- subset(summ, barcode_arrangement=="barcode31")
head(bc31)

bc32 <- subset(summ, barcode_arrangement=="barcode32")
head(bc32)

library(ggplot2)

plot_bc31<-ggplot(bc31, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc31")

plot_bc32<-ggplot(bc32, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc32")

grid.arrange(plot_bc31, plot_bc32,
             nrow = 2)

###############################
# compare barcodes 26 and 28 

bc26 <- subset(summ, barcode_arrangement=="barcode26")
head(bc26)

bc28 <- subset(summ, barcode_arrangement=="barcode28")
head(bc28)

library(ggplot2)

plot_bc26<-ggplot(bc26, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc26")

plot_bc28<-ggplot(bc28, aes(x=sequence_length_template)) + 
  geom_histogram(binwidth=30)+
  xlim(0,2500)+
  ggtitle("bc28")

library(gridExtra)

grid.arrange(plot_bc25, plot_bc27, plot_bc26, plot_bc28,
             nrow = 2, ncol = 2)
