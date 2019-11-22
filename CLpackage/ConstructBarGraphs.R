#72 rows per sample size in insurance
library(ggplot2)
library(lattice)
hailfinderCLResults <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\hailfinder_resultsCL.csv")
insuranceCLResults <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\insurance_resultsCL.csv")
asiaCLResults <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\asia_resultsCL.csv")
alarmCLResults <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\alarm_resultsCL.csv")
hailfinderMMHC_resultsCL <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\hailfinderMMHC_resultsCL.csv")
insuranceMMHC_resultsCL <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\insuranceMMHC_resultsCL.csv")
asiaMMHC_resultsCL <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\asiaMMHC_resultsCL.csv")
alarmMMHC_resultsCL <- read.csv("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\alarmMMHC_resultsCL.csv")
rownames(hailfinderCLResults) <- hailfinderCLResults$label
rownames(insuranceCLResults) <- insuranceCLResults$label
rownames(asiaCLResults) <- asiaCLResults$label
rownames(alarmCLResults) <- alarmCLResults$label
rownames(hailfinderMMHC_resultsCL) <- hailfinderMMHC_resultsCL$label
rownames(insuranceMMHC_resultsCL) <- insuranceMMHC_resultsCL$label
rownames(asiaMMHC_resultsCL) <- asiaMMHC_resultsCL$label
rownames(alarmMMHC_resultsCL) <- alarmMMHC_resultsCL$label

sizes <- c(500, 1000, 5000, 10000, 15000, 20000)
sizesNames <- c("500", "1000", "5000", "10000", "15000", "20000")
pvals <- c(0.05, 0.01, 0.005)
scores <- c("bde", "aic", "bic", "loglik")
starts <- c(3,4)
steps <- c(2,3,4)

vect <- data.frame("StepSize" = integer(32),
                   "SampleSize" = integer(32),
                   "SHD" = integer(32))



plotData <- function(pval,start,score)
{
  combos<-expand.grid(sizes,pvals,steps,starts,scores)
  
  first <- combos[combos$Var2 == pval & combos$Var4 == start & combos$Var5 == score,]
  firstAppend <- combos[combos$Var2 == pval & combos$Var4 == start & combos$Var5 == score,]
  
  colnames(first)[colnames(first)=="Var1"] <- "SampleSize"
  colnames(first)[colnames(first)=="Var2"] <- "Pvals"
  colnames(first)[colnames(first)=="Var3"] <- "StepSize"
  colnames(first)[colnames(first)=="Var4"] <- "Start"
  colnames(first)[colnames(first)=="Var5"] <- "Score"
  
  colnames(firstAppend)[colnames(firstAppend)=="Var1"] <- "SampleSize"
  colnames(firstAppend)[colnames(firstAppend)=="Var2"] <- "Pvals"
  colnames(firstAppend)[colnames(firstAppend)=="Var3"] <- "StepSize"
  colnames(firstAppend)[colnames(firstAppend)=="Var4"] <- "Start"
  colnames(firstAppend)[colnames(firstAppend)=="Var5"] <- "Score"
 
  first <- rbind(first,firstAppend) #print(merge(first, firstAppend))
  #first$SHD <- integer(32)
  print(first)
  
  
  rowMatchFirst <- paste("asia",first$Pvals,first$SampleSize,first$StepSize,first$Start,first$Score, sep = "_")
  print(rowMatchFirst)
  first$SHD <- asiaCLResults[rowMatchFirst,"SHD"]
  first$StepSize[19:36] <- "MMHC"
  first$SHD[19:36] <- asiaMMHC_resultsCL[rowMatchFirst, "SHD"]
  print(first)
  title = paste("P val:",pval,"Start size:", start, "Score:",score)
  
  file = paste0("Pval",pval,"Startsize", start, "Score",score)
  filename = paste0("C:\\Users\\TristansLaptop\\Documents\\SummerResearchMCNAIR2019\\SummerResearchMCNAIR2019\\clTests\\plots\\", file,".pdf")
  pdf(file=filename)
  #colnames(first)[colnames(first) == "MMHC"] <- "StepSize"
  library(ggplot2)
  
  # Faceting
  plot <- ggplot(first, aes(y=SHD, x = StepSize, fill = factor(StepSize))) +
    #aes(y=SHD, x = MMHC, fill = factor(MMHC)) + 
    geom_bar( stat="identity", width = 0.7, position = "dodge") +
    facet_wrap(~SampleSize, nrow = 1, strip.position = "bottom") +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          panel.border = element_rect(colour = "black", fill=NA),
          aspect.ratio = 2/1, axis.text = element_text(colour = 1, size = 13),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16),
          legend.position="bottom") +
    scale_fill_discrete(labels = c("Step size 2", "Step size 3", "Step size 4", "MMHC")) +
    coord_cartesian(ylim=c(0,65)) +
    scale_y_continuous(breaks=seq(0,65,10), expand = c(0, 0)) +
    xlab("Sample Size")
  
  
  print(plot)
  
  dev.off()
}

i<-1
for(pval in pvals){
  for(start in starts){
    for(score in scores){
      plotData(pval,start,score)
      i<-i+1
    }
  }
}
print(i)


