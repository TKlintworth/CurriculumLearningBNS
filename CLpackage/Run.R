Run <- function()
{
  source('~/SummerResearchMCNAIR2019/SummerResearchMCNAIR2019/CL.R')
  sizes <- c(500, 1000, 5000)#, 10000, 15000, 20000)
  #sizes <- c(500, 1000)
  pvals <- c(0.05, 0.01, 0.005)
  scores <- c("bde", "aic", "bic", "loglik")
  starts <- c(3,4)
  steps <- c(2,3,4)
  q <- 1
  #How many rows will this data frame have? length(sizes)*length(pvals)*length(scores)*length(starts)*length(steps)
  rows = length(sizes) * length(pvals) * length(scores) * length(starts) * length(steps)
  df <- data.frame(label = c(rows),
                   TP = integer(rows),
                   FP = integer(rows),
                   FN = integer(rows),
                   SHD = integer(rows))
                   #MMHC = integer(rows))
  # pre-allocate space
  for(i in 1:rows){
    df$label[i] <- "NA"
    df$TP[i] <- 0
    df$FP[i] <- 0
    df$FN[i] <- 0
    df$SHD[i] <- 0
    #df$MMHC[i] <- 0
  }
  print(df)
  #print(str(df))

  
  # load the data.
  data(asia)


  true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

  
  for(size in sizes){
    for(pval in pvals){
      for(score in scores){
        for(step in steps){
          for(start in starts){
            learned <- mmhc(asia[sample(nrow(asia), size), ])#mmhc(alarm)#recoverWholeBN(data = alarm, initSize = start, stepSize = step, size = size, pval = pval, score = score)
            diff <- compare(true, learned)
            SHD <- shd(learned, true)
            df[q, "label"] <- as.character(paste("asia",pval,size,step,start,score,sep = "_"))
            df[q, "TP"] <- diff$tp
            df[q, "FP"] <- diff$fp
            df[q, "FN"] <- diff$fn
            df[q, "SHD"] <- SHD
            #df[q, "MMHC"] <- shd(learned, mmhc(true))#shd(mmhc(as.data.frame(learned)),true)
            q <- q + 1
            print(df)
            #print(paste("asia",pval,size,step,start,score,sep = "_"))
          }
        }
      }
    }
  }
  write.csv(df, file = "asiaMMHC_resultsCL.csv", row.names = FALSE)
  #print("shd(learned, true):")
  #print(shd(learned,true))
  
  
}
