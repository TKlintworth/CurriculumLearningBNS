
#######################
# File: CL.R
# AUthor: Tristan Klintworth
# Purpose: Uses functions from the bnlearn R package (http://www.bnlearn.com/) and infotheo R package (https://cran.r-project.org/web/packages/infotheo/index.html) 
#   and reproduces **without implementing the penalty function**
# the algorithm from the paper [insert paper citation here]
# Functions:
#   -getDataset():
#   -ConstructCurriculum():
#######################


#######################
# Author: Tristan Klintworth
# Code adapted from pseudocode in Learning Bayesian network structures under incremental construction curricula (Yanpeng Zhao, Yetian Chen, Kewei Tu, Jin Tian)
# Parameter numeric: FALSE by default gives names of the variables, passing TRUE gives numeric indices
#######################

#Entering initSize as 1 causes an error
ConstructCurriculum <- function(dataset, initSize, stepSize = 4, numeric=FALSE)
{
  if(stepSize %% 1 != 0)
  {
    print("Please enter a positive integer step size")
  }
  else
  {
    ##if numeric argument passed in as true, give numeric indices of curriculum
    if(numeric == TRUE)
    {
      ### Initialization 
      
      V <- colnames(dataset)
      MI <- mutinformation(dataset)
      diag(MI) <- 0
      aveMI <- (rowMeans(MI))
      xStar <- which.max(aveMI)
      X1 <- as.vector(names(xStar))
      
      for (i in 1:(initSize-1))
      {
        Ynames <- colnames(MI[ , !colnames(MI) %in% X1])
        Ynums <- (MI[ , !colnames(MI) %in% X1])
        
        aveMis <- numeric(length = length(Ynames))
        aveMisNames <- character(length = length(Ynames))
        
        for (y in Ynames)
        {
          ySums = 0
          for (x in X1)
          {
            ySums = ySums + MI[y,x]
          }
          aveMis[i] <- ySums/length(X1)
          aveMisNames[i] <- y
          i <- i + 1
        }
        Ystar <- aveMisNames[which.max(aveMis)]
        X1 <- c(X1, Ystar)
      }
      #end
      
      #m <- ceiling((n-s)/t + 1)
      m <- ceiling(((length(V)-initSize)/stepSize) + 1)
      Ynames <- colnames(MI[ , !colnames(MI) %in% X1])
      
      #this initializes curriculum to a list containing m vectors
      curriculum <- vector("list",m)
      curriculum[[1]] <- X1
      
      for (i in 2:m)
      {
        curriculum[[i]] <- curriculum[[i-1]]
        
        for(j in 1:stepSize)
        {
          #will only evaluate if Y is non empty
          
          if(length(Ynames))
          {
            aveMis <- c()#numeric(length = length(Ynames))
            aveMisNames <- c()#character(length = length(Ynames))
            currentX <- c()
            
            for (y in Ynames)
            {
              k <- 1
              ySums = 0
              for (x in curriculum[[i]])
              {
                ySums = ySums + MI[y,x]
              }
              
              #aveMis[k] <- ySums/length(curriculum[[i]])
              aveMis <- c(aveMis, (ySums/length(curriculum[[i]])))
              aveMisNames <- c(aveMisNames, y)
              k <- k + 1
            }
            #Ystar <- which(V==aveMisNames[which.max(aveMis)])
            Ystar <- aveMisNames[which.max(aveMis)]
            Ynames <- Ynames[!(Ynames %in% Ystar)]
            
            currentX = c(currentX, Ystar)
            #curriculum[[i]] <- currentX
            curriculum[[i]] <- c(curriculum[[i]], currentX)
            j <- j + 1
          }
        }
      }
      for(i in 1:length(curriculum))
      {
        curriculum[[i]] <- (match(curriculum[[i]], colnames(dataset)))-1
        
      }
      return(curriculum)
    }
    else
    {
      ### Initialization 
      
      V <- colnames(dataset)
      MI <- mutinformation(dataset)
      diag(MI) <- 0
      aveMI <- (rowMeans(MI))
      xStar <- which.max(aveMI)
      X1 <- as.vector(names(xStar))
      
      for (i in 1:(initSize-1))
      {
        Ynames <- colnames(MI[ , !colnames(MI) %in% X1])
        Ynums <- (MI[ , !colnames(MI) %in% X1])
        
        aveMis <- numeric(length = length(Ynames))
        aveMisNames <- character(length = length(Ynames))
        
        for (y in Ynames)
        {
          ySums = 0
          for (x in X1)
          {
            ySums = ySums + MI[y,x]
          }
          aveMis[i] <- ySums/length(X1)
          aveMisNames[i] <- y
          i <- i + 1
        }
        Ystar <- aveMisNames[which.max(aveMis)]
        X1 <- c(X1, Ystar)
      }
      #end Initialization
      
      #m <- ceiling((n-s)/t + 1)
      m <- ceiling(((length(V)-initSize)/stepSize) + 1)
      Ynames <- colnames(MI[ , !colnames(MI) %in% X1])
      
      #this initializes curriculum to a list containing m vectors
      curriculum <- vector("list",m)
      curriculum[[1]] <- X1
      
      for (i in 2:m)
      {
        curriculum[[i]] <- curriculum[[i-1]]
        
        for(j in 1:stepSize)
        {
          #will only evaluate if Y is non empty
          
          if(length(Ynames))
          {
            aveMis <- c()#numeric(length = length(Ynames))
            aveMisNames <- c()#character(length = length(Ynames))
            currentX <- c()
            
            for (y in Ynames)
            {
              k <- 1
              ySums = 0
              for (x in curriculum[[i]])
              {
                ySums = ySums + MI[y,x]
              }
              aveMis <- c(aveMis, (ySums/length(curriculum[[i]])))
              aveMisNames <- c(aveMisNames, y)
              k <- k + 1
            }
            Ystar <- aveMisNames[which.max(aveMis)]
            Ynames <- Ynames[!(Ynames %in% Ystar)]
            
            currentX = c(currentX, Ystar)
            curriculum[[i]] <- c(curriculum[[i]], currentX)
            j <- j + 1
          }
        }
      }
      return(curriculum)
    }
  }
}

######################################################
# Author: Mohammad Ali Javidian & Tristan Klintworth
# Purpose: Group the data (D1,D2,..,Di) based on all different configs and 
#   find group Di with the largest number of duplicated rows
# Returns subset of data matching Di
# Params:
#     dataset- a .csv file contains the whole data
#     restOFvars- should be a list of column variables (i.e., the rest of the vars in the current iteration)
###############################################################################

getDataset <- function(dataset, restOFvars){
  data<-dataset[,restOFvars]
  ### partitioning and counting the number of rows in a data frame in R based on group [duplicated rows]
  ds <- aggregate(list(numdup=rep(1,nrow(data))), data, length)
  #print(ds)
  #print(ds[which(ds$numdup==max(ds$numdup)),!names(ds) %in% c("numdup")])
  maxD <- ds[which(ds$numdup==max(ds$numdup)),!names(ds) %in% c("numdup")]
  #print(maxD)
  ### if there are more than one max, returns the sub dataset corresponding
  ### to the first D_i with largest number of duplicated rows
  indices <- which(apply(data, 1, function(x) all(x == maxD[1,])))
  return(dataset[indices,])
}



###########################################################################
# Author: Tristan Klintworth
# Purpose: Construct the full Bayesian Network using curriculum learning
# Description : 
# Parameters : 
#    data - data to learn BN
#    initSize - integer to be used to initialize the curriculum
#    stepSize - integer to be used as step size in ConstructCurriculum function call
#    numeric - boolean value to be used to ConstructCurriculum function call to determine the type of output
#    alpha - custom alpha value to be used in hill climbing algorithm
#    size - how many observations in the dataset to use (useful to subset a large dataset)
#    score -  label of the network score to be used in the algorithm
#         * the Bayesian Information Criterion score (bic), which is equivalent to the Minimum Description Length (MDL) and is also known as Schwarz Information Criterion.
#         * the logarithm of the Bayesian Dirichlet equivalent score (bde), a score equivalent Dirichlet posterior density.
#         * the Akaike Information Criterion score (aic).
#         * the multinomial log-likelihood (loglik) score, which is equivalent to the entropy measure used in Weka.
# Each Iteration until you've recovered the full BN : 
#   MI ---> Order(V)
#   stepsize := 4 and start with 4 vars
#   find the correct Di by calling the getDataset function
#   call mmpc, passing in Di and the original dataset (?) in order to recover a skeleton of the BN for the current Di
#   call hc, passing in Di and the dataset (?) in order to orient the direction of the previous learned skeleton based on a score function
############################################################################


recoverWholeBN <- function(data, initSize, stepSize = 4, numeric = FALSE, pval = 0.05, size, score)
{
  #Get subset of the data as specified by the user
  dataSize <- data[1:size,]
  print(dataSize)
  
  # Get the full curriculum list including all steps
  curriculum <- ConstructCurriculum(data, initSize = initSize, stepSize = stepSize, numeric = numeric)
  print("Full curriculum list")
  print(curriculum)
  # Run the MMPC algorithm to generate the parents and children ( PC ) set Si for each Xi .
  S <- mmpc(data)
  Smat <- amat(S)
  #print(S)
  #print("Full adjacency matrix of the dataset")
  #print(Smat[])
  
  #sub matrix should only be the variables in Xi
  #same as subS <- Smat[curriculum[[1]],curriculum[[1]]]
  subS <- lapply(curriculum, function(x) Smat[x,x])
  #print("All submatrices corresponding to the adjacency matrices of Xi's from curriculum list")
  print(subS[[1]])
  
  #The first variable set 
  X1 <- curriculum[[1]]
  
  # Manipulate the first variables to be able to be passed to bnlearn's model2network function
  # This is in order to create an empty network with just the initialization variables
  initNet<-toString(unlist(X1))
  initNet<-gsub(', ',"][",initNet)
  initNet<- paste0("[",initNet,"]")
  print(paste("Nodes in initial net, using variables in X(1) : ", initNet))
  
  # Should work exactly the same as a call with a string in this format : plot(model2network("[VALV][MINV][VLNG][PVS]"))
  G0 <- model2network(initNet)
  plot(G0)
  # OR, i wasted my time doing the above, bnlearn has a function for this: plot(empty.graph(X1)) works
  
  D1 <- getDataset(data, names(data[,!(names(data) %in% c(curriculum[[1]]))]))
  write.csv(D1,file = "D1.csv",row.names = FALSE)
  
  G1 <- hc(D1[,curriculum[[1]]], score = score, alpha = pval) #blacklist = blacklists[[1]])#, blacklist = !(arcs(S1))) G1 <- hc(D1[,X1], blacklist = blacklists[[1]])
  plot(G1)
  print(G1)
  print(amat(G1))
  
  # Build list of blacklists for each Xi
  bl <- c()
  blacklists <- c()
  for(x in subS){
    #print(x)
    for(i in 1:(nrow(x)-1)){
      for(j in (i+1):(nrow(x))){
        if(x[i,j] == 0){
          # Put [i,j] into the blacklist
          bl <- rbind(bl,c(rownames(x)[i],colnames(x)[j]))
        }
        if(x[j,i] == 0){
          # Put [j,i] into the blacklist
          bl <- rbind(bl,c(colnames(x)[j],rownames(x)[i]))
        }
      }
    }
    blacklists <- c(blacklists, list(bl))
  }
  #print("blacklists corresponding to each iteration")
  print(blacklists[[1]])
  
  
  # Perform every iteration
  to <- length(curriculum)
  offset <- 1
  for(i in (1 + offset):to)
  {
    print(i)
    
    Di <- getDataset(data, names(data[,!(names(data) %in% c(curriculum[[i]]))]))
    
    # this gives the "new variables" - the start point of the search is set to the partial network 
    # learned in the previous stage plus the new variables without any edge attached,
    newVars <- curriculum[[i]][!curriculum[[i]] %in% curriculum[[i-1]]]
    print(newVars)
    newVars <- empty.graph(newVars)
    print(modelstring(newVars))
    # Important.. use to add the 3 unconnected nodes
    if(i > 2)
    {
      start1 <- model2network(paste(modelstring(Gi),modelstring(newVars),sep = ""))
      #print(paste(modelstring(Gi),modelstring(newVars),sep = ""))
      #print(start1)
    }
    else #For the 2nd iteration, first G was built outside of the loop
    {
      start1 <- model2network(paste(modelstring(G1),modelstring(newVars),sep = ""))
      print(paste(modelstring(G1),modelstring(newVars),sep = ""))
      print(start1)
      print(amat(start1))
    }
    
    Gi <- hc(Di[,curriculum[[i]]], start = start1, blacklist = blacklists[[i]], score = score, alpha = pval) #D2[,curriculum[[2]]] blacklist = bl) #,start = G1) whitelist = arcs(S2)
    plot(Gi)
    
    print(Gi)
  }
  return(Gi)

  
}
