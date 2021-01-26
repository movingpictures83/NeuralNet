dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

#Neural Nets
library(randomForest)
library(rpart)
library(rattle)
library(caret)
#library(mlr)
input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  print("READING INPUT FILES...");
  pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }
  t1 <<- read.table(paste(pfix, toString(parameters["training",2]), sep=""), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  t2 <<- read.table(paste(pfix, toString(parameters["clinical",2]), sep=""), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  print("DONE");
  prefix <<- paste(pfix, toString(parameters["prefix", 2]), sep="");
  joinby <<- toString(parameters["joinby", 2])
  classcol <<- toString(parameters["classcol", 2])
  target <<- toString(parameters["target", 2])
  threshold <<- toString(parameters["threshold", 2])
  id <<- toString(parameters["id", 2])
  myX <<- toString(parameters["x", 2])
  trainmethodC <<- toString(parameters["trainControl", 2])
  trainmethod <<- toString(parameters["train", 2])
  #t1 <- as.data.frame(t(t1), stringsAsFactors=FALSE)
#colnames(t1)[1] = "CEL"
# rownames(t1) <- substring(rownames(t1), 2, length(rownames(t1)))
# write.table(t1, file = "shortRMA.csv", sep = ",", row.names = FALSE, col.names = TRUE)
#x <- as.data.frame(merge(t1, t2, by ="CEL", stringsAsFactors=FALSE))

#studyID = unique(x$STUDYID)
}

run <- function() {
   t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
   colnames(t1)[1] <<- joinby
   train_set_size <<- ncol(t1)
   class_index <<- grep(classcol, colnames(t2))
   x <<- as.data.frame(merge(t1, t2, by =joinby, stringsAsFactors=FALSE))
   studyID <<- unique(as.character(unlist(x[id])))
}

output <- function(outputfile) {
for(virus in c(target)){
  
  v1 <- x[x[,id]==virus,]
times = unique(as.numeric(unlist(v1[myX])))
  #times = unique(v1$TIMEHOURS)
  maxAc = 0
  for(t in times){
    print(t)
    #for(threshold in c(0.0002, 0.0005, 0.0007, 0.0009))
    for(thresh in c(as.numeric(threshold)))
    {
      t.x = v1[v1[,myX]==t,]
      t3 = read.csv(paste(prefix,virus,"_",t,".csv",sep=""),header = TRUE)
      # newV = t.x[,c(1,as.numeric(unlist(t3[1]))+1,22279:dim(t.x)[2])]
      
      resCon = sapply(t3[2], function(x) x > thresh)
      fin = t3[2][resCon]
      print(train_set_size+1)
      print(class_index)
      newV = t.x[,c(as.numeric(unlist(t3[1][[1]][1:length(fin)]))+1,(train_set_size+1):dim(t.x)[2])]
      
      set.seed(1283)
      datX = data.matrix(newV[,1:length(fin)])
      rf.label = as.factor(newV[,(length(fin)+class_index)])
      #colnames(datX) <- colnames(newV)[,length(fin)];
      cv.folds <- createMultiFolds(rf.label, k=10, times = 10)
      print("TRAIN CONTROL...")
      fit  = trainControl(method = trainmethodC, number = 10, repeats = 10, index = cv.folds)
      #print(fin)
      #print(length(fin))
      #print(t3)
      #print(resCon)
      #print(datX)
      #print(newV)
      print("TRAINING...")
      res = train(x = datX, y = rf.label, method = trainmethod, tuneLength = 3, ntree = 1000, trControl = fit)
      df = c(max(res$results$Accuracy), t, virus, length(fin), thresh)
      write(df, file = outputfile, sep = ",", append=TRUE)
      if(max(res$results$Accuracy)>maxAc)
      {
        maxAc = max(res$results$Accuracy)
        bestFit = fit
        bestRes = res
        bestTime = t
        bestThreshold = thresh
        bestNumberOfFeatures = length(fin)
      }
    }
  }
  print(virus)
  print(maxAc)
  print(bestTime)
  print(bestThreshold)
  print(bestNumberOfFeatures)
}
}


