budget19 <- read.csv("/Users/haoxiangyang/Desktop/Git/PERT/src/budget_19.csv",header = TRUE)
diffData <- budget19[,1:9];
ubData <- budget19[,10:18];
lbData <- budget19[,19:27];
sizeList <- c(10,20,50,75,100,200,300,400,500);
logSize <- log(sizeList,10);
atSize <- logSize*10;
nameSize <- c("log 10","log 20","log 50","log 75","log 100","log 200","log 300","log 400","log 500");

png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/case19_budget_diff.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,4,2,2));
diffPlot <- boxplot(diffData$X10,diffData$X20,diffData$X50,diffData$X75,diffData$X100,
        diffData$X200,diffData$X300,diffData$X400, diffData$X500, at = atSize,names = nameSize);

meanDiff <- colMeans(diffData)
points(atSize,meanDiff,pch = 20, col = "#377EB8");
diffLinFit <- lm(meanDiff~atSize)
lines(atSize, diffLinFit$coefficients[2]*atSize + diffLinFit$coefficients[1])
dev.off();

ubLowerCI = colMeans(ubData) - 1.96 * apply(ubData, 2, sd) / sqrt(20);
ubUpperCI = colMeans(ubData) + 1.96 * apply(ubData, 2, sd) / sqrt(20);
ubMean = colMeans(ubData);
lbLowerCI = colMeans(lbData) - 1.96 * apply(lbData, 2, sd) / sqrt(20);
lbUpperCI = colMeans(lbData) + 1.96 * apply(lbData, 2, sd) / sqrt(20);
lbMean = colMeans(lbData);

png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/case19_budget_ublb.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,4,2,2));
plotCI(ubMean, ui = ubUpperCI, li = ubLowerCI, err = "y", pch = 16,
       col = "#377EB8", xlab = "", ylab = "Upper Bound", ylim = c(100,160),slty = "solid", scol = "#377EB8",
       lwd = 3,axes = F)
axis(side = 2)  
axis(side = 1, at = 1:9,  ## add custom x-axis
     label = sizeList)
plotCI(lbMean, ui = lbUpperCI, li = lbLowerCI, err = "y", pch = 16,
       col = "#E41A1C", xlab = NA, ylab = NA,slty = "solid", scol = "#E41A1C",
       lwd = 3,add = TRUE, axes = F)
legend("topright",c("Upper Bound","Lower Bound"),col = c("#377EB8","#E41A1C"),pch = 20);
dev.off();

values11 <- read.csv("/Users/haoxiangyang/Desktop/Git/PERT/src/values_11.csv",header = TRUE)
altLowerCI = colMeans(values11) - 1.96 * apply(values11, 2, sd) / sqrt(20);
altUpperCI = colMeans(values11) + 1.96 * apply(values11, 2, sd) / sqrt(20);
altMean = colMeans(values11);

png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/case11_altValues.png", width= 10,height=6,units = 'in',res = 300);
plotCI(altMean, ui = altUpperCI, li = altLowerCI, err = "y", pch = 16,
       col = "#377EB8", xlab = "", ylab = "Upper Bound", ylim = c(95,115),slty = "solid", scol = "#377EB8",
       lwd = 3,axes = F)
axis(side = 2)  
axis(side = 1, at = 1:5,  ## add custom x-axis
     label = colnames(values11))
dev.off();
tt1 <- t.test(values11$DET,values11$FULL,alternative = "greater",paired = TRUE)
tt2 <- t.test(values11$EXP,values11$FULL,alternative = "greater",paired = TRUE)
tt3 <- t.test(values11$dOnly,values11$FULL,alternative = "greater",paired = TRUE)
tt4 <- t.test(values11$HOnly,values11$FULL,alternative = "greater",paired = TRUE)

time11 <- read.csv("/Users/haoxiangyang/Desktop/Git/PERT/src/time_11.csv",header = FALSE)
time11 <- read.csv("/Users/haoxiangyang/Desktop/Git/PERT/src/time_14.csv",header = FALSE)
png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/case14_time.png", width= 10,height=6,units = 'in',res = 300);
plot(time11$V1,time11$V2,pch = 20, col = "#377EB8", xlab = "sample size", ylab = "Time (sec)",
     ylim = c(0,10000), cex = 2)
points(time11$V1,time11$V3,pch = 20, col = "#E41A1C", cex = 2)
legend("topleft",c("Decomposition","Extensive"),col = c("#377EB8","#E41A1C"),pch = 20);
dev.off();
dev.off();


#=========================================================================================
# New plots
# Value plot
valueGraph <- read.csv("/Users/haoxiangyang/Desktop/PERT_tests/results/value_figure_title.csv",header = TRUE)
detList <- valueGraph$DET/valueGraph$FULL;
expList <- valueGraph$EXP/valueGraph$FULL;
donlyList <- valueGraph$dOnly/valueGraph$FULL;
honlyList <- valueGraph$HOnly/valueGraph$FULL;
fullList <- valueGraph$FULL/valueGraph$FULL;
test <- rbind(detList,expList,donlyList,honlyList,fullList)
png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/graphValues_col.png", width= 10,height=6,units = 'in',res = 300);
barplot(test,beside=T,legend.text = c("DET","EXP","dOnly","HOnly","FULL"),
        args.legend = list(x=25,y = 2,bty = "n",cex = 1.2),main = "Upper bound estimation of alternative solutions",
        ylab="Upper bound",xpd=FALSE,ylim = c(0.8,2),xlab = "Test Cases",
        names.arg = c("Case 11","Case 14","Case 19","Case 35"),col=brewer.pal(n = 5, name = "RdBu"),
        cex.lab=1.2,cex.axis = 1.2,cex.names = 1.2)
dev.off();

# budget plot
ubmean <- read.csv("/Users/haoxiangyang/Desktop/PERT_tests/results/ubOut.csv",header = FALSE);
lbmean <- read.csv("/Users/haoxiangyang/Desktop/PERT_tests/results/lbOut.csv",header = FALSE);
ubsd <- read.csv("/Users/haoxiangyang/Desktop/PERT_tests/results/ubsdOut.csv",header = FALSE);
lbsd <- read.csv("/Users/haoxiangyang/Desktop/PERT_tests/results/lbsdOut.csv",header = FALSE);

sizeList <- c(10,20,50,100,200,500);

fileList <- c(11,14,19,35);
for (j in 1:4){
  par(mar = c(5,4,2,2));
  ubmeanS <- c();
  ubuS <- c();
  ublS <- c();
  lbmeanS <- c();
  lbuS <- c();
  lblS <- c();
  for (i in 1:6){
    ubmeanS <- c(ubmeanS,ubmean[j,i]);
    ubuS <- c(ubuS,ubmean[j,i] + ubsd[j,i]);
    ublS <- c(ublS,ubmean[j,i] - ubsd[j,i]);
    lbmeanS <- c(lbmeanS,lbmean[j,i]);
    lbuS <- c(lbuS,lbmean[j,i] + lbsd[j,i]);
    lblS <- c(lblS,lbmean[j,i] - lbsd[j,i]);
  }
  outString <- paste("/Users/haoxiangyang/Desktop/Git/PERT/Writeup/graphBudgets_",fileList[j],".png",sep = "")
  png(file = outString, width= 10,height=6,units = 'in',res = 300);
  plotCI(lbmeanS, ui = lbuS, li = lblS, err = "y", pch = 16,
         col = "#377EB8", xlab = "Sample size", ylab = "Value of Bounds",slty = "solid", scol = "#377EB8",
         lwd = 3,axes = F, cex.lab=1.2, main = paste("Case ",fileList[j]), cex.main = 1.5)
  axis(side = 2, cex.axis = 1.2)
  axis(side = 1, at = 1:6,  ## add custom x-axis
       label = sizeList, cex.axis = 1.2)
  plotCI(ubmeanS, ui = ubuS, li = ublS, err = "y", pch = 16,
         col = "#E41A1C", xlab = NA, ylab = NA,slty = "solid", scol = "#E41A1C",
         lwd = 3,add = TRUE, axes = F)
  legend("topright",c("Upper Bound","Lower Bound"),col = c("#E41A1C","#377EB8"),pch = 20, cex = 1.2);
  dev.off();
}

outString <- "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/graphBudgets.png"
png(file = outString, width= 22,height=14,units = 'in',res = 300);
par(mfrow=c(2,2));
for (j in 1:4){
  par(mar = c(5,5,2.5,2.5));
  ubmeanS <- c();
  ubuS <- c();
  ublS <- c();
  lbmeanS <- c();
  lbuS <- c();
  lblS <- c();
  for (i in 1:6){
    ubmeanS <- c(ubmeanS,ubmean[j,i]);
    ubuS <- c(ubuS,ubmean[j,i] + ubsd[j,i]);
    ublS <- c(ublS,ubmean[j,i] - ubsd[j,i]);
    lbmeanS <- c(lbmeanS,lbmean[j,i]);
    lbuS <- c(lbuS,lbmean[j,i] + lbsd[j,i]);
    lblS <- c(lblS,lbmean[j,i] - lbsd[j,i]);
  }
  plotCI(lbmeanS, ui = lbuS, li = lblS, err = "y", pch = 20,
         col = "#377EB8", xlab = "Sample size", ylab = "Value of Bounds",slty = "solid", scol = "#377EB8",
         lwd = 5,axes = F, cex.lab=2, main = paste("Case ",fileList[j]), cex.main = 3)
  axis(side = 2, cex.axis = 2)
  axis(side = 1, at = 1:6,  ## add custom x-axis
       label = sizeList, cex.axis = 2)
  plotCI(ubmeanS, ui = ubuS, li = ublS, err = "y", pch = 20,
         col = "#E41A1C", xlab = NA, ylab = NA,slty = "solid", scol = "#E41A1C",
         lwd = 5,add = TRUE, axes = F)
  legend("topright",c("Upper Bound","Lower Bound"),col = c("#E41A1C","#377EB8"),pch = 20, cex = 2);
}
dev.off();
