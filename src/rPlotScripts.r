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
png(file = "/Users/haoxiangyang/Desktop/Git/PERT/Writeup/case11_time.png", width= 10,height=6,units = 'in',res = 300);
plot(time11$V1,time11$V2,pch = 20, col = "#377EB8", xlab = "sample size", ylab = "log(Time (sec))",
     ylim = c(0,3500), cex = 2)
points(time11$V1,time11$V3,pch = 20, col = "#E41A1C", cex = 2)
legend("topleft",c("Decomposition","Extensive"),col = c("#377EB8","#E41A1C"),pch = 20);
dev.off();