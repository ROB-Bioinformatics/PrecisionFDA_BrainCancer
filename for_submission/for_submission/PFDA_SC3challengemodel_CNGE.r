clin <- read.table("data_files\\sc3_Phase1_CN_GE_Outcome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
phen <- read.table("data_files\\sc3_Phase1_CN_GE_Phenotype.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
##Transform the phenotypic data into numeric coding for each of the phenotypes
phen[phen==""] <- NA
phen[phen==" "] <- NA
phen[phen=="UNKNOWN"] <- NA
phen[phen=="UNCLASSIFIED"] <- NA
##Numeric coding for the WHO grading: 1-4
phen$WHO_GRADING[phen$WHO_GRADING=="I"] <- 1
phen$WHO_GRADING[phen$WHO_GRADING=="II"] <- 2
phen$WHO_GRADING[phen$WHO_GRADING=="III"] <- 3
phen$WHO_GRADING[phen$WHO_GRADING=="IV"] <- 4
##Numeric coding for the Race, Sex and Cancer type parameters: 1=yes; 0=no
phenrace <- subset(phen, select=c("PATIENTID","RACE"))
phenrace[,2] <- paste("Race",phenrace[,2], sep="_")
phenrace <- reshape2::dcast(phenrace, PATIENTID ~ RACE)
phensex <- subset(phen, select=c("PATIENTID","SEX"))
phensex[,2] <- paste("Sex",phensex[,2], sep="_")
phensex <- reshape2::dcast(phensex, PATIENTID ~ SEX)
phencancer <- subset(phen, select=c("PATIENTID","CANCER_TYPE"))
phencancer[,2] <- paste("CancerType",phencancer[,2], sep="_")
phencancer <- reshape2::dcast(phencancer, PATIENTID ~ CANCER_TYPE)
phenwho <- phen[,c(1,4)]
phencode <- dplyr::left_join(phensex,phenrace)
phencode <- dplyr::left_join(phencode,phencancer)
##Join the clinical and transformed phenotypic data
clinphen <- dplyr::inner_join(clin,phencode)
rownames(clinphen) <- clinphen[,1]
clinphen <- clinphen[,-1]
clinphen[is.na(clinphen)] <- 0
clinphen[clinphen != 0] <- 1
clinphen$Sex_NA[clinphen$Sex_NA == 1] <- -1
clinphen$Race_NA[clinphen$Race_NA == 1] <- -1
clinphen$CancerType_NA[clinphen$CancerType_NA == 1] <- -1
clinphen <- data.frame(PATIENTID=rownames(clinphen),clinphen)
clinphen <- dplyr::left_join(phenwho,clinphen)
rownames(clinphen) <- clinphen[,1]
clinphen <- clinphen[,-1]
clinphen[is.na(clinphen)] <- 0
for(i in 1:ncol(clinphen))
clinphen[,i] <- as.numeric(clinphen[,i])
##Remove effects of NAs
clinphen[,3] <- clinphen[,3] + clinphen[,5]
clinphen[,4] <- clinphen[,4] + clinphen[,5]
clinphen[,6] <- clinphen[,6] + clinphen[,8]
clinphen[,7] <- clinphen[,7] + clinphen[,8]
clinphen[,9] <- clinphen[,9] + clinphen[,8]
clinphen[,10] <- clinphen[,10] + clinphen[,13]
clinphen[,11] <- clinphen[,11] + clinphen[,13]
clinphen[,12] <- clinphen[,12] + clinphen[,13]
clinphen[,14] <- clinphen[,14] + clinphen[,13]
##Get rid of NAs
clinphen <- clinphen[,-c(5,8,13)]
##Use GLM method for selecting relevant phenotypic parameters associated with the survival status
f1 <- as.formula(paste("SURVIVAL_STATUS ~", paste(colnames(clinphen)[c(1,3:ncol(clinphen))], collapse= "+")))
glmcnge <- glm(f1, data=clinphen, family=gaussian())
summary(glmcnge)
siggenes <- rownames(subset(as.data.frame(summary(glmcn)[12]), coefficients.Pr...t.. < 0.05))[-1]
clindataphen <- subset(clinphen, select=colnames(clinphen) %in% c("SURVIVAL_STATUS",siggenes))
clindataphen <- data.frame(PATIENTID=rownames(clindataphen),clindataphen)

##correlation of significant phenotypes with GE
dataexp <- read.table("data_files\\sc3_Phase1_CN_GE_FeatureMatrix.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
clindataexp <- dplyr::inner_join(clindataphen,dataexp)
rownames(clindataexp) <- clindataexp[,1]
clindataexp <- clindataexp[,-1]
corge <- cor(clindataexp[,3:ncol(clindataexp)], clindataexp[,1:2], method="pearson")
corgecut <- subset(as.data.frame(corge), (SURVIVAL_STATUS > 0.4 | SURVIVAL_STATUS < -0.4| CancerType_GBM > 0.4 | CancerType_GBM < -0.4) & ((SURVIVAL_STATUS > 0 & CancerType_GBM > 0) | (SURVIVAL_STATUS < 0 & CancerType_GBM < 0)))
clindataexp <- subset(clindataexp, select=colnames(clindataexp) %in% c("SURVIVAL_STATUS",siggenes,rownames(corgecut)))

##Elastic Net 
library(glmnet)
library(pROC)
n <- nrow(clindataexp)
sample <- sample(seq(n), size = n * 0.6, replace = FALSE)
trainexp <- clindataexp[sample,]
testexp <- clindataexp[-sample,]
surv1exp <- subset(trainexp, SURVIVAL_STATUS==1)
surv0exp <- subset(trainexp, SURVIVAL_STATUS==0)
n1 <- nrow(surv1exp)
n0 <- nrow(surv0exp)
##Randomly choose n0 number of (with survival status=0) from n1 number of samples (with survival status=1) for the training set, 
##so that the number of samples in each group (survival status 0 or 1) will be the same in the training set
sample1 <- sample(seq(n1), size = n0, replace = FALSE)
surv1expdata <- surv1exp[sample1,]
#The rest of samples with survival status=1 can be used in the testing
surv1testdata <- surv1exp[-sample1,]
trainexp <- rbind(surv1expdata,surv0exp)
testexp <- rbind(testexp,surv1testdata)
##Do not run..The following lines describe the model training and parameter optimization.. The optimized pre-built model is saved as CNmodel.RData (that we loaded below for testing)
#a <- seq(0.1, 0.9, 0.05)
##Optimize lambda and alpha hyperparameters by 10-fold cross-validation
#search <- foreach(i = a, .combine = rbind) %dopar% {
#  cv <- glmnet::cv.glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "binomial", nfold = 10, type.measure = "class", parallel = TRUE, alpha = i)
#  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
#}
#cv3 <- search[search$cvm == min(search$cvm), ]
#md3 <- glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
##Load the pre-built elastic net model
data(CNGEmodel)
pathselect <- coef(md3)
pathpar <- as.matrix(pathselect)
selectlinc <- unique(rownames(subset(as.data.frame(pathpar), s0!=0))[-1])
write.table(subset(as.data.frame(pathpar), s0!=0),"sc3_Phase1_CNGE_ElNet_parameters.txt", sep="\t", quote=FALSE)
predicttest <- as.numeric(predict(md3, as.matrix(trainexp[,2:ncol(trainexp)]), type = "class")[,1])
classtest <- as.numeric(trainexp[,1])
elnetroct <- roc(as.numeric(classtest), as.numeric(predicttest))
png("ElNet_SC3modelCNGEROC_Allsampletrain.png", width=2800, height=2600, res=300)
plot(elnetroct)
dev.off()

##Data frame with prediction and true class
predictclass <- data.frame(PATIENTID=rownames(clindataexp),Prediction=predicttest, Trueclass=classtest)
write.table(predictclass, "sc3_Phase1_CNGE_predictions.txt", sep="\t", quote=FALSE, row.names=FALSE)


##Create confusion matrix
tp=0
tn=0
fp=0
fn=0
for(i in 1:length(predicttest)){
if(predicttest[i]==1 & classtest[i]==1)
tp <- tp+1
if(predicttest[i]==1 & classtest[i]==0)
fp <- fp+1
if(predicttest[i]==0 & classtest[i]==0)
tn <- tn+1
if(predicttest[i]==0 & classtest[i]==1)
fn <- fn+1
}

percentp <- tp/(tp+fn)
percentn <- tn/(tn+fp)

confusion <- matrix(c(tp,fp,fn,tn), nrow=2, ncol=2)
rownames(confusion) <- c("True-positive","True-negative")
colnames(confusion) <- c("Prediction-positive","Prediction-negative")
write.table(confusion, "sc3_Phase1_CNGE_confusion-matrix.txt", sep="\t", quote=FALSE)

