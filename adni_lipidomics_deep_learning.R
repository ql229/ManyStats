### ManyStats analysis : DeepLearning
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY
## Diagnosis group prediction using clinical variables in the ADNI cohort. https://ida.loni.usc.edu/login.jsp?search=true
## Here, Tensorflow KERAS API was used to create a MLP model. Alzheimer's Disease Neuroimaging Initiative
## Useful video  https://www.youtube.com/watch?v=hd81EH1g1bE

#### One time code.
devtools::install_github("rstudio/tfestimators")
devtools::install_github("rstudio/tensorflow")
devtools::install_github("rstudio/keras")
devtools::install_github("rstudio/tfdatasets")
install_keras()

### Start here
library(keras)

# Read data
data1 <- read.csv(file.choose(), header = T,stringsAsFactors = F)
selectedVars.ind <- c("TOTAL13","APOE4","SPARE_AD","BMI","EDULEVEL") ## We select these variables.

## numerical coding of the variables.
data1$DIAG[which(data1$DIAG=="CN")] <- 0
data1$DIAG[which(data1$DIAG=="LMCI")] <- 2
data1$DIAG[which(data1$DIAG=="AD")] <- 1

# create the input matrix.
data <- data1[,selectedVars.ind]
naind <- which(sapply(1:nrow(data), function(x)  { length(c(grep("N/A",data[x,]),which(is.na(data1$DIAG[x])==T))) })==0) # get rid of NA and #N/A
data <- data[naind,]
data1$DIAG <- as.numeric(data1$DIAG)  #
yvar <- data1$DIAG[naind]
for (i in 1:ncol(data)) {
  data[,i] <- as.numeric(data[,i])
}
data <- as.matrix(data)
data <- scale(data,T)
dimnames(data) <- NULL # remove default names.
data <- normalize(data)

# Data partitioning
set.seed(123) # results repeatable.
ind <- sample(2, nrow(data), replace = T, prob = c(0.7, 0.3))
training <- data[ind==1,]
test <- data[ind==2,]
trainingtarget <- yvar[ind==1]
testtarget <- yvar[ind==2]

# One Hot Encoding
trainLabels <- to_categorical(trainingtarget) #
testLabels <- to_categorical(testtarget)

# create sequential model
model <- keras_model_sequential()
model %>%
  layer_dense(units=100, activation = "relu", input_shape = c(ncol(data))) %>%
  layer_dropout(rate = 0.6) %>% ## very important to ovoid overfitting.
  layer_dense(units=50, activation = "tanh") %>%
  layer_dropout(rate = 0.6) %>%
  layer_dense(units = 3, activation = "softmax")
summary(model)

adam <- optimizer_adam(lr = 0.0005 , beta_1=0.9, beta_2 = 0.999, epsilon=1e-08, decay = 0.0) ## Adam optimizer is better than others for classification problems.

# Compile
model %>%
  compile(loss = 'categorical_crossentropy',
          optimizer = adam,
          metrics = 'accuracy')

# Fit model
history <- model %>%
  fit(training,
      trainLabels,
      epoch = 100,
      validation_split = 0.40) # more the better.
plot(history)

# Model evaluation
model1 <- model %>%
  evaluate(test, testLabels)

# Prediction & confusion matrix
prob <- model %>%
  predict_proba(test)

pred <- model %>%
  predict_classes(test)

table1 <- table(Predicted = pred, Actual = testtarget) # confusion matrix.


