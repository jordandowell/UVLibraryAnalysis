## ========================================================================== ##
## LEGACY MANUSCRIPT MODEL-TESTING SCRIPT
## ========================================================================== ##
## This file preserves the supervised PCA/LDA and Naive Bayes experiments used
## during manuscript development. It remains available for reproducibility but
## is not part of the current compound-first Betula workflow.
##
## Important limitations
##   * Some exploratory fits assess the same library observations used to fit.
##   * Poor class confusion matrices motivated the newer open-set approach.
##   * This script does not generate the current practical Betula decisions.
##
## Run RUN_ANALYSIS.R for the current analysis. The historical calculations
## below are intentionally left unchanged so the manuscript work is auditable.
##
## Data contract
##   This script reads the authenticated standards from library/ and historical
##   manuscript-era sample inputs referenced below. It is not genus-generic and
##   does not produce the landscape abundance matrix.

# Packages required by the manuscript-era workflow.
library(MASS)
library(factoextra)
library(tidyverse)
library(dplyr)
library(wavethresh)
library(sgolay)

source("00_library_inputs.R")
library_bundle <- load_library_inputs()
AbsLib <- library_bundle$spectra
# Standardize manuscript-era figure labels while retaining NFP as an acronym.
AbsLib$Compound <- ifelse(
  tolower(trimws(AbsLib$Compound)) == "p-cresol", "p-Cresol",
  paste0(toupper(substr(trimws(AbsLib$Compound), 1, 1)),
         substr(trimws(AbsLib$Compound), 2, nchar(trimws(AbsLib$Compound))))
)
AbsLib$Short_Class <- ifelse(
  toupper(trimws(AbsLib$Short_Class)) == "NFP", "NFP",
  tools::toTitleCase(tolower(trimws(AbsLib$Short_Class)))
)

# Historical min-max normalization removes absolute-intensity differences.
# Select wavelength headers explicitly so custom libraries may contain more or
# fewer metadata columns without changing the manuscript-era calculation.
library_wavelength_columns <- which(is.finite(suppressWarnings(as.numeric(
  sub("nm$", "", names(AbsLib), ignore.case = TRUE)
))))
SmoothedNorm_Lib <- as.data.frame(AbsLib[, library_wavelength_columns]) %>%
  rowwise() %>%
  mutate(across(everything(), ~ (.x - min(c_across(everything()), na.rm = TRUE)) /
                  (max(c_across(everything()), na.rm = TRUE) - min(c_across(everything()), na.rm = TRUE)))) %>%
  ungroup()

norm_data_lib <- as.matrix(SmoothedNorm_Lib)
my_filter = 3
Smoothed_Lib <- sgolayfilt(norm_data_lib, p = my_filter, n = my_filter + 3 - my_filter%%2, m = 0, ts = 1)

extract_wavelet_features <- function(spectra_matrix, wavelet="la8") {

  spectra_matrix <- as.matrix(spectra_matrix)
  n_samples <- nrow(spectra_matrix)

  feature_list <- list()

  for (i in 1:n_samples) {

    spectrum <- as.numeric(spectra_matrix[i, ])

    # wavethresh requires length = power of 2
    n <- length(spectrum)
    next_pow2 <- 2^ceiling(log2(n))

    if (n != next_pow2) {
      spectrum <- c(spectrum, rep(0, next_pow2 - n))
    }

    # Perform discrete wavelet transform
    wt <- wd(spectrum, filter.number=8, family="DaubLeAsymm")

    n_levels <- wt$nlevels

    energy <- numeric(n_levels)
    mean_coef <- numeric(n_levels)
    var_coef <- numeric(n_levels)

    # Extract features per scale
    for (lvl in 0:(n_levels-1)) {

      coef <- accessD(wt, lvl)

      energy[lvl+1] <- sum(coef^2)
      mean_coef[lvl+1] <- mean(coef)
      var_coef[lvl+1] <- var(coef)
    }

    features <- c(energy, mean_coef, var_coef)

    feature_list[[i]] <- features
  }

  feature_matrix <- do.call(rbind, feature_list)

  colnames(feature_matrix) <- c(
    paste0("energy_L", 1:n_levels),
    paste0("mean_L", 1:n_levels),
    paste0("var_L", 1:n_levels)
  )

  return(feature_matrix)
}

# ---------------------------
# RUN FEATURE EXTRACTION
# ---------------------------

wave_features <- extract_wavelet_features(Smoothed_Lib)
wave_features_cleaned <- wave_features[, -c(17)]
final_wavefeatures <- cbind(AbsLib$Short_Class, AbsLib$Compound, wave_features_cleaned)
colnames(final_wavefeatures)[1] <- "Known_Class"
colnames(final_wavefeatures)[2] <- "Compound"

TestingData_Library <- as.data.frame(final_wavefeatures)

NumericSmoothed <-TestingData_Library[, -c(1:2)] %>%
  mutate(across(where(is.character), as.numeric))
NumericSmoothed <- cbind(TestingData_Library[, c(1:2)], NumericSmoothed)


##Generate principal component data ###########################################

PCA_Data <- prcomp(NumericSmoothed[, -c(1:2)], center = TRUE, scale. = TRUE)

pca_scores <- PCA_Data$x


## Test PCA-LDA Models #########################################################

nPCs <- 6

LDA_Data <- as.data.frame(cbind(TestingData_Library$Known_Class, pca_scores[, c(1:nPCs)]))
LDA_Data <-LDA_Data[, -1] %>%
  mutate(across(where(is.character), as.numeric))

LDA_Data <- as.data.frame(cbind(TestingData_Library$Known_Class, LDA_Data))

colnames(LDA_Data)[1] <- "Known_Class"

# 2. Fit the LDA model using a formula (Target ~ Predictors)
lda_result_6PC <- lda(Known_Class ~ ., data = LDA_Data)

# 3. View model output (Shows prior probabilities and coefficients)
print(lda_result_6PC)

# 4. Project data onto the new linear discriminants
lda_pred_6PC <- predict(lda_result_6PC)
lda_scores_6PC <- lda_pred_6PC$x

# Generate the confusion matrix
conf_matrix_6PC <- table(Actual = LDA_Data$Known_Class, Predicted = lda_pred_6PC$class)


per_class_metrics_noloocv_6PCs <- calculate_multiclass_metrics(conf_matrix_6PC)

# Store results for each LOOCV Iteration
PCALDA_Results_noloocv_6PC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_noloocv_6PC <- data.frame(
  Recall <- sum(per_class_metrics_noloocv_6PCs$TP)/(sum(per_class_metrics_noloocv_6PCs$TP) + sum(per_class_metrics_noloocv_6PCs$FN)),
  Precision <- sum(per_class_metrics_noloocv_6PCs$TP)/(sum(per_class_metrics_noloocv_6PCs$TP) + sum(per_class_metrics_noloocv_6PCs$FP)),
  Specif <- sum(per_class_metrics_noloocv_6PCs$TN)/(sum(per_class_metrics_noloocv_6PCs$TN) + sum(per_class_metrics_noloocv_6PCs$FP)),
  Acc <- sum(per_class_metrics_noloocv_6PCs$TP + per_class_metrics_noloocv_6PCs$TN)/(sum(per_class_metrics_noloocv_6PCs$TP) + sum(per_class_metrics_noloocv_6PCs$FN) + sum(per_class_metrics_noloocv_6PCs$TN) + sum(per_class_metrics_noloocv_6PCs$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_noloocv_6PC$Predicted_Class), actual = as.factor(PCALDA_Results_noloocv_6PC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_noloocv_6PC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                 "MCC", "Fmeasure")
print(PCALDA_noloocv_6PC)



## Leave One Out Cross Validation ##############################################

for (i in 1:n) {

  ## Split data
  train_data <- NumericSmoothed[-i, ]
  test_data  <- NumericSmoothed[i, , drop = FALSE]

  removed_class <- NumericSmoothed[i, "Compound"]

  subsetted_data <- train_data[train_data$Compound != removed_class,]

  ## Separate predictors and response
  train_x <- subsetted_data[, -c(1:2)]
  train_y <- as.factor(subsetted_data$Known_Class)

  test_x <- test_data[, -1]

  ## PCA on training data only
  pca_model <- prcomp(train_x,
                      center = TRUE,
                      scale. = TRUE)

  ## Project training and test data
  train_scores <- predict(pca_model, train_x)[, 1:nPCs, drop = FALSE]
  test_scores  <- predict(pca_model, test_x)[, 1:nPCs, drop = FALSE]

  ## Build LDA training data
  lda_train <- data.frame(
    Known_Class = train_y,
    train_scores
  )

  ## Fit LDA
  lda_model <- lda(Known_Class ~ ., data = lda_train)

  ## Predict training data
  train_predictions <- predict(
    lda_model,
    newdata = data.frame(train_scores)
  )$class

  train_conf_matrix <- table(
    predicted = train_predictions,
    actual = train_y
  )

  loocv_accuracy[i] <-
    sum(diag(train_conf_matrix)) / sum(train_conf_matrix)

  ## Predict held-out sample
  predicted_removed[i] <- as.character(
    predict(
      lda_model,
      newdata = data.frame(test_scores)
    )$class
  )
}

predicted_removed


conf_matrix_6PC_LOOCV <- table(Actual = LDA_Data$Known_Class, Predicted = predicted_removed)

per_class_metrics_loocv_6PCs <- calculate_multiclass_metrics(conf_matrix_6PC_LOOCV)
print(per_class_metrics_loocv_6PCs)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_6PC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_6PC <- data.frame(
  Recall <- sum(per_class_metrics_loocv_6PCs$TP)/(sum(per_class_metrics_loocv_6PCs$TP) + sum(per_class_metrics_loocv_6PCs$FN)),
  Precision <- sum(per_class_metrics_loocv_6PCs$TP)/(sum(per_class_metrics_loocv_6PCs$TP) + sum(per_class_metrics_loocv_6PCs$FP)),
  Specif <- sum(per_class_metrics_loocv_6PCs$TN)/(sum(per_class_metrics_loocv_6PCs$TN) + sum(per_class_metrics_loocv_6PCs$FP)),
  Acc <- sum(per_class_metrics_loocv_6PCs$TP + per_class_metrics_loocv_6PCs$TN)/(sum(per_class_metrics_loocv_6PCs$TP) + sum(per_class_metrics_loocv_6PCs$FN) + sum(per_class_metrics_loocv_6PCs$TN) + sum(per_class_metrics_loocv_6PCs$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_6PC$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_6PC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_6PC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                          "MCC", "Fmeasure")
print(PCALDA_LOOCV_6PC)

### testing model using 11 PCs (explaining 99% of the total variance) ##########


nPCs <- 11

LDA_Data <- as.data.frame(cbind(TestingData_Library$Known_Class, pca_scores[, c(1:nPCs)]))
LDA_Data <-LDA_Data[, -1] %>%
  mutate(across(where(is.character), as.numeric))

LDA_Data <- as.data.frame(cbind(TestingData_Library$Known_Class, LDA_Data))

colnames(LDA_Data)[1] <- "Known_Class"

# 2. Fit the LDA model using a formula (Target ~ Predictors)
lda_result_11PC <- lda(Known_Class ~ ., data = LDA_Data)

# 3. View model output (Shows prior probabilities and coefficients)
print(lda_result_11PC)

# 4. Project data onto the new linear discriminants
lda_pred_11PC <- predict(lda_result_11PC)
lda_scores_11PC <- lda_pred_11PC$x


# Generate the confusion matrix
conf_matrix_11PC <- table(Actual = LDA_Data$Known_Class, Predicted = lda_pred_11PC$class)


per_class_metrics_noloocv_11PCs <- calculate_multiclass_metrics(conf_matrix_11PC)

# Store results for each LOOCV Iteration
PCALDA_Results_noloocv_11PC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = lda_pred_11PC$class
  #Model_Accuracy = loocv_accuracy
)

PCALDA_noloocv_11PC <- data.frame(
  Recall <- sum(per_class_metrics_noloocv_11PCs$TP)/(sum(per_class_metrics_noloocv_11PCs$TP) + sum(per_class_metrics_noloocv_11PCs$FN)),
  Precision <- sum(per_class_metrics_noloocv_11PCs$TP)/(sum(per_class_metrics_noloocv_11PCs$TP) + sum(per_class_metrics_noloocv_11PCs$FP)),
  Specif <- sum(per_class_metrics_noloocv_11PCs$TN)/(sum(per_class_metrics_noloocv_11PCs$TN) + sum(per_class_metrics_noloocv_11PCs$FP)),
  Acc <- sum(per_class_metrics_noloocv_11PCs$TP + per_class_metrics_noloocv_11PCs$TN)/(sum(per_class_metrics_noloocv_11PCs$TP) + sum(per_class_metrics_noloocv_11PCs$FN) + sum(per_class_metrics_noloocv_11PCs$TN) + sum(per_class_metrics_noloocv_11PCs$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_noloocv_11PC$Predicted_Class), actual = as.factor(PCALDA_Results_noloocv_11PC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_noloocv_11PC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                  "MCC", "Fmeasure")
print(PCALDA_noloocv_11PC)


## Leave one out cross validation for 11 PCs ###################################

for (i in 1:n) {

  ## Split data
  train_data <- NumericSmoothed[-i, ]
  test_data  <- NumericSmoothed[i, , drop = FALSE]

  removed_class <- NumericSmoothed[i, "Compound"]

  subsetted_data <- train_data[train_data$Compound != removed_class,]

  ## Separate predictors and response
  train_x <- subsetted_data[, -c(1:2)]
  train_y <- as.factor(subsetted_data$Known_Class)

  test_x <- test_data[, -1]

  ## PCA on training data only
  pca_model <- prcomp(train_x,
                      center = TRUE,
                      scale. = TRUE)

  ## Project training and test data
  train_scores <- predict(pca_model, train_x)[, 1:nPCs, drop = FALSE]
  test_scores  <- predict(pca_model, test_x)[, 1:nPCs, drop = FALSE]

  ## Build LDA training data
  lda_train <- data.frame(
    Known_Class = train_y,
    train_scores
  )

  ## Fit LDA
  lda_model <- lda(Known_Class ~ ., data = lda_train)

  ## Predict training data
  train_predictions <- predict(
    lda_model,
    newdata = data.frame(train_scores)
  )$class

  train_conf_matrix <- table(
    predicted = train_predictions,
    actual = train_y
  )

  loocv_accuracy[i] <-
    sum(diag(train_conf_matrix)) / sum(train_conf_matrix)

  ## Predict held-out sample
  predicted_removed[i] <- as.character(
    predict(
      lda_model,
      newdata = data.frame(test_scores)
    )$class
  )
}

confusion_matrix_PCALDA_11PC<- table(
  predicted = predicted_removed,
  actual = actual_removed
)

per_class_metrics_loocv_11PCs <- calculate_multiclass_metrics(confusion_matrix_PCALDA_11PC)
print(per_class_metrics_loocv_11PCs)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_11PC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_11PC <- data.frame(
  Recall <- sum(per_class_metrics_loocv_11PCs$TP)/(sum(per_class_metrics_loocv_11PCs$TP) + sum(per_class_metrics_loocv_11PCs$FN)),
  Precision <- sum(per_class_metrics_loocv_11PCs$TP)/(sum(per_class_metrics_loocv_11PCs$TP) + sum(per_class_metrics_loocv_11PCs$FP)),
  Specif <- sum(per_class_metrics_loocv_11PCs$TN)/(sum(per_class_metrics_loocv_11PCs$TN) + sum(per_class_metrics_loocv_11PCs$FP)),
  Acc <- sum(per_class_metrics_loocv_11PCs$TP + per_class_metrics_loocv_11PCs$TN)/(sum(per_class_metrics_loocv_11PCs$TP) + sum(per_class_metrics_loocv_11PCs$FN) + sum(per_class_metrics_loocv_11PCs$TN) + sum(per_class_metrics_loocv_11PCs$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_11PC$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_11PC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_11PC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                "MCC", "Fmeasure")
print(PCALDA_LOOCV_11PC)


complete_LDADAta <- rbind(PCALDA_noloocv_6PC, PCALDA_LOOCV_6PC, PCALDA_noloocv_11PC, PCALDA_LOOCV_11PC)
rownames(complete_LDADAta)[1:4] <- c("PCA-LDA (90% var)", "LOOCV PCA-LDA (90% var)", "PCA-LDA (99% var)", "LOOCV PCA-LDA (99% var)")




##### NBC-KDE Model Testing ####################################################


#initially test model to confirm it can make strong clusters using the entire library
nbckde_model_use <- naive_bayes(Known_Class~., data = NumericSmoothed,
                             usekernel = TRUE, usepoisson = FALSE, laplace = 0.5
)
NumericSmoothed
library_predictions <- predict(nbckde_model_use, data = NumericSmoothed)

confusion_matrix_library <- table(predicted = library_predictions, actual = NumericSmoothed$Known_Class)
print(confusion_matrix_library)


per_class_metrics_NBCKDE_noloocv <- calculate_multiclass_metrics(confusion_matrix_library)


# Store results for each LOOCV Iteration
PCALDA_Result_NBCKDE <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = library_predictions
  #Model_Accuracy = loocv_accuracy
)

PCALDA_NBCKDE <- data.frame(
  Recall <- sum(per_class_metrics_NBCKDE_noloocv$TP)/(sum(per_class_metrics_NBCKDE_noloocv$TP) + sum(per_class_metrics_NBCKDE_noloocv$FN)),
  Precision <- sum(per_class_metrics_NBCKDE_noloocv$TP)/(sum(per_class_metrics_NBCKDE_noloocv$TP) + sum(per_class_metrics_NBCKDE_noloocv$FP)),
  Specif <- sum(per_class_metrics_NBCKDE_noloocv$TN)/(sum(per_class_metrics_NBCKDE_noloocv$TN) + sum(per_class_metrics_NBCKDE_noloocv$FP)),
  Acc <- sum(per_class_metrics_NBCKDE_noloocv$TP + per_class_metrics_NBCKDE_noloocv$TN)/(sum(per_class_metrics_NBCKDE_noloocv$TP) + sum(per_class_metrics_NBCKDE_noloocv$FN) + sum(per_class_metrics_NBCKDE_noloocv$TN) + sum(per_class_metrics_NBCKDE_noloocv$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Result_NBCKDE$Predicted_Class), actual = as.factor(PCALDA_Result_NBCKDE$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_NBCKDE) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                   "MCC", "Fmeasure")
print(PCALDA_NBCKDE)



n <- nrow(NumericSmoothed)

# Storage vectors
loocv_accuracy <- numeric(n)
predicted_removed <- character(n)
actual_removed <- as.character(NumericSmoothed$Known_Class)

#LOOCV Proceedure
for (i in 1:n) {

  ## Split data
  train_data <- NumericSmoothed[-i, ]
  test_data  <- NumericSmoothed[i, , drop = FALSE]

  removed_class <- NumericSmoothed[i, "Compound"]

  subsetted_data <- train_data[train_data$Compound != removed_class,]

  ## Separate predictors and response
  train_x <- subsetted_data[, -c(1:2)]
  train_y <- as.factor(subsetted_data$Known_Class)

  # Train model
  gnb_model <- naive_bayes(
    train_y~.,
    data = subsetted_data,
    usekernel = TRUE,
    usepoisson = FALSE,
    laplace = 0.1
  )

  # Store the true class of the removed sample
  actual_removed[i] <- as.character(test_data$Known_Class)

  # Predict training data to compute model accuracy
  train_predictions <- predict(
    gnb_model,
    newdata = train_x,
    type = "class"
  )

  train_conf_matrix_NBCKDE <- table(
    predicted = train_predictions,
    actual = train_y
  )

  loocv_accuracy[i] <- sum(diag(train_conf_matrix_NBCKDE)) / sum(train_conf_matrix_NBCKDE)

  # Predict the removed sample
  test_x <- test_data[, -c(1,2), drop = FALSE]

  predicted_removed[i] <- as.character(
    predict(
      gnb_model,
      newdata = test_x,
      type = "class"
    )
  )
}

# Define all possible classes
all_classes <- levels(as.factor(NumericSmoothed$Known_Class))

# Convert to factors with identical levels
predicted_removed <- factor(predicted_removed, levels = all_classes)
actual_removed    <- factor(actual_removed, levels = all_classes)

# n × n confusion matrix
confusion_matrix_NBCKDE <- table(
  Predicted = predicted_removed,
  Actual = actual_removed
)

confusion_matrix_NBCKDE

per_class_metrics_NBCKDE <- calculate_multiclass_metrics(confusion_matrix_NBCKDE)


# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_NBCKDE <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_NBCKDE <- data.frame(
  Recall <- sum(per_class_metrics_NBCKDE$TP)/(sum(per_class_metrics_NBCKDE$TP) + sum(per_class_metrics_NBCKDE$FN)),
  Precision <- sum(per_class_metrics_NBCKDE$TP)/(sum(per_class_metrics_NBCKDE$TP) + sum(per_class_metrics_NBCKDE$FP)),
  Specif <- sum(per_class_metrics_NBCKDE$TN)/(sum(per_class_metrics_NBCKDE$TN) + sum(per_class_metrics_NBCKDE$FP)),
  Acc <- sum(per_class_metrics_NBCKDE$TP + per_class_metrics_NBCKDE$TN)/(sum(per_class_metrics_NBCKDE$TP) + sum(per_class_metrics_NBCKDE$FN) + sum(per_class_metrics_NBCKDE$TN) + sum(per_class_metrics_NBCKDE$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_NBCKDE$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_NBCKDE$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_NBCKDE) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                 "MCC", "Fmeasure")
print(PCALDA_LOOCV_NBCKDE)




### Gaussian Naive Bayes Classifier Test #######################################

#initially test model to confirm it can make strong clusters using the entire library
GNBC_model_use <- gaussian_naive_bayes(as.matrix(NumericSmoothed[, -c(1:2)]), NumericSmoothed$Known_Class, data = NumericSmoothed,
)
NumericSmoothed
library_predictions_GNBC <- predict(GNBC_model_use, data = NumericSmoothed)

confusion_matrix_library_GNBC <- table(predicted = library_predictions_GNBC, actual = NumericSmoothed$Known_Class)
print(confusion_matrix_library)


per_class_metrics_GNBC_noloocv <- calculate_multiclass_metrics(confusion_matrix_library_GNBC)


# Store results for each LOOCV Iteration
PCALDA_Result_GNBC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = library_predictions_GNBC
  #Model_Accuracy = loocv_accuracy
)

PCALDA_GNBC <- data.frame(
  Recall <- sum(per_class_metrics_GNBC_noloocv$TP)/(sum(per_class_metrics_GNBC_noloocv$TP) + sum(per_class_metrics_GNBC_noloocv$FN)),
  Precision <- sum(per_class_metrics_GNBC_noloocv$TP)/(sum(per_class_metrics_GNBC_noloocv$TP) + sum(per_class_metrics_GNBC_noloocv$FP)),
  Specif <- sum(per_class_metrics_GNBC_noloocv$TN)/(sum(per_class_metrics_GNBC_noloocv$TN) + sum(per_class_metrics_GNBC_noloocv$FP)),
  Acc <- sum(per_class_metrics_GNBC_noloocv$TP + per_class_metrics_GNBC_noloocv$TN)/(sum(per_class_metrics_GNBC_noloocv$TP) + sum(per_class_metrics_GNBC_noloocv$FN) + sum(per_class_metrics_GNBC_noloocv$TN) + sum(per_class_metrics_GNBC_noloocv$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Result_GNBC$Predicted_Class), actual = as.factor(PCALDA_Result_GNBC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_GNBC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                   "MCC", "Fmeasure")
print(PCALDA_GNBC)




#LOOCV Proceedure
for (i in 1:n) {

  ## Split data
  train_data <- NumericSmoothed[-i, ]
  test_data  <- NumericSmoothed[i, , drop = FALSE]

  removed_class <- test_data$Compound

  ## Remove all replicates of the test compound
  subsetted_data <- train_data[train_data$Compound != removed_class, ]

  ## Training predictors and response
  train_x <- as.matrix(subsetted_data[, -c(1,2)])
  storage.mode(train_x) <- "double"

  train_y <- factor(subsetted_data$Known_Class)

  ## Train model
  gnb_model <- gaussian_naive_bayes(
    x = train_x,
    y = train_y
  )

  ## Store true class
  actual_removed[i] <- as.character(test_data$Known_Class)

  ## Training predictions
  train_predictions <- predict(
    gnb_model,
    newdata = train_x
  )

  train_conf_matrix_NBCKDE <- table(
    predicted = train_predictions,
    actual = train_y
  )

  loocv_accuracy[i] <- sum(diag(train_conf_matrix_NBCKDE)) /
    sum(train_conf_matrix_NBCKDE)

  ## Test predictors
  test_x <- as.matrix(test_data[, -c(1,2), drop = FALSE])
  storage.mode(test_x) <- "double"

  ## Ensure identical column names
  colnames(test_x) <- colnames(train_x)

  ## Predict left-out sample
  predicted_removed[i] <- as.character(
    predict(
      gnb_model,
      newdata = test_x
    )
  )
}

# Define all possible classes
all_classes <- levels(as.factor(NumericSmoothed$Known_Class))

# Convert to factors with identical levels
predicted_removed <- factor(predicted_removed, levels = all_classes)
actual_removed    <- factor(actual_removed, levels = all_classes)

# n × n confusion matrix
confusion_matrix_GNBC <- table(
  Predicted = predicted_removed,
  Actual = actual_removed
)

confusion_matrix_GNBC

per_class_metrics_GNBC <- calculate_multiclass_metrics(confusion_matrix_GNBC)


# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_GNBC <- data.frame(
  Sample = length(NumericSmoothed$Known_Class),
  Actual_Class = NumericSmoothed$Known_Class,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_GNBC <- data.frame(
  Recall <- sum(per_class_metrics_GNBC$TP)/(sum(per_class_metrics_GNBC$TP) + sum(per_class_metrics_GNBC$FN)),
  Precision <- sum(per_class_metrics_GNBC$TP)/(sum(per_class_metrics_GNBC$TP) + sum(per_class_metrics_GNBC$FP)),
  Specif <- sum(per_class_metrics_GNBC$TN)/(sum(per_class_metrics_GNBC$TN) + sum(per_class_metrics_GNBC$FP)),
  Acc <- sum(per_class_metrics_GNBC$TP + per_class_metrics_GNBC$TN)/(sum(per_class_metrics_GNBC$TP) + sum(per_class_metrics_GNBC$FN) + sum(per_class_metrics_GNBC$TN) + sum(per_class_metrics_GNBC$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_GNBC$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_GNBC$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_GNBC) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                   "MCC", "Fmeasure")
print(PCALDA_LOOCV_GNBC)


complete_NBC_Results <- rbind(PCALDA_GNBC, PCALDA_LOOCV_GNBC, PCALDA_NBCKDE, PCALDA_LOOCV_NBCKDE)
rownames(complete_NBC_Results)[1:4] <- c("GNBC", "LOOCV GNBC", "NBC-KDE", "LOOCV NBC-KDE")

complete_NBC_Results





final_data <- rbind(complete_LDADAta, complete_NBC_Results)

final_data

write.csv(final_data, "Updated_ModelMetrics_Sad.csv")

















n_iter <- 1000

all_classes <- unique(NumericSmoothed$Known_Class)

predicted_removed <- vector("list", n_iter)
actual_removed <- vector("list", n_iter)
train_accuracy <- numeric(n_iter)

for (iter in seq_len(n_iter)) {

  ## Randomly choose one compound from each class
  removed_compounds <- sapply(all_classes, function(cl) {

    compounds <- unique(
      NumericSmoothed$Compound[
        NumericSmoothed$Known_Class == cl
      ]
    )

    sample(compounds, 1)
  })

  ## Training and test sets
  train_data <- NumericSmoothed[
    !(NumericSmoothed$Compound %in% removed_compounds),
  ]

  test_data <- NumericSmoothed[
    NumericSmoothed$Compound %in% removed_compounds,
  ]

  ## Training predictors
  train_x <- as.matrix(train_data[, -c(1,2)])
  storage.mode(train_x) <- "double"

  train_y <- factor(train_data$Known_Class)

  ## Train model
  gnb_model <- gaussian_naive_bayes(
    x = train_x,
    y = train_y
  )

  ## Training accuracy
  train_pred <- predict(gnb_model, newdata = train_x)

  train_cm <- table(train_pred, train_y)

  train_accuracy[iter] <-
    sum(diag(train_cm)) / sum(train_cm)

  ## Test predictors
  test_x <- as.matrix(test_data[, -c(1,2)])
  storage.mode(test_x) <- "double"
  colnames(test_x) <- colnames(train_x)

  ## Test predictions
  predicted_removed[[iter]] <-
    as.character(predict(gnb_model, newdata = test_x))

  actual_removed[[iter]] <-
    as.character(test_data$Known_Class)
}

predicted_removed <- unlist(predicted_removed)
actual_removed <- unlist(actual_removed)

levels_all <- levels(factor(NumericSmoothed$Known_Class))

confusion_matrix <- table(
  Predicted = factor(predicted_removed, levels = levels_all),
  Actual = factor(actual_removed, levels = levels_all)
)

confusion_matrix

per_class_metrics_GNBC_loocv_1perclass <- calculate_multiclass_metrics(confusion_matrix)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_GNBC_1perclass <- data.frame(
  Actual_Class = actual_removed,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_GNBC_1perclass <- data.frame(
  Recall <- sum(per_class_metrics_GNBC_loocv_1perclass$TP)/(sum(per_class_metrics_GNBC_loocv_1perclass$TP) + sum(per_class_metrics_GNBC_loocv_1perclass$FN)),
  Precision <- sum(per_class_metrics_GNBC_loocv_1perclass$TP)/(sum(per_class_metrics_GNBC_loocv_1perclass$TP) + sum(per_class_metrics_GNBC_loocv_1perclass$FP)),
  Specif <- sum(per_class_metrics_GNBC_loocv_1perclass$TN)/(sum(per_class_metrics_GNBC_loocv_1perclass$TN) + sum(per_class_metrics_GNBC_loocv_1perclass$FP)),
  Acc <- sum(per_class_metrics_GNBC_loocv_1perclass$TP + per_class_metrics_GNBC_loocv_1perclass$TN)/(sum(per_class_metrics_GNBC_loocv_1perclass$TP) + sum(per_class_metrics_GNBC_loocv_1perclass$FN) + sum(per_class_metrics_GNBC_loocv_1perclass$TN) + sum(per_class_metrics_GNBC_loocv_1perclass$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_GNBC_1perclass$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_GNBC_1perclass$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_GNBC_1perclass) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                 "MCC", "Fmeasure")
print(PCALDA_LOOCV_GNBC_1perclass)




n_iter <- 1000

all_classes <- unique(NumericSmoothed$Known_Class)

predicted_removed <- vector("list", n_iter)
actual_removed <- vector("list", n_iter)
train_accuracy <- numeric(n_iter)





for (iter in seq_len(n_iter)) {

  ## Randomly choose one compound from each class
  removed_compounds <- sapply(all_classes, function(cl) {

    compounds <- unique(
      NumericSmoothed$Compound[
        NumericSmoothed$Known_Class == cl
      ]
    )

    sample(compounds, 1)
  })

  ## Training and test sets
  train_data <- NumericSmoothed[
    !(NumericSmoothed$Compound %in% removed_compounds),
  ]

  test_data <- NumericSmoothed[
    NumericSmoothed$Compound %in% removed_compounds,
  ]

  ## Training predictors
  train_x <- as.matrix(train_data[, -c(1,2)])
  storage.mode(train_x) <- "double"

  train_y <- factor(train_data$Known_Class)

  ## Train model
  gnb_model <- naive_bayes(
     train_y~.,
     data = train_data,
     use_kernel = TRUE,
     use_poisson = FALSE
  )

  ## Training accuracy
  train_pred <- predict(gnb_model, newdata = train_x)

  train_cm <- table(train_pred, train_y)

  train_accuracy[iter] <-
    sum(diag(train_cm)) / sum(train_cm)

  ## Test predictors
  test_x <- as.matrix(test_data[, -c(1,2)])
  storage.mode(test_x) <- "double"
  colnames(test_x) <- colnames(train_x)

  ## Test predictions
  predicted_removed[[iter]] <-
    as.character(predict(gnb_model, newdata = test_x))

  actual_removed[[iter]] <-
    as.character(test_data$Known_Class)
}

predicted_removed <- unlist(predicted_removed)
actual_removed <- unlist(actual_removed)

levels_all <- levels(factor(NumericSmoothed$Known_Class))

confusion_matrix <- table(
  Predicted = factor(predicted_removed, levels = levels_all),
  Actual = factor(actual_removed, levels = levels_all)
)

confusion_matrix

per_class_metrics_NBCKDE_loocv_1perclass <- calculate_multiclass_metrics(confusion_matrix)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_NBCKDE_1perclass <- data.frame(
  Actual_Class = actual_removed,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_NBCKDE_1perclass <- data.frame(
  Recall <- sum(per_class_metrics_NBCKDE_loocv_1perclass$TP)/(sum(per_class_metrics_NBCKDE_loocv_1perclass$TP) + sum(per_class_metrics_NBCKDE_loocv_1perclass$FN)),
  Precision <- sum(per_class_metrics_NBCKDE_loocv_1perclass$TP)/(sum(per_class_metrics_NBCKDE_loocv_1perclass$TP) + sum(per_class_metrics_NBCKDE_loocv_1perclass$FP)),
  Specif <- sum(per_class_metrics_NBCKDE_loocv_1perclass$TN)/(sum(per_class_metrics_NBCKDE_loocv_1perclass$TN) + sum(per_class_metrics_NBCKDE_loocv_1perclass$FP)),
  Acc <- sum(per_class_metrics_NBCKDE_loocv_1perclass$TP + per_class_metrics_NBCKDE_loocv_1perclass$TN)/(sum(per_class_metrics_NBCKDE_loocv_1perclass$TP) + sum(per_class_metrics_NBCKDE_loocv_1perclass$FN) + sum(per_class_metrics_NBCKDE_loocv_1perclass$TN) + sum(per_class_metrics_NBCKDE_loocv_1perclass$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_NBCKDE_1perclass$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_NBCKDE_1perclass$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_NBCKDE_1perclass) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                           "MCC", "Fmeasure")
print(PCALDA_LOOCV_NBCKDE_1perclass)

print(PCALDA_LOOCV_GNBC_1perclass)





nPCs <- 6

all_classes <- unique(NumericSmoothed$Known_Class)

predicted_removed <- vector("list", n_iter)
actual_removed <- vector("list", n_iter)
train_accuracy <- numeric(n_iter)

for (iter in seq_len(n_iter)) {

  ## Randomly choose one compound from each class
  removed_compounds <- sapply(all_classes, function(cl) {

    compounds <- unique(
      NumericSmoothed$Compound[
        NumericSmoothed$Known_Class == cl
      ]
    )

    sample(compounds, 1)
  })

  ## Training and test sets
  train_data <- NumericSmoothed[
    !(NumericSmoothed$Compound %in% removed_compounds), ]

  test_data <- NumericSmoothed[
    NumericSmoothed$Compound %in% removed_compounds, ]

  ## Predictor matrices
  train_x <- as.matrix(train_data[, -c(1,2)])
  test_x  <- as.matrix(test_data[, -c(1,2)])

  storage.mode(train_x) <- "double"
  storage.mode(test_x)  <- "double"

  train_y <- factor(train_data$Known_Class)

  ## -------------------------
  ## PCA on training data only
  ## -------------------------
  pca_model <- prcomp(
    train_x,
    center = TRUE,
    scale. = TRUE
  )

  n_pc <- nPCs

  ## PCA scores
  train_scores <- pca_model$x[, 1:n_pc, drop = FALSE]

  ## Project test data
  test_scores <- predict(
    pca_model,
    newdata = test_x
  )[, 1:n_pc, drop = FALSE]

  ## -------------------------
  ## LDA
  ## -------------------------
  lda_model <- lda(
    x = train_scores,
    grouping = train_y
  )

  ## Training accuracy
  train_pred <- predict(
    lda_model,
    newdata = train_scores
  )$class

  train_cm <- table(
    train_pred,
    train_y
  )

  train_accuracy[iter] <-
    sum(diag(train_cm)) / sum(train_cm)

  ## Test predictions
  predicted_removed[[iter]] <-
    as.character(
      predict(
        lda_model,
        newdata = test_scores
      )$class
    )

  actual_removed[[iter]] <-
    as.character(test_data$Known_Class)
}

predicted_removed <- unlist(predicted_removed)
actual_removed <- unlist(actual_removed)


confusion_matrix <- table(
  Predicted = factor(predicted_removed, levels = levels_all),
  Actual = factor(actual_removed, levels = levels_all)
)

predicted_removed
confusion_matrix

per_class_metrics_6PCs_loocv_1perclass <- calculate_multiclass_metrics(confusion_matrix)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_6PCs_1perclass <- data.frame(
  Actual_Class = actual_removed,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_6PCs_1perclass <- data.frame(
  Recall <- sum(per_class_metrics_6PCs_loocv_1perclass$TP)/(sum(per_class_metrics_6PCs_loocv_1perclass$TP) + sum(per_class_metrics_6PCs_loocv_1perclass$FN)),
  Precision <- sum(per_class_metrics_6PCs_loocv_1perclass$TP)/(sum(per_class_metrics_6PCs_loocv_1perclass$TP) + sum(per_class_metrics_6PCs_loocv_1perclass$FP)),
  Specif <- sum(per_class_metrics_6PCs_loocv_1perclass$TN)/(sum(per_class_metrics_6PCs_loocv_1perclass$TN) + sum(per_class_metrics_6PCs_loocv_1perclass$FP)),
  Acc <- sum(per_class_metrics_6PCs_loocv_1perclass$TP + per_class_metrics_6PCs_loocv_1perclass$TN)/(sum(per_class_metrics_6PCs_loocv_1perclass$TP) + sum(per_class_metrics_6PCs_loocv_1perclass$FN) + sum(per_class_metrics_6PCs_loocv_1perclass$TN) + sum(per_class_metrics_6PCs_loocv_1perclass$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_6PCs_1perclass$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_6PCs_1perclass$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_6PCs_1perclass) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                             "MCC", "Fmeasure")
print(PCALDA_LOOCV_6PCs_1perclass)


print(PCALDA_LOOCV_GNBC_1perclass)

print(PCALDA_LOOCV_NBCKDE_1perclass)




nPCs <- 11

all_classes <- unique(NumericSmoothed$Known_Class)

predicted_removed <- vector("list", n_iter)
actual_removed <- vector("list", n_iter)
train_accuracy <- numeric(n_iter)

for (iter in seq_len(n_iter)) {

  ## Randomly choose one compound from each class
  removed_compounds <- sapply(all_classes, function(cl) {

    compounds <- unique(
      NumericSmoothed$Compound[
        NumericSmoothed$Known_Class == cl
      ]
    )

    sample(compounds, 1)
  })

  ## Training and test sets
  train_data <- NumericSmoothed[
    !(NumericSmoothed$Compound %in% removed_compounds), ]

  test_data <- NumericSmoothed[
    NumericSmoothed$Compound %in% removed_compounds, ]

  ## Predictor matrices
  train_x <- as.matrix(train_data[, -c(1,2)])
  test_x  <- as.matrix(test_data[, -c(1,2)])

  storage.mode(train_x) <- "double"
  storage.mode(test_x)  <- "double"

  train_y <- factor(train_data$Known_Class)

  ## -------------------------
  ## PCA on training data only
  ## -------------------------
  pca_model <- prcomp(
    train_x,
    center = TRUE,
    scale. = TRUE
  )

  n_pc <- nPCs

  ## PCA scores
  train_scores <- pca_model$x[, 1:n_pc, drop = FALSE]

  ## Project test data
  test_scores <- predict(
    pca_model,
    newdata = test_x
  )[, 1:n_pc, drop = FALSE]

  ## -------------------------
  ## LDA
  ## -------------------------
  lda_model <- lda(
    x = train_scores,
    grouping = train_y
  )

  ## Training accuracy
  train_pred <- predict(
    lda_model,
    newdata = train_scores
  )$class

  train_cm <- table(
    train_pred,
    train_y
  )

  train_accuracy[iter] <-
    sum(diag(train_cm)) / sum(train_cm)

  ## Test predictions
  predicted_removed[[iter]] <-
    as.character(
      predict(
        lda_model,
        newdata = test_scores
      )$class
    )

  actual_removed[[iter]] <-
    as.character(test_data$Known_Class)
}

predicted_removed <- unlist(predicted_removed)
actual_removed <- unlist(actual_removed)


confusion_matrix <- table(
  Predicted = factor(predicted_removed, levels = levels_all),
  Actual = factor(actual_removed, levels = levels_all)
)

predicted_removed
confusion_matrix

per_class_metrics_11PCs_loocv_1perclass <- calculate_multiclass_metrics(confusion_matrix)

# Store results for each LOOCV Iteration
PCALDA_Results_LOOCV_11PCs_1perclass <- data.frame(
  Actual_Class = actual_removed,
  Predicted_Class = predicted_removed
  #Model_Accuracy = loocv_accuracy
)

PCALDA_LOOCV_11PCs_1perclass <- data.frame(
  Recall <- sum(per_class_metrics_11PCs_loocv_1perclass$TP)/(sum(per_class_metrics_11PCs_loocv_1perclass$TP) + sum(per_class_metrics_11PCs_loocv_1perclass$FN)),
  Precision <- sum(per_class_metrics_11PCs_loocv_1perclass$TP)/(sum(per_class_metrics_11PCs_loocv_1perclass$TP) + sum(per_class_metrics_11PCs_loocv_1perclass$FP)),
  Specif <- sum(per_class_metrics_11PCs_loocv_1perclass$TN)/(sum(per_class_metrics_11PCs_loocv_1perclass$TN) + sum(per_class_metrics_11PCs_loocv_1perclass$FP)),
  Acc <- sum(per_class_metrics_11PCs_loocv_1perclass$TP + per_class_metrics_11PCs_loocv_1perclass$TN)/(sum(per_class_metrics_11PCs_loocv_1perclass$TP) + sum(per_class_metrics_11PCs_loocv_1perclass$FN) + sum(per_class_metrics_11PCs_loocv_1perclass$TN) + sum(per_class_metrics_11PCs_loocv_1perclass$FP)),
  MCC <- mcc(preds = as.factor(PCALDA_Results_LOOCV_11PCs_1perclass$Predicted_Class), actual = as.factor(PCALDA_Results_LOOCV_11PCs_1perclass$Actual_Class)),
  Fmeasure <- 2*(Precision*Recall)/(Precision + Recall)
)

colnames(PCALDA_LOOCV_11PCs_1perclass) <- c("Recall", "Precision", "Specificity", "Accuracy",
                                           "MCC", "Fmeasure")

KFold_Final <- rbind(PCALDA_LOOCV_6PCs_1perclass, PCALDA_LOOCV_11PCs_1perclass, PCALDA_LOOCV_GNBC_1perclass, PCALDA_LOOCV_NBCKDE_1perclass)

rownames(KFold_Final)[1:4] <- c("PCA-LDA KFold (90% var)", "PCA-LDA KFold (99% var)", "GNBC KFold", "NBC-KDE KFold")


RealFinal_Dataset_VerySad <- rbind(final_data, KFold_Final)
