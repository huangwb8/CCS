

setGeneric("scaller", function(object, ...) {
  standardGeneric("scaller")
})


#' @rdname CCS-method.scaller
#' @title CCS method: scaller
#' @description \code{scaller} method for \code{CCS} class
#' @param ... Parameters for \code{\link[xgboost]{xgboost}}.
#' @param best_iteration if \code{best_iteration = NULL}(default), cross validation would be used to decide the best iteration round.
#' @param scaller.type One of `xgboost` and `dl`
#' @param net Character (such as `'v20250706'`,`'v20250707'`) or any `torch` net you define.
#' @inheritParams CCSPublicParams
#' @inheritParams GSClassifier::modelData
#' @import xgboost
#' @importFrom GSClassifier modelData
#' @importFrom pROC multiclass.roc roc
#' @importFrom irr kappa2
#' @import tidyr
#' @details
#' The CCS object must be modified by \code{\link{CCS-method.cluster}} or \code{\link{CCS-method.optimizeCluster}}.
#' @seealso \code{\link{ccs}};\code{\link[xgboost]{xgboost}};
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' resCCS <- scaller(resCCS)
#' @exportMethod scaller
setMethod(
  "scaller",
  signature(object='CCS'),
  function(
    object,
    model.dir = NULL,
    Prop = 0.8,
    best_iteration = NULL,
    seed = 456,
    cover = FALSE,
    rm.zero = TRUE,
    scaller.type = c('xgboost','dl')[1],
    net = 'v20250706',
    numCores = 6,
    verbose = TRUE,
    ...
  ){

    # Test
    if(F){
      library(luckyBase)
      Plus.library(c('xgboost', 'GSClassifier','pROC','irr','tidyr'))
      model.dir <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240106'
      path_resCCS <- paste0(model.dir,'/resCCS.rds')
      object <- readRDS(path_resCCS)
      # Renew CCS manually
      classifier <- list(
        "13" = c("13","15","20"),
        "6" = c("48","36","38","3"),
        "45" = c("45","23"),
        "75" = c("75","73"),
        "26" = c("26","25","58","65","27","24"),
        "34" = c("34", "31"),
        "40" = c("56", "43","40"),
        "55" = c("55", "57"),
        "16" = c("16", "66"),
        "82" = c("80", "82"),
        "44" = c("44", "37"),
        "76" = c("76","77","42"),
        "4" = c("4","10"),
        "8" = c("8","17"),
        "18" = c("18","9"),
        "7" = c("7","1","60","5","61"),
        "35" = c("35", "67", "47","41","83","64","22","14","11","12", "69"),
        "70" = c("70", "59")
      )
      object  <- CCS::manualCluster(object , classifier)
      numCores = 6; seed = 456; Prop = 0.8
      verbose = TRUE; cover = TRUE; rm.zero = TRUE
    }

    # CCS normalization: 0-based for xgboost/DL models
    set.seed(seed); seeds <- sample(1:10000, 10, replace = T)
    y2 <- object@Data$CCS
    if(rm.zero){
      target_sample <- !y2 %in% 0
      y2 <- y2[target_sample]
      if(verbose) LuckyVerbose('scaller: Remove ', sum(!target_sample), ' zero value in CCS subtypes...')
    }
    y2_uniqueSubtype <- unique(y2)
    cluster_translator <- data.frame(
      raw = y2_uniqueSubtype,
      adjust = 0:(length(y2_uniqueSubtype) - 1),
      stringsAsFactors = F
    )
    y2_name <- names(y2)
    y2 <- convert(y2, 'raw', 'adjust', cluster_translator) %>% as.integer()
    names(y2) <- y2_name

    # Data
    dat <- data.frame(ID=y2_name, subtype = y2, stringsAsFactors = F)
    dat_modelData <- GSClassifier::modelData(
      design = dat,
      id.col = "ID",
      variable = c("subtype"),
      Prop = Prop,
      seed = seeds[1])
    id_train <- as.character(dat_modelData$Data$Train[,'ID'])
    id_valid <- as.character(dat_modelData$Data$Valid[,'ID'])
    res2 <- object@Data$Probability$d1[y2_name,]
    res2_train <- res2[id_train,] # dim(res2_train)
    y2_train <- y2[id_train]
    res2_valid <- res2[id_valid,] # dim(res2_valid)
    y2_valid <- y2[id_valid]

    # Self-defined parameters
    if(verbose) LuckyVerbose('scaller: adjust parameters...')
    params <- object@Repeat$params
    params_name <- names(params)
    args <- list(...)
    args_name <- names(args)
    for(i in 1:length(args_name)){ # i=1
      args_name_i <- args_name[i]
      if(args_name_i %in% params_name){
        if(args[[args_name_i]] != params[[args_name_i]]){
          if(verbose) LuckyVerbose('Change `', args_name_i, '` from ', params[[args_name_i]], ' to ', args[[args_name_i]], '...', levels = 3)
        }
      } else {
        if(verbose) LuckyVerbose('Add new parameter `', args_name_i, '` = ', args[[args_name_i]], levels = 3)
      }
      params[[args_name_i]] <- args[[args_name_i]]
    }
    params$batch_size <- ifelse(is.null(params$batch_size), length(unique(c(y2_train,y2_valid))) * 3, params$batch_size)

    # Parameters
    params_xg <- params[-match(c('n','sampSize','ptail','nround.mode'), names(params)) %>% .[!is.na(.)]]
    params_xg3 <- params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg)) %>% .[!is.na(.)]]
    params_xg3[['nthread']] <- numCores


    # Model path
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }
    path_scaller <- paste0(model.dir, '/scaller.rds')


    # Check old scaller
    if(!is.null(object@Data[['scaller']]) & !cover){
      stop('scaller: With scaller object in the CCS object. Ignore! You can set `cover = TRUE` to make `scaller` run and cover the old result.')
    }
    if(verbose) LuckyVerbose('scaller: Build subtype caller ...')


    # Subtype Caller
    if((!file.exists(path_scaller)) | cover){

      if(scaller.type == 'xgboost'){
        scaller <- scaller_xgboost(
          x_train = res2_train,
          y_train = y2_train,
          best_iteration = best_iteration,
          params = params,
          seeds = seeds,
          numCores = numCores,
          verbose = verbose
        )
      } else if(scaller.type == 'dl'){
        scaller <- scaller_dl(
          x_train = res2_train,
          y_train = y2_train,
          x_valid = res2_valid,
          y_valid = y2_valid,
          net = net,
          first.width = ifelse(is.null(params$first.width),512, params$first.width),
          num_epochs = params$nrounds,
          batch_size = params$batch_size,
          learning_rate = params$eta,
          seed = seeds[1]
        )

        # Save luz model
        path_scaller <- generate_scaller_path(model.dir)
        luz_save(scaller, path_scaller)
      } else {
        stop('Not supported `scaller.type!`. You should choose one of `dl` and `xgboost`.')
      }
    } else {
      if(verbose) LuckyVerbose('scaller: The subtype caller exists. Use it!')
      scaller <- readRDS(path_scaller)
      params_xg3 <- object@Data[['scaller.parameters']][['Params']]
    }
    object@Data[['scaller']] <- scaller


    # Validation
    if(verbose) LuckyVerbose('scaller: validate `scaller` in the internal cohort...')
    y2_valid_pred <- predict(scaller, res2_valid)
    params_xg4 <- params_xg3
    if(scaller.type == 'xgboost'){
      params_xg4[['best_iteration']] <- best_iteration
    } else if(scaller.type == 'dl'){
      params_xg4[['net']] <- net
      params_xg4[['num_epochs']] <- params$nrounds
      params_xg4[['batch_size']] <- params$batch_size
      params_xg4[['learning_rate']] <- params$eta
      params_xg4[['seed']] <- seeds[1]
      params_xg4[c('nthread','eta','subsample','gamma','alpha','max_depth','colsample_bytree','min_child_weight')] <- NULL
      # params_xg4 <- list(
      #   net = net,
      #   num_epochs = params$nrounds,
      #   batch_size = params$batch_size,
      #   learning_rate = params$eta,
      #   seed = seeds[1]
      # )
    } else {
      stop('Not supported `scaller.type!`. You should choose one of `dl` and `xgboost`.')
    }


    # Performance
    res.all <- CCS:::compareRealPred(y2_valid, y2_valid_pred, cluster_translator)
    if(verbose) LuckyVerbose('scaller: Performance -- multi_auc=', round(res.all$multi_auc[1], 4), '; accuracy=', round(res.all$accuracy[1], 4), '; kappa=', round(res.all$kappa[1], 4))

    # Output
    object@Data[['scaller.type']] <- scaller.type
    object@Data[['scaller.performance']] <- res.all
    object@Data[['scaller.parameters']] <- list(
      Prop = Prop,
      Seed = seed,
      cluster_translator = cluster_translator,
      Params = params_xg4
    )
    if(verbose) LuckyVerbose('scaller: All done!')
    return(object)
  }
)


#' @description Scaller via XGboost
#' @importFrom xgboost xgb.DMatrix xgb.cv xgboost
#' @importFrom tidyr `%>%`
scaller_xgboost <- function(
    x_train,
    y_train,
    best_iteration,
    params,
    seeds,
    numCores,
    verbose
){
  nSubtype <- length(unique(y_train))
  dtrain <- xgb.DMatrix(x_train, label = y_train)

  # Parameters
  params_xg <- params[-match(c('n','sampSize','ptail','nround.mode'), names(params)) %>% .[!is.na(.)]]
  params_xg3 <- params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg)) %>% .[!is.na(.)]]
  params_xg3[['nthread']] <- numCores

  # find best interation
  if(is.null(best_iteration)){
    # xgboost cv for best interation exploration
    if(verbose) LuckyVerbose('scaller: xgboost cv for best interation exploration...')
    set.seed(seeds[2])
    cvRes <- xgb.cv(
      params = params_xg3,
      data = dtrain,
      nrounds=params_xg$nrounds,
      nfold=params_xg$nfold,
      num_class = nSubtype,
      early_stopping_rounds=10,
      objective = "multi:softmax",
      verbose = ifelse(verbose, 1, 0)
    )
    best_iteration <- cvRes$best_iteration
    # mymusic()
  } else {
    if(verbose) LuckyVerbose('scaller: use self-defined best interaction round...')
  }

  # xgboost via best interation
  if(verbose) LuckyVerbose('scaller: xgboost model based on the best interation...')
  set.seed(seeds[2])
  scaller <- xgboost(
    params = params_xg3,
    data = dtrain,
    nrounds = best_iteration,
    num_class = nSubtype,
    objective = "multi:softmax",
    verbose = ifelse(verbose, 1, 0)
  )
  return(scaller)
}


#' @name CCSDeepClassifier-class
#' @title S4 Class for a Custom Deep Learning Classifier
#' @description This class encapsulates a trained deep learning model and the feature names required for prediction. It is designed to provide a standardized interface for deep learning models, particularly those trained with the 'luz' package.
#' @slot model A fitted model object. While typically a 'luz_module_fitted' object from the 'luz' package, this slot is defined as 'ANY' to accommodate various model structures.
#' @slot feature_names A character vector specifying the names and order of the features used to train the model. This ensures that new data for prediction has the correct structure.
#' @author Weibin Huang <hwb2012@@qq.com>
#' @exportClass CCSDeepClassifier
#' @keywords classes
setClass(
  "CCSDeepClassifier",
  slots = c(
    model = "ANY",
    feature_names = "character"
  ),
  prototype = list(
    model = list(),
    feature_names = character(0)
  )
)


#' @description Scaller via Deep network
#' @param x_train A numeric matrix of training features (samples x features).
#' @param y_train A numeric or integer vector of training labels.
#' @param x_valid A numeric matrix of validation features.
#' @param y_valid A numeric or integer vector of validation labels.
#' @param num_epochs The maximum number of epochs to train.
#' @param batch_size The number of samples per batch. Given your hardware, a larger size (e.g., 128, 256) is recommended.
#' @param learning_rate The initial learning rate for the optimizer.
#' @import torch
#' @import luz
#' @return A `luz_module_fitted` object, which contains the trained model, training history, and the best model's weights.
scaller_dl <- function(
    x_train, y_train, x_valid, y_valid,
    net = 'v20250706',
    first.width = 512,
    num_epochs = 100,
    batch_size = 128,
    learning_rate = 1e-3,
    seed = 2025
) {

  # Test
  if(F){
    library(luckyBase)
    np <- c('luz','torch','yardstick','dplyr','tidyr','tibble'); Plus.library(np)
    x_train = res2_train
    y_train = y2_train
    x_valid = res2_valid
    y_valid = y2_valid
    num_epochs = 1000
    batch_size = length(unique(c(y2_train,y2_valid))) * 3
    learning_rate = 1e-3
    seed = 2025
  }

  # 0. Randomization for reproducibility
  if(T){
    set.seed(seed)
    torch_manual_seed(seed)

    # 1) 全局 deterministic 算法（高版本才有）
    if (exists("torch_use_deterministic_algorithms", mode = "function")) {
      torch_use_deterministic_algorithms(TRUE, warn_only = FALSE)
    }

    # 2) cuDNN 相关开关 —— 兼容不同命名
    if (exists("torch_backends_cudnn_set_deterministic", mode = "function")) {
      torch_backends_cudnn_set_deterministic(TRUE)
      torch_backends_cudnn_set_benchmark(FALSE)
    } else if (exists("torch_backends_cudnn_deterministic", mode = "function")) {
      # 早期版本（<=0.6）是直接暴露变量
      torch_backends_cudnn_deterministic(TRUE)
      torch_backends_cudnn_benchmark(FALSE)
    } else {
      message("[Warn] 当前 torch 版本无法显式设置 cuDNN deterministic。")
    }

    # 3) 避免部分 cuBLAS 非确定性 kernel（CUDA ≥10.2 有效）
    if (cuda_is_available())
      Sys.setenv(CUBLAS_WORKSPACE_CONFIG = ":16:8")
  }

  # 1. Determine device (CPU or GPU)
  device <- if (cuda_is_available()) "cuda" else "cpu"
  cat("Using device:", device, "\n")

  # 2. Data Preprocessing: Convert R objects to torch tensors
  all_labels <- c(y_train, y_valid)
  unique_labels <- sort(unique(all_labels))
  num_classes <- length(unique_labels)
  feature_names <- intersect(colnames(x_train), colnames(x_valid))
  x_train <- x_train[, feature_names, drop=FALSE]
  x_valid <- x_valid[, feature_names, drop=FALSE]
  y_train_adj <- y_train + 1; y_valid_adj <- y_valid + 1 # 1-based vector
  x_train_tensor <- torch_tensor(x_train, dtype = torch_float(), device = device)
  y_train_tensor <- torch_tensor(y_train_adj, dtype = torch_long(), device = device)
  x_valid_tensor <- torch_tensor(x_valid, dtype = torch_float(), device = device)
  y_valid_tensor <- torch_tensor(y_valid_adj, dtype = torch_long(), device = device)

  # 3. Create torch datasets and dataloaders
  train_ds <- tensor_dataset(x_train_tensor, y_train_tensor)
  valid_ds <- tensor_dataset(x_valid_tensor, y_valid_tensor)
  train_dl <- dataloader(train_ds, batch_size = batch_size, shuffle = TRUE)
  valid_dl <- dataloader(valid_ds, batch_size = batch_size, shuffle = FALSE, drop_last = TRUE)

  # 4. Define the Neural Network Architecture
  input_features <- ncol(x_train)
  if(is.character(net)){
    net <- get_torch_net(input_features, num_classes, first.width = first.width, net_version = net)
  }

  # 5. Train the Neural Network
  scaller <- net %>%
    setup(
      loss = nn_cross_entropy_loss(),
      optimizer = optim_adam,
      metrics = list(
        accuracy = luz_metric_accuracy(),
        auc = luz_metric_multiclass_auroc(average = 'macro')
      )
    ) %>%
    set_opt_hparams(lr = learning_rate) %>%
    fit(
      train_dl,
      epochs = num_epochs,
      valid_data = valid_dl,
      verbose = TRUE,
      callbacks = list(
        luz_callback_metrics(),
        luz_callback_early_stopping(
          # you have to use "valid_xxx" as a mointor object
          # monitor = "valid_auc", mode = "max",
          monitor = "valid_loss", mode = "min",
          patience = 10,
        )
      )
    )

  # Output
  l <- new(
    'CCSDeepClassifier',
    model = scaller,
    feature_names = feature_names
  )
}


#' @description Preparation of torch net
#' @import torch
get_torch_net <- function(
    input_features,
    num_classes,
    first.width = 512,
    net_version = 'v20250706'
){

  if(net_version == 'v20250706'){
    net <- nn_module(
      "DeepClassifier",
      initialize = function() {
        self$net <- nn_sequential(
          # Layer 1
          nn_linear(input_features, first.width),
          nn_batch_norm1d(first.width),
          nn_relu(),
          nn_dropout(p = 0.3),

          # Layer 2
          nn_linear(first.width, 256),
          nn_batch_norm1d(256),
          nn_relu(),
          nn_dropout(p = 0.25),

          # Layer 3
          nn_linear(256, 128),
          nn_batch_norm1d(128),
          nn_relu(),
          nn_dropout(p = 0.2),

          # Output layer
          nn_linear(128, num_classes)
        )

        # Xavier initialization (optional but recommended)
        self$net$apply(function(m) {
          if (inherits(m, "nn_linear")) {
            nn_init_xavier_uniform_(m$weight)
            if (!is.null(m$bias)) nn_init_constant_(m$bias, 0)
          }
        })
      },
      forward = function(x) {
        self$net(x)
      }
    )

  } else if(net_version == 'v20250707') {

    # ResNet-style MLP

    # Residual Block
    ResidualBlock <- nn_module(
      "ResidualBlock",
      initialize = function(features) {
        self$layers <- nn_sequential(
          nn_linear(features, features),
          nn_batch_norm1d(features),
          nn_relu(),
          nn_dropout(0.3),
          nn_linear(features, features),
          nn_batch_norm1d(features)
        )
        self$relu <- nn_relu()
      },
      forward = function(x) {
        # The core idea of a residual connection:
        # The output of the layers is added to the original input.
        # This "skip connection" allows the model to easily learn an identity function.
        output <- self$layers(x) + x
        self$relu(output)
      }
    )

    # Main classifier using these blocks
    net <- nn_module(
      "DeepClassifier",
      initialize = function() {

        # An initial layer to project input features to the hidden dimension
        self$entry_layer <- nn_sequential(
          nn_linear(input_features, first.width),
          nn_batch_norm1d(first.width),
          nn_relu()
        )

        # Stack multiple residual blocks
        self$residual_blocks <- nn_sequential(
          ResidualBlock(first.width),
          ResidualBlock(first.width)
          # ResidualBlock(first.width) # You can easily add more blocks here to make the network deeper
        )

        # Final layer for classification
        self$output_layer <- nn_linear(first.width, num_classes)
      },

      forward = function(x) {
        x <- self$entry_layer(x)
        x <- self$residual_blocks(x)
        self$output_layer(x) # Output logits
      }
    )

  } else {
    stop(net_version, ': This version is still not supported!')
  }

  return(net)
}


#' @description Prediction function for `scaller_dl`
#' @param scaller A `luz_module_fitted` object, as returned by `luz::fit()`.
#' @param new_data A numeric matrix or data frame with the same features as the training data.
#' @param batch_size Integer. The number of samples to process in each batch. Adjusting this can help manage memory usage, especially on GPUs.
#' @return A numeric vector of predicted class labels. The labels are returned in their original format (e.g., 0-indexed), assuming they were converted to 1-indexed for training.
#' @import torch
#' @import luz
#' @export
predict.CCSDeepClassifier <- function(
    scaller,
    new_data,
    batch_size = 256
) {

  # Test
  if(F){
    new_data <- as.matrix(X_CCSprobability[sampleIndex, feature_names])
    batch_size = 256
  }

  object <- scaller@model
  feature_names <- scaller@feature_names
  new_data <- new_data[,feature_names,drop=FALSE]

  # 1. Set model to evaluation mode
  object$model$eval()

  # 2. Infer device from the model
  device <- if (cuda_is_available()) "cuda" else "cpu"

  # 3. Preprocess the new data
  # Convert the R data frame/matrix to a torch tensor.
  # It's crucial to ensure the data is in matrix format.
  x_tensor <- torch_tensor(
    as.matrix(new_data),
    dtype = torch_float(),
    device = device
  )

  # 4. Create a dataloader for prediction
  # For prediction, we only have features (x), not labels (y).
  # `shuffle = FALSE` is essential to keep the prediction order consistent with the input `new_data`.
  pred_ds <- tensor_dataset(x_tensor)
  pred_dl <- dataloader(pred_ds, batch_size = batch_size, shuffle = FALSE)

  # 5. Get raw model outputs (logits) using the base luz predict function. This performs the prediction in batches.
  raw_preds_tensor <- predict(object, pred_dl)

  # 6. Post-process the predictions
  # `torch_argmax` finds the index of the maximum value along dimension 2 (the class dimension). This index corresponds to the predicted class. The result is 1-indexed (e.g., class 1, 2, 3, ...).
  predicted_indices <- torch_argmax(raw_preds_tensor, dim = 2)

  # 7. Convert to an R vector
  # Move the tensor from the GPU (if applicable) to the CPU and then convert to a standard R vector.
  predicted_classes_r <- as_array(predicted_indices$cpu())

  # 8. Adjust for 0-based indexing
  predicted_labels <- predicted_classes_r - 1

  return(predicted_labels)
}



#' @title Beta-version: Calculate Permutation Feature Importance for a luz Model
#' @description This function computes feature importance for a model trained with the 'luz' package using the permutation-based method. It measures the drop in model performance when a single feature's values are randomly shuffled.
#' @param fitted_luz_model A `luz_module_fitted` object, the result of `luz::fit()`.
#' @param feature_names A character vector of feature names. The order must match the columns of the input tensor `X`.
#' @param metric_name A string specifying the metric to monitor. This must be one of the names reported by `luz` in the validation metrics (e.g., "valid_loss", "valid_multiclass_auroc", "valid_accuracy").
#' @param mode A string, either "max" or "min" or "auto". This defines whether a higher value of the metric is better ("max") or a lower value is better ("min"). If "auto", the function infers the mode by checking if "loss" is in `metric_name`.
#' @param plot A logical value. If `TRUE` (the default), a ggplot of the feature importances is printed.
#' @param ... Additional arguments (currently not used, for future extension).
#' @return
#' A `data.frame` containing the 'Feature' and its calculated 'Importance',
#' sorted from most to least important. If `plot = TRUE`, the data frame is
#' returned invisibly.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
importance_dl <- function(fitted_luz_model,
                          feature_names,
                          metric_name = "valid_loss",
                          mode = "auto",
                          plot = TRUE,
                          ...) {

  # -----------------------------------------------------------------
  # 1. Input Validation and Setup
  # -----------------------------------------------------------------
  # Check for required packages
  if (!requireNamespace("luz", quietly = TRUE) ||
      !requireNamespace("torch", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("forcats", quietly = TRUE)) {
    stop("This function requires 'luz', 'torch', 'dplyr', 'ggplot2', and 'forcats'. Please install them.")
  }

  # Validate input object type
  if (!inherits(fitted_luz_model, "luz_module_fitted")) {
    stop("'fitted_luz_model' must be an object of class 'luz_module_fitted'.")
  }

  # Automatically determine the mode if set to "auto"
  if (mode == "auto") {
    mode <- if (grepl("loss", metric_name, ignore.case = TRUE)) "min" else "max"
    cat(sprintf("Auto-detected mode as '%s' for metric '%s'.\n", mode, metric_name))
  } else if (!mode %in% c("min", "max")) {
    stop("'mode' must be one of 'min', 'max', or 'auto'.")
  }

  # Extract the validation dataloader
  valid_dl <- fitted_luz_model$records$valid_dl[[1]]
  if (is.null(valid_dl)) {
    stop("Could not find a validation dataloader in the fitted model. Importance requires validation data.")
  }

  # -----------------------------------------------------------------
  # 2. Extract Data and Calculate Baseline Performance
  # -----------------------------------------------------------------
  # This is a robust way to get all data from a dataloader
  all_data <- luz::dataloader_get_all_data(valid_dl)
  X_valid <- all_data[[1]]
  y_valid <- all_data[[2]]

  # Check if feature_names length matches data
  num_features <- ncol(X_valid)
  if (length(feature_names) != num_features) {
    stop(sprintf("Length of 'feature_names' (%d) does not match the number of features in the data (%d).",
                 length(feature_names), num_features))
  }

  cat("Calculating baseline performance on validation set...\n")
  # Use luz::evaluate for a consistent and simple way to get metrics
  baseline_metrics <- luz::evaluate(fitted_luz_model, valid_dl)

  if (!metric_name %in% names(baseline_metrics)) {
    stop(sprintf("Metric '%s' not found. Available metrics are: %s",
                 metric_name, paste(names(baseline_metrics), collapse = ", ")))
  }
  baseline_score <- as.numeric(baseline_metrics[[metric_name]])
  cat(sprintf("Baseline %s: %.4f\n", metric_name, baseline_score))

  # -----------------------------------------------------------------
  # 3. Permutation Loop
  # -----------------------------------------------------------------
  importances <- numeric(num_features)
  names(importances) <- feature_names

  # Set model to evaluation mode
  fitted_luz_model$model$eval()

  # Use a progress bar for better user experience
  pb <- utils::txtProgressBar(min = 0, max = num_features, style = 3)

  for (i in 1:num_features) {
    # Create a copy of the validation data to modify
    X_permuted <- X_valid$clone()

    # Permute (shuffle) the i-th column
    # A random permutation of indices ensures we shuffle correctly
    perm_indices <- torch::torch_randperm(nrow(X_permuted), device = X_permuted$device)
    X_permuted[, i] <- X_permuted[perm_indices, i]

    # Create a temporary dataloader for the permuted data
    permuted_ds <- torch::tensor_dataset(X_permuted, y_valid)
    permuted_dl <- torch::dataloader(permuted_ds, batch_size = valid_dl$batch_size)

    # Get performance on the permuted data
    permuted_metrics <- luz::evaluate(fitted_luz_model, permuted_dl)
    permuted_score <- as.numeric(permuted_metrics[[metric_name]])

    # Calculate importance based on the mode
    if (mode == "max") {
      # For metrics like AUC/Accuracy, importance is the drop in score
      importances[i] <- baseline_score - permuted_score
    } else {
      # For metrics like Loss, importance is the increase in score
      importances[i] <- permuted_score - baseline_score
    }

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # -----------------------------------------------------------------
  # 4. Format and Visualize Results
  # -----------------------------------------------------------------
  importance_df <- data.frame(
    Feature = names(importances),
    Importance = as.numeric(importances)
  ) %>%
    # Use fct_reorder to sort features by importance for plotting
    dplyr::mutate(Feature = forcats::fct_reorder(Feature, Importance))

  if (plot) {
    # Create a lollipop chart for an elegant look
    p <- ggplot2::ggplot(importance_df, ggplot2::aes(x = Importance, y = Feature)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, y = Feature, xend = Importance, yend = Feature),
        color = "grey50", linewidth = 0.8
      ) +
      ggplot2::geom_point(
        ggplot2::aes(color = Importance > 0), size = 4, show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00")) +
      ggplot2::labs(
        title = "Permutation Feature Importance",
        subtitle = sprintf("Based on change in '%s' on the validation set", metric_name),
        x = paste("Importance (Performance Drop)", if(mode == "min") "(Higher is worse)" else ""),
        y = "Feature"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 12, color = "gray30"),
        axis.title = ggplot2::element_text(face = "bold")
      )

    print(p)
    return(invisible(importance_df))
  }

  return(importance_df)
}
