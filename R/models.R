#' @title Fit model
#'
#' @description
#' Fits various types of regression models including generalized linear models,
#' regularized models, Bayesian models, and machine learning approaches.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param method A character string indicating the method to fit the model.
#' This can take any of the following choices:
#'
#' * *`srm`* - Sparse Regression Model
#'
#' * *`glm`* - Generalized Linear Model using *`stats::glm`*
#'
#' * *`glmnet`*, *`cv.glmnet`* - Regularized GLM using *`glmnet::glmnet`*
#'
#' * *`brms`* - Bayesian Regression using *`brms`* package
#'
#' * *`xgb`* - Gradient Boosting using *`xgboost`*
#'
#' * *`bagging_ridge`* - Bagging Ridge Regression via scikit-learn
#'
#' * *`bayesian_ridge`* - Bayesian Ridge Regression via scikit-learn
#'
#' @param family A description of the error distribution and link function.
#' See *`stats::family`* for details.
#'
#' @param alpha The elasticnet mixing parameter, between 0 and 1.
#' See *`glmnet::glmnet`* for details.
#'
#' @param ... Additional parameters passed to the underlying model fitting function.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Goodness of fit measures
#' * *`coefs`* - Fitted coefficients
#'
#' @export
fit_model <- function(
    formula,
    data,
    method = c(
      "srm",
      "glm",
      "glmnet",
      "cv.glmnet",
      "brms",
      "xgb",
      "bagging_ridge",
      "bayesian_ridge"
    ),
    family = gaussian,
    alpha = 1,
    ...
  ) {
  # Match args
  method <- match.arg(method)
  result <- switch(method,
    "srm" = fit_srm(formula, data, ...),
    "glm" = fit_glm(formula, data, family = family, ...),
    "glmnet" = fit_glmnet(formula, data, family = family, alpha = alpha, ...),
    "cv.glmnet" = fit_cvglmnet(formula, data, family = family, alpha = alpha, ...),
    "brms" = fit_brms(formula, data, family = family, ...),
    "xgb" = fit_xgb(formula, data, ...),
    "bagging_ridge" = fit_bagging_ridge(formula, data, alpha = alpha, ...),
    "bayesian_ridge" = fit_bayesian_ridge(formula, data, ...)
  )

  return(result)
}

#' @title Fit a sparse regression model
#'
#' @description
#' Fits a sparse regression model using specialized algorithms for high-dimensional data.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param ... Additional parameters passed to the underlying sparse regression function.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Goodness of fit measures
#' * *`coefs`* - Fitted coefficients with sparse structure
#'
#' @export
fit_srm <- function(
    formula,
    data,
    ...) {
  model_mat <- stats::model.matrix(
    formula,
    data = data
  )[, -1]
  response <- data[[formula[[2]]]]

  fit <- fit_sparse_regression(
    model_mat,
    response,
    ...
  )
  fit_inf <- print(fit)
  lambda <- fit_inf$lambda[which.max(fit_inf$suppSize)]
  gamma <- fit_inf$gamma[which.max(fit_inf$suppSize)]

  y_pred <- predict(
    fit,
    newx = model_mat,
    gamma = gamma,
    lambda = lambda
  )
  gof <- tibble::tibble(
    rsq = r2(response, y_pred)
  )
  coefs <- as.vector(
    coef(
      fit,
      lambda = lambda,
      gamma = gamma
    )
  )[-1]
  coefs <- normalization(coefs, method = "sum")
  if (length(coefs) != ncol(model_mat)) {
    coefs <- rep(0, ncol(model_mat))
  }
  coefs <- tibble::as_tibble(
    cbind(
      as.data.frame(colnames(model_mat)),
      as.data.frame(coefs)
    )
  )
  colnames(coefs) <- c("term", "estimate")

  return(list(gof = gof, coefs = coefs))
}

#' @title Cross-validation for sparse regression models
#'
#' @description
#' Performs cross-validation for sparse regression models to select optimal parameters.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param ... Additional parameters passed to the cross-validation function.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Cross-validated goodness of fit measures
#' * *`coefs`* - Coefficients selected through cross-validation
#'
#' @export
fit_cvsrm <- function(
    formula,
    data,
    ...) {
  model_mat <- stats::model.matrix(formula, data = data)[, -1]
  response <- data[[formula[[2]]]]

  regulators_num <- ncol(model_mat)
  fit <- try(
    fit_sparse_regression(
      model_mat, response,
      cross_validation = TRUE,
      seed = 1,
      n_folds = 10
    )
  )

  if (any(class(fit) == "try-error")) {
    if (verbose) message("Cross validation error, used fit instead.")
    fit <- try(
      fit_sparse_regression(
        model_mat, response,
        penalty = penalty,
        algorithm = algorithm,
        regulators_num = regulators_num
      )
    )
    if (any(class(fit) == "try-error")) {
      return(rep(0, ncol(model_mat)))
    }
    fit_inf <- print(fit)
    lambda <- fit_inf$lambda[which.max(fit_inf$suppSize)]
    gamma <- fit_inf$gamma[which.max(fit_inf$suppSize)]
  } else {
    gamma <- fit$fit$gamma[which(
      unlist(lapply(
        fit$cvMeans, min
      )) == min(unlist(lapply(fit$cvMeans, min)))
    )]
    lambda_list <- dplyr::filter(print(fit), gamma == gamma)
    if (regulators_num %in% lambda_list$regulators_num) {
      lambda <- lambda_list$regulators_num[which(
        lambda_list$regulators_num == regulators_num
      )]
    } else {
      lambda <- min(lambda_list$lambda)
    }
  }

  y_pred <- predict(
    fit,
    newx = model_mat,
    gamma = gamma,
    lambda = lambda
  )
  gof <- tibble::tibble(
    rsq = r2(response, y_pred)
  )
  coefs <- as.vector(
    coef(
      fit,
      lambda = lambda,
      gamma = gamma
    )
  )[-1]
  coefs <- normalization(coefs, method = "sum")
  if (length(coefs) != ncol(model_mat)) {
    coefs <- rep(0, ncol(model_mat))
  }
  coefs <- tibble::as_tibble(
    cbind(
      as.data.frame(colnames(model_mat)),
      as.data.frame(coefs)
    )
  )

  colnames(coefs) <- c("term", "estimate")

  return(list(gof = gof, coefs = coefs))
}

#' @title Fit generalized linear model
#'
#' @description
#' Fits a generalized linear model using standard maximum likelihood estimation.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param family A description of the error distribution and link function.
#' See *`stats::family`* for details.
#'
#' @param ... Additional parameters passed to *`stats::glm`*.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Goodness of fit measures including R-squared
#' * *`coefs`* - Fitted coefficients with standard errors and p-values
#'
#' @export
fit_glm <- function(
    formula,
    data,
    family = gaussian,
    ...) {
  fit <- suppressWarnings(
    stats::glm(
      formula,
      data = data,
      family = family,
      ...
    )
  )
  s <- summary(fit)
  gof <- tibble::tibble(
    rsq = with(s, 1 - deviance / null.deviance)
  )
  coefs <- tibble::as_tibble(s$coefficients, rownames = "term")
  colnames(coefs) <- c("term", "estimate", "std_err", "statistic", "pval")
  return(list(gof = gof, coefs = coefs))
}

#' @title Fit regularized generalized linear model
#'
#' @description
#' Fits a regularized GLM using elastic net penalties via glmnet.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param family A description of the error distribution and link function.
#' See *`stats::family`* for details.
#'
#' @param alpha The elasticnet mixing parameter, between 0 and 1:
#' * *`alpha = 1`*: Lasso penalty
#' * *`alpha = 0`*: Ridge penalty
#' * *`0 < alpha < 1`*: Elastic net penalty
#'
#' @param ... Additional parameters passed to *`glmnet::glmnet`*.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Goodness of fit measures including lambda and R-squared
#' * *`coefs`* - Regularized coefficient estimates
#'
#' @export
fit_glmnet <- function(
    formula,
    data,
    family = gaussian,
    alpha = 0.5,
    ...) {
  fit <- glmnetUtils::glmnet(
    formula,
    data = data,
    family = family,
    alpha = alpha,
    ...
  )
  class(fit) <- "glmnet"
  which_max <- which(fit$dev.ratio > max(fit$dev.ratio) * 0.95)[1]
  lambda_choose <- fit$lambda[which_max]
  gof <- tibble::tibble(
    lambda = lambda_choose,
    rsq = fit$dev.ratio[which_max],
    alpha = alpha
  )
  coefs <- tibble::as_tibble(
    as.matrix(coef(fit, s = lambda_choose)),
    rownames = "term"
  )
  colnames(coefs) <- c("term", "estimate")

  return(list(gof = gof, coefs = coefs))
}

#' @title Cross-validation for regularized generalized linear models
#'
#' @description
#' Performs cross-validation to select optimal regularization parameters for GLMnet models.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param family A description of the error distribution and link function.
#' See *`stats::family`* for details.
#'
#' @param alpha The elasticnet mixing parameter, between 0 and 1.
#' See *`glmnet::glmnet`* for details.
#'
#' @param ... Additional parameters passed to *`glmnet::cv.glmnet`*.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Cross-validated performance metrics
#' * *`coefs`* - Coefficients selected by cross-validation
#'
#' @export
fit_cvglmnet <- function(
    formula,
    data,
    family = gaussian,
    alpha = 0.5,
    ...) {
  fit <- glmnetUtils::cv.glmnet(
    formula,
    data = data,
    family = family,
    alpha = alpha,
    ...
  )
  class(fit) <- "cv.glmnet"
  which_max <- fit$index["1se", ]
  gof <- tibble::tibble(
    lambda = fit$lambda.1se,
    rsq = fit$glmnet.fit$dev.ratio[which_max],
    alpha = alpha
  )
  coefs <- tibble::as_tibble(
    as.matrix(coef(fit)),
    rownames = "term"
  )
  colnames(coefs) <- c("term", "estimate")
  return(list(gof = gof, coefs = coefs))
}

#' @title Fit a Bayesian regression model with brms and Stan
#'
#' @description
#' Fits a Bayesian regression model using the brms interface to Stan.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param family A description of the error distribution and link function.
#' See *`stats::family`* for details.
#'
#' @param prior The prior distribution specification for coefficients.
#' Default *`normal(0,1)`* provides ridge-like regularization.
#' See *`brms::set_prior`* for details.
#'
#' @param ... Additional parameters passed to *`brms::brm`*.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Bayesian R-squared and model fit metrics
#' * *`coefs`* - Posterior estimates with uncertainty intervals
#'
#' @export
fit_brms <- function(
    formula,
    data,
    family = gaussian,
    prior = brms::prior(normal(0, 1)),
    ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("The brms package is required to use brms bayesian regression models.")
  }
  fit <- suppressMessages(
    brms::brm(
      formula,
      data = data,
      family = family,
      prior = prior,
      silent = TRUE,
      refresh = 0,
      ...
    )
  )
  gof <- tibble::tibble(
    rsq = as.matrix(brms::bayes_R2(fit))[, "Estimate"]
  )
  coefs <- tibble::as_tibble(
    as.matrix(
      brms::fixef(fit, probs = c(0.05, 0.95))
    ),
    rownames = "term"
  )
  colnames(coefs) <- c("term", "estimate", "est_error", "q5", "q95")
  coefs$pval <- bayestestR::p_map(fit)$p_MAP
  return(list(gof = gof, coefs = coefs))
}

#' @title Fit a gradient boosting regression model with XGBoost
#'
#' @description
#' Fits a gradient boosted tree model using the XGBoost implementation.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param params A list of XGBoost parameters including:
#' * *`max_depth`* - Maximum tree depth
#' * *`eta`* - Learning rate
#' * *`objective`* - Loss function to optimize
#'
#' @param nrounds Maximum number of boosting iterations.
#'
#' @param nthread Number of parallel threads used (-1 for all cores).
#'
#' @param ... Additional parameters passed to *`xgboost::xgboost`*.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Model performance metrics
#' * *`coefs`* - Feature importance measures
#'
#' @export
fit_xgb <- function(
    formula,
    data,
    params = list(
      max_depth = 3,
      eta = 0.01,
      objective = "reg:squarederror"
    ),
    nrounds = 1000,
    nthread = -1,
    ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("The xgboost package is required to use xgb models.")
  }
  model_mat <- stats::model.matrix(formula, data = data)
  response <- data[[formula[[2]]]]
  fit <- xgboost::xgboost(
    data = model_mat,
    label = response,
    verbose = 0,
    params = params,
    nrounds = nrounds,
    nthread = nthread,
    ...
  )
  y_pred <- predict(fit, newdata = model_mat)
  gof <- tibble::tibble(
    rsq = r2(response, y_pred)
  )
  coefs <- tibble::as_tibble(
    as.data.frame(
      xgboost::xgb.importance(
        model = fit
      )
    )
  )
  colnames(coefs) <- c("term", "gain", "cover", "frequency")

  return(list(gof = gof, coefs = coefs))
}

#' @title Fit a bagging ridge regression model via scikit-learn
#'
#' @description
#' Fits an ensemble of ridge regression models using scikit-learn's BaggingRegressor.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param alpha Positive float specifying the regularization strength.
#'
#' @param solver Algorithm to use for optimization:
#' * *`auto`* - Automatically chosen based on data
#' * *`svd`* - Singular Value Decomposition
#' * *`cholesky`* - Cholesky decomposition
#' * *`lsqr`* - Least Squares
#' * *`sparse_cg`* - Conjugate Gradient for sparse matrices
#' * *`sag`*, *`saga`* - Stochastic Average Gradient descent
#'
#' @param bagging_number Number of base estimators in ensemble.
#'
#' @param n_jobs Number of parallel jobs (-1 for all cores).
#'
#' @param p_method Method for calculating p-values:
#' * *`t`* - Student's t-test
#' * *`wilcox`* - Wilcoxon test
#'
#' @param ... Additional parameters passed to scikit-learn.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Ensemble model performance metrics
#' * *`coefs`* - Aggregated coefficient estimates with p-values
#'
#' @export
fit_bagging_ridge <- function(
    formula,
    data,
    alpha = 1,
    solver = "auto",
    bagging_number = 200L,
    n_jobs = 1,
    p_method = c("wilcox", "t"),
    ...) {
  p_method <- match.arg(p_method)
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The reticulate package is required to use bagging ridge models.")
  }
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")
  sklearn <- reticulate::import("sklearn")

  model_mat <- stats::model.matrix(formula, data = data)[, -1]
  if (is.null(ncol(model_mat))) {
    stop("The bagging ridge model requires at least two variables.")
  }
  response <- data[[formula[[2]]]]

  model <- sklearn$ensemble$BaggingRegressor(
    base_estimator = sklearn$linear_model$Ridge(
      alpha = alpha,
      solver = solver,
      random_state = as.integer(123),
      ...
    ),
    n_estimators = as.integer(bagging_number),
    bootstrap = TRUE,
    max_features = 0.8,
    n_jobs = as.integer(n_jobs),
    verbose = FALSE
  )
  model <- model$fit(model_mat, response)

  idx_features <- do.call(cbind, model$estimators_features_) + 1
  coefs_features <- sapply(model$estimators_, function(x) x$coef_)
  coefs <- t(sapply(1:bagging_number, function(i) {
    coefs <- stats::setNames(rep(NaN, ncol(model_mat)), colnames(model_mat))
    if (ncol(model_mat) > 2) {
      coefs[idx_features[, i]] <- coefs_features[, i]
    }
    if (ncol(model_mat) == 2) {
      coefs[idx_features[, i]] <- coefs_features[i]
    }
    return(coefs)
  }))
  p <- switch(p_method,
    "t" = apply(
      coefs, 2, function(x) {
        stats::t.test(x[!is.nan(x)])$p.value
      }
    ),
    "wilcox" = apply(
      coefs, 2, function(x) {
        stats::wilcox.test(x[!is.nan(x)])$p.value
      }
    )
  )
  coefs <- tibble::tibble(
    term = colnames(model_mat),
    estimate = colMeans(coefs, na.rm = TRUE),
    pval = p,
    neglog10p = -log10(ifelse(is.na(p), 1, p))
  )
  y_pred <- model_mat %*% matrix(coefs$estimate)
  gof <- tibble::tibble(
    rsq = r2(response, y_pred)
  )

  return(list(gof = gof, coefs = coefs))
}

#' @title Fit a Bayesian ridge regression model via scikit-learn
#'
#' @description
#' Fits a Bayesian ridge regression model using scikit-learn's implementation.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param ... Additional parameters passed to scikit-learn's BayesianRidge.
#'
#' @return A list containing two data frames:
#' * *`gof`* - Model performance metrics
#' * *`coefs`* - Coefficient estimates with uncertainty measures
#'
#' @export
fit_bayesian_ridge <- function(
    formula,
    data,
    ...) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The reticulate package is required to use bayesian ridge models.")
  }
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")
  sklearn <- reticulate::import("sklearn")

  model_mat <- stats::model.matrix(formula, data = data)[, -1]
  if (is.null(ncol(model_mat))) {
    stop("The bayesian ridge model requires at least two variables.")
  }
  response <- data[[formula[[2]]]]

  model <- sklearn$linear_model$BayesianRidge(...)
  model <- model$fit(model_mat, response)

  coefs <- model$coef_
  coef_var <- diag(model$sigma_)
  p <- stats::pnorm(q = 0, mean = abs(coefs), sd = sqrt(coef_var)) * 2
  coefs <- tibble::tibble(
    term = colnames(model_mat),
    estimate = coefs,
    est_variance = coef_var,
    pval = p,
    neglog10p = -log10(ifelse(is.na(p), 1, p))
  )
  gof <- tibble::tibble(
    rsq = model$score(model_mat, response)
  )

  return(list(gof = gof, coefs = coefs))
}
