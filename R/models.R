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
#' * *`susie`* - SuSiE regression using *`susieR`* package
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
#' * *`metrics`* - Goodness of fit measures
#' * *`coefficients`* - Fitted coefficients
#'
#' @export
#' @examples
#' data("example_matrix")
#' df <- as.data.frame(example_matrix)
#' fit_model(g1 ~ ., data = df, method = "srm")
#' fit_model(g1 ~ ., data = df, method = "glm")
#' fit_model(g1 ~ ., data = df, method = "glmnet")
#' fit_model(g1 ~ ., data = df, method = "cv.glmnet")
#' fit_model(g1 ~ ., data = df, method = "brms")
#' fit_model(g1 ~ ., data = df, method = "xgb")
#' fit_model(g1 ~ ., data = df, method = "susie")
fit_model <- function(
    formula,
    data,
    method = c(
      "srm",
      "cvsrm",
      "glm",
      "glmnet",
      "cv.glmnet",
      "brms",
      "xgb",
      "susie"
    ),
    family = gaussian,
    alpha = 1,
    ...) {
  method <- match.arg(method)
  result <- switch(
    EXPR = method,
    "srm" = fit_srm2(
      formula, data, ...
    ),
    "cvsrm" = fit_cvsrm(
      formula, data, ...
    ),
    "glm" = fit_glm(
      formula, data,
      family = family, ...
    ),
    "glmnet" = fit_glmnet(
      formula, data,
      family = family, alpha = alpha, ...
    ),
    "cv.glmnet" = fit_cvglmnet(
      formula, data,
      family = family, alpha = alpha, ...
    ),
    "brms" = fit_brms(
      formula, data,
      family = family, ...
    ),
    "xgb" = fit_xgb(
      formula, data, ...
    ),
    "susie" = fit_susie(
      formula, data, ...
    )
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
#' @param verbose If `TRUE`, show warning messages.
#'
#' @param ... Additional parameters passed to the underlying sparse regression function.
#'
#' @return A list containing two data frames:
#' * *`metrics`* - Goodness of fit measures
#' * *`coefficients`* - Fitted coefficients with sparse structure
#'
#' @export
#' @examples
#' data("example_matrix")
#' df <- as.data.frame(example_matrix)
#' fit_srm2(g1 ~ ., data = df)
fit_srm2 <- function(
    formula,
    data,
    verbose = TRUE,
    ...) {
  model_mat <- stats::model.matrix(
    formula,
    data = data
  )[, -1]
  response <- data[[formula[[2]]]]

  result <- fit_srm(
    x = model_mat,
    y = response,
    verbose = verbose,
    ...
  )

  result$metrics <- tibble::tibble(
    r_squared = result$metrics$r_squared
  )

  result$coefficients <- tibble::tibble(
    variable = result$coefficients$variable,
    coefficient = result$coefficients$coefficient
  )

  return(result)
}

#' @title Cross-validation for sparse regression models
#'
#' @description
#' Fits a sparse regression model with cross-validation.
#'
#' @md
#' @inheritParams fit_srm2
#'
#' @return A list containing two data frames:
#' * `metrics` - Goodness of fit measures
#' * `coefficients` - Fitted coefficients with sparse structure
#'
#' @export
#' @examples
#' data("example_matrix")
#' df <- as.data.frame(example_matrix)
#' fit_cvsrm(g1 ~ ., data = df)
fit_cvsrm <- function(
    formula,
    data,
    ...) {
  model_mat <- stats::model.matrix(
    formula,
    data = data
  )[, -1]
  response <- data[[formula[[2]]]]

  result <- fit_srm(
    x = model_mat,
    y = response,
    cross_validation = TRUE,
    ...
  )
  result$metrics <- tibble::tibble(
    r_squared = result$metrics$r_squared
  )
  result$coefficients <- tibble::tibble(
    variable = result$coefficients$variable,
    coefficient = result$coefficients$coefficient
  )

  return(result)
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
#' * *`metrics`* - Goodness of fit measures including R-squared
#' * *`coefficients`* - Fitted coefficients with standard errors and p-values
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
  metrics <- tibble::tibble(
    r_squared = with(s, 1 - deviance / null.deviance)
  )
  coefficients <- tibble::as_tibble(
    s$coefficients,
    rownames = "variable"
  )
  colnames(coefficients) <- c(
    "variable", "coefficient", "std_err", "statistic", "pval"
  )
  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
}

#' @title Fit regularized generalized linear model
#'
#' @description
#' Fits a regularized generalized linear model using elastic net penalties via glmnet.
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
#' * *`metrics`* - Goodness of fit measures including lambda and R-squared
#' * *`coefficients`* - Regularized coefficient estimates
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
  metrics <- tibble::tibble(
    lambda = lambda_choose,
    r_squared = fit$dev.ratio[which_max],
    alpha = alpha
  )
  coefficients <- tibble::as_tibble(
    as.matrix(coef(fit, s = lambda_choose)),
    rownames = "variable"
  )
  colnames(coefficients) <- c("variable", "coefficient")

  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
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
#' * *`metrics`* - Cross-validated performance metrics
#' * *`coefficients`* - Coefficients selected by cross-validation
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
  metrics <- tibble::tibble(
    lambda = fit$lambda.1se,
    r_squared = fit$glmnet.fit$dev.ratio[which_max],
    alpha = alpha
  )
  coefficients <- tibble::as_tibble(
    as.matrix(coef(fit)),
    rownames = "variable"
  )
  colnames(coefficients) <- c("variable", "coefficient")
  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
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
#' * *`metrics`* - Bayesian R-squared and model fit metrics
#' * *`coefficients`* - Posterior estimates with uncertainty intervals
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
  metrics <- tibble::tibble(
    r_squared = as.matrix(brms::bayes_R2(fit))[, "Estimate"]
  )
  coefficients <- tibble::as_tibble(
    as.matrix(
      brms::fixef(fit, probs = c(0.05, 0.95))
    ),
    rownames = "variable"
  )
  colnames(coefficients) <- c(
    "variable", "coefficient", "est_error", "q5", "q95"
  )
  coefficients$pval <- bayestestR::p_map(fit)$p_MAP
  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
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
#' * *`metrics`* - Model performance metrics
#' * *`coefficients`* - Feature importance measures
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
  metrics <- tibble::tibble(
    r_squared = r_square(response, y_pred)
  )
  coefficients <- tibble::as_tibble(
    as.data.frame(
      xgboost::xgb.importance(
        model = fit
      )
    )
  )
  colnames(coefficients) <- c(
    "variable", "gain", "cover", "frequency"
  )

  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
}

#' @title Fit a SuSiE regression model
#'
#' @description
#' Fits a Sum of Single Effects (SuSiE) regression model using the susieR package.
#'
#' @md
#' @param formula An object of class *`formula`* with a symbolic description
#' of the model to be fitted.
#'
#' @param data A *`data.frame`* containing the variables in the model.
#'
#' @param ... Additional parameters passed to *`susieR::susie`*.
#'
#' @return A list containing two data frames:
#' * *`metrics`* - Goodness of fit measures
#' * *`coefficients`* - Fitted coefficients
#'
#' @export
fit_susie <- function(
    formula,
    data,
    ...) {
  model_mat <- stats::model.matrix(
    formula,
    data = data
  )[, -1]
  response <- data[[formula[[2]]]]

  fit <- susieR::susie(
    X = model_mat,
    y = response,
    ...
  )

  y_pred <- predict(fit)
  r2 <- 1 - sum((response - y_pred)^2) / sum((response - mean(response))^2)

  metrics <- tibble::tibble(
    r_squared = r2
  )

  coefs <- coef(fit)[-1]
  coefficients <- tibble::tibble(
    variable = colnames(model_mat),
    coefficient = coefs
  )

  return(
    list(
      model = fit,
      metrics = metrics,
      coefficients = coefficients
    )
  )
}
