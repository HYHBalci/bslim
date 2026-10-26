# Fixed external-score orthogonalisation is exactly the ordinary linear sampler
# applied to Z - pi. These are small native integration tests, not convergence
# or frequentist-coverage tests.
.fixed_ortho_data <- function() {
    n <- 45L
    u <- seq(-1.8, 1.8, length.out = n)
    X <- cbind(x1 = u, x2 = sin(2.1 * u), x3 = cos(1.3 * u))
    pi <- plogis(0.2 + 0.5 * X[, 1] - 0.3 * X[, 2])
    z <- as.numeric(seq_len(n) %% 3L != 0L)
    y <- 2 + X[, 1]^2 + (z - pi) * (0.7 + 0.3 * X[, 2]) +
        0.2 * sin(seq_len(n))
    train <- seq_len(36L)
    test <- 37L:n
    list(X_train = X[train, ], Z_train = z[train], y_train = y[train],
         propensity_train = pi[train], X_test = X[test, ],
         Z_test = z[test], propensity_test = pi[test])
}

.fixed_ortho_args <- function(num_gfr = 2L, seed = 1709L) {
    c(.fixed_ortho_data(), list(
        num_gfr = num_gfr, num_burnin = 3L, num_mcmc = 5L,
        general_params = list(
            random_seed = seed, verbose = FALSE, standardize = TRUE,
            standardize_cov = FALSE, sample_global_prior = "hc-hs",
            unlink = TRUE, gibbs = TRUE, hn_scale = 1, p0 = 0,
            orthogonalize = TRUE, propensity_covariate = "none",
            propensity_seperate = "none", adaptive_coding = FALSE,
            sample_sigma2_global = TRUE, sigma_residual = 0,
            interaction_rule = "continuous", save_output = TRUE,
            regularize_ATE = FALSE, simple_prior = FALSE, use_ncp = FALSE,
            robust = FALSE, keep_gfr = FALSE, keep_burnin = FALSE,
            keep_every = 1L, num_chains = 1L),
        prognostic_forest_params = list(num_trees = 3L, min_samples_leaf = 3L,
                                       max_depth = 3L),
        treatment_effect_forest_params = list(num_trees = 0L)))
}

.fixed_ortho_fit <- function(args = .fixed_ortho_args()) {
    # The legacy native sampler writes a debug log in the working directory.
    # Keep that generated file in this R session's temporary directory.
    withr::local_dir(tempdir())
    do.call(bcf_linear_probit, args)
}

test_that("fixed orthogonalisation requires valid external training scores", {
    args <- .fixed_ortho_args()
    bad_scores <- list(NULL, numeric(0), args$propensity_train[-1L],
                       c(NA_real_, args$propensity_train[-1L]),
                       c(NaN, args$propensity_train[-1L]),
                       c(Inf, args$propensity_train[-1L]),
                       c(-0.01, args$propensity_train[-1L]),
                       c(1.01, args$propensity_train[-1L]),
                       cbind(args$propensity_train, args$propensity_train),
                       as.character(args$propensity_train))
    for (bad in bad_scores) {
        invalid <- args
        invalid["propensity_train"] <- list(bad)
        expect_error(.fixed_ortho_fit(invalid), "propensity_train")
    }
})

test_that("fixed orthogonalisation requires valid external test scores", {
    args <- .fixed_ortho_args()
    bad_scores <- list(NULL, numeric(0), args$propensity_test[-1L],
                       c(NA_real_, args$propensity_test[-1L]),
                       c(Inf, args$propensity_test[-1L]),
                       c(-0.01, args$propensity_test[-1L]),
                       c(1.01, args$propensity_test[-1L]),
                       cbind(args$propensity_test, args$propensity_test),
                       as.character(args$propensity_test))
    for (bad in bad_scores) {
        invalid <- args
        invalid["propensity_test"] <- list(bad)
        expect_error(.fixed_ortho_fit(invalid), "propensity_test")
    }
})

test_that("fixed orthogonalisation requires raw binary treatment", {
    args <- .fixed_ortho_args()
    for (field in c("Z_train", "Z_test")) {
        bad_treatments <- list(args[[field]] - 0.5, args[[field]][-1L],
                               c(NA_real_, args[[field]][-1L]),
                               c(Inf, args[[field]][-1L]),
                               cbind(args[[field]], args[[field]]))
        for (bad in bad_treatments) {
            invalid <- args
            invalid[[field]] <- bad
            expect_error(.fixed_ortho_fit(invalid), field)
        }
    }
})

test_that("fixed orthogonalisation matches manual residualisation in GFR and MCMC", {
    for (gfr in c(0L, 2L)) {
        args <- .fixed_ortho_args(num_gfr = gfr)
        original_pi_train <- args$propensity_train
        original_pi_test <- args$propensity_test
        orthogonal <- .fixed_ortho_fit(args)
        manual_args <- args
        manual_args$Z_train <- args$Z_train - args$propensity_train
        # The ordinary linear sampler's historical test-set return uses an
        # inactive tau forest. Fit training rows only, then reconstruct held-out
        # predictions using its coefficient samples and prognostic forest.
        manual_args$X_test <- NULL
        manual_args$Z_test <- NULL
        manual_args$propensity_test <- NULL
        # Supplied train scores are unused when both propensity settings are
        # none, but avoid a legacy NULL-to-native conversion in no-GFR fits.
        manual_args$general_params$orthogonalize <- FALSE
        manual <- .fixed_ortho_fit(manual_args)

        # A hidden propensity update would consume RNG and break this equality.
        for (field in c("alpha", "Beta", "Beta_int", "Tau", "Tau_int",
                        "Tau_glob", "sigma2_global_samples", "mu_hat_train")) {
            expect_equal(orthogonal[[field]], manual[[field]], tolerance = 1e-12,
                         info = paste("GFR", gfr, field))
        }
        # Reconstruct the ordinary sampler's linear predictions from every
        # coefficient draw; its historical cached tau fields are last-draw only.
        manual_train <- predict_linear_bcf(manual, args$X_train,
                                            manual_args$Z_train)
        manual_test <- predict_linear_bcf(manual, args$X_test,
                                           args$Z_test - args$propensity_test)
        expect_equal(orthogonal$tau_hat_train, manual_train$tau_hat,
                     tolerance = 1e-12)
        expect_equal(orthogonal$y_hat_train, manual_train$y_hat,
                     tolerance = 1e-12)
        expect_equal(orthogonal$tau_hat_test, manual_test$tau_hat,
                     tolerance = 1e-12)
        expect_equal(orthogonal$mu_hat_test, manual_test$mu_hat,
                     tolerance = 1e-12)
        expect_equal(orthogonal$y_hat_test, manual_test$y_hat,
                     tolerance = 1e-12)
        expect_identical(args$propensity_train, original_pi_train)
        expect_identical(args$propensity_test, original_pi_test)
        expect_true(orthogonal$model_params$orthogonalize)
        expect_identical(orthogonal$model_params$orthogonalization, "fixed-external")
        expect_true(orthogonal$model_params$hchs_variance_scope)
        expect_identical(orthogonal$model_params$variance_update,
                         "hchs-full-conditional-20260927")
        expect_identical(orthogonal$model_params$propensity_covariate, "none")
        expect_false(orthogonal$model_params$internal_propensity_model)
        expect_false(orthogonal$model_params$adaptive_coding)
        expect_null(orthogonal$theta_samples)
        expect_null(orthogonal$bart_propensity_model)
        expect_identical(as.numeric(orthogonal$propensity_train), original_pi_train)
        expect_identical(as.numeric(orthogonal$propensity_test), original_pi_test)
        expect_identical(as.numeric(orthogonal$treatment_coding_train),
                         args$Z_train - original_pi_train)
        expect_identical(as.numeric(orthogonal$treatment_coding_test),
                         args$Z_test - original_pi_test)
        expect_equal(orthogonal$model_params$num_prognostic_covariates,
                     ncol(args$X_train))
        expect_true(all(is.finite(orthogonal$sigma2_global_samples)))
        expect_true(all(orthogonal$sigma2_global_samples > 0))
        expect_true(all(is.finite(orthogonal$Tau_glob)))
        expect_true(all(orthogonal$Tau_glob > 0))
        expect_equal(dim(orthogonal$y_hat_train), c(36L, 5L))
        expect_equal(dim(orthogonal$y_hat_test), c(9L, 5L))
        expect_true(all(is.finite(orthogonal$y_hat_train)))
        expect_true(all(is.finite(orthogonal$y_hat_test)))

        train_prediction <- predict_linear_bcf(orthogonal, args$X_train,
                                               args$Z_train, original_pi_train)
        test_prediction <- predict_linear_bcf(orthogonal, args$X_test,
                                              args$Z_test, original_pi_test)
        expect_equal(orthogonal$y_hat_train,
                     train_prediction$mu_hat + train_prediction$tau_hat *
                         (args$Z_train - original_pi_train), tolerance = 1e-12)
        expect_equal(orthogonal$y_hat_test,
                     test_prediction$mu_hat + test_prediction$tau_hat *
                         (args$Z_test - original_pi_test), tolerance = 1e-12)
        expect_equal(train_prediction$y_hat, orthogonal$y_hat_train,
                     tolerance = 1e-12)
        expect_equal(test_prediction$y_hat, orthogonal$y_hat_test,
                     tolerance = 1e-12)
    }
})

test_that("orthogonalisation excludes propensity from mu even with default options", {
    args <- .fixed_ortho_args()
    baseline <- .fixed_ortho_fit(args)
    # Defaults request adaptive coding and propensity in mu. Fixed external
    # residualisation takes precedence, rather than silently combining codings.
    args$general_params$propensity_covariate <- NULL
    args$general_params$adaptive_coding <- NULL
    fit <- .fixed_ortho_fit(args)
    expect_identical(fit$model_params$propensity_covariate, "none")
    expect_false(fit$model_params$adaptive_coding)
    expect_equal(fit$Beta, baseline$Beta, tolerance = 1e-12)
    expect_equal(fit$mu_hat_train, baseline$mu_hat_train, tolerance = 1e-12)
    prediction <- predict_linear_bcf(fit, args$X_test, args$Z_test,
                                     args$propensity_test)
    changed_pi <- pmin(0.99, args$propensity_test + 0.1)
    changed <- predict_linear_bcf(fit, args$X_test, args$Z_test, changed_pi)
    expect_equal(changed$mu_hat, prediction$mu_hat, tolerance = 0)
    expect_equal(changed$tau_hat, prediction$tau_hat, tolerance = 0)
    expect_equal(changed$y_hat - prediction$y_hat,
                 -(changed_pi - args$propensity_test) * prediction$tau_hat,
                 tolerance = 1e-12)
})

test_that("new-data orthogonal predictions validate external propensity", {
    args <- .fixed_ortho_args()
    fit <- .fixed_ortho_fit(args)
    bad_scores <- list(NULL, numeric(0), args$propensity_test[-1L],
                       c(NA_real_, args$propensity_test[-1L]),
                       c(NaN, args$propensity_test[-1L]),
                       c(Inf, args$propensity_test[-1L]),
                       c(-0.01, args$propensity_test[-1L]),
                       c(1.01, args$propensity_test[-1L]),
                       cbind(args$propensity_test, args$propensity_test),
                       as.character(args$propensity_test))
    for (bad in bad_scores) {
        expect_error(predict_linear_bcf(fit, args$X_test, args$Z_test, bad),
                     "[Pp]ropensity")
    }
    endpoint_scores <- rep(c(0, 1), length.out = nrow(args$X_test))
    endpoint_prediction <- predict_linear_bcf(fit, args$X_test, args$Z_test,
                                              endpoint_scores)
    expect_equal(endpoint_prediction$y_hat,
                 endpoint_prediction$mu_hat + endpoint_prediction$tau_hat *
                     (args$Z_test - endpoint_scores), tolerance = 1e-12)
    expect_error(predict_linear_bcf(fit, args$X_test, args$Z_test - 0.5,
                                   args$propensity_test), "Z")
    expect_error(predict_linear_bcf(fit, args$X_test, args$Z_test[-1L],
                                   args$propensity_test), "Z")
    expect_equal(predict.bslim_bcfmodel(fit, args$X_test, args$Z_test,
                                 args$propensity_test),
                 predict_linear_bcf(fit, args$X_test, args$Z_test,
                                    args$propensity_test), tolerance = 1e-12)
})

test_that("fixed orthogonalisation retains aligned GFR predictions", {
    for (mcmc in c(0L, 3L)) {
        args <- .fixed_ortho_args(num_gfr = 2L)
        args$num_burnin <- 0L
        args$num_mcmc <- mcmc
        args$general_params$keep_gfr <- TRUE
        # A one-column matrix is accepted as well as a numeric vector.
        args$propensity_train <- as.matrix(args$propensity_train)
        args$propensity_test <- as.matrix(args$propensity_test)
        fit <- .fixed_ortho_fit(args)
        prediction <- predict_linear_bcf(fit, args$X_test, args$Z_test,
                                         args$propensity_test)
        expect_equal(dim(fit$y_hat_test), c(9L, 2L + mcmc))
        expect_equal(prediction$y_hat, fit$y_hat_test, tolerance = 1e-12)
        expect_equal(prediction$tau_hat, fit$tau_hat_test, tolerance = 1e-12)
        expect_equal(fit$y_hat_test,
                     fit$mu_hat_test + fit$tau_hat_test *
                         (args$Z_test - as.numeric(args$propensity_test)),
                     tolerance = 1e-12)
    }
})

test_that("fixed orthogonalisation can fit without evaluation observations", {
    args <- .fixed_ortho_args()
    args$X_test <- NULL
    args$Z_test <- NULL
    args$propensity_test <- NULL
    fit <- .fixed_ortho_fit(args)
    expect_null(fit$y_hat_test)
    expect_null(fit$propensity_test)
    expect_null(fit$treatment_coding_test)
    expect_equal(dim(fit$y_hat_train), c(36L, 5L))
})

test_that("fixed orthogonalisation rejects unsupported output/transform modes", {
    args <- .fixed_ortho_args()
    for (unsupported in list(TRUE, 1, NA)) {
        args$general_params$standardize_cov <- unsupported
        expect_error(.fixed_ortho_fit(args), "standardize_cov")
    }
    args$general_params$standardize_cov <- FALSE
    for (unsupported in list(TRUE, 1, NA)) {
        args$general_params$keep_burnin <- unsupported
        expect_error(.fixed_ortho_fit(args), "keep_burnin")
    }
})

test_that("fixed orthogonalisation aligns multiple thinned chains without debug output", {
    args <- .fixed_ortho_args(num_gfr = 2L)
    args$general_params$num_chains <- 2L
    args$general_params$keep_every <- 2L
    args$general_params$save_output <- FALSE
    fit <- .fixed_ortho_fit(args)
    expect_equal(dim(fit$orthogonal_coefficients), c(10L, 7L))
    expected_coefficients <- cbind(as.vector(t(fit$alpha)),
        matrix(aperm(fit$Beta, c(2, 1, 3)), nrow = 10L),
        matrix(aperm(fit$Beta_int, c(2, 1, 3)), nrow = 10L))
    expect_equal(unname(fit$orthogonal_coefficients), unname(expected_coefficients),
                 tolerance = 0)
    expect_equal(dim(fit$y_hat_test), c(9L, 10L))
    expect_true(all(is.finite(fit$orthogonal_coefficients)))
    prediction <- predict_linear_bcf(fit, args$X_test, args$Z_test,
                                     args$propensity_test)
    expect_equal(prediction$y_hat, fit$y_hat_test, tolerance = 1e-12)
    expect_equal(prediction$tau_hat, fit$tau_hat_test, tolerance = 1e-12)
    treated <- predict_linear_bcf(fit, args$X_test, rep(1, nrow(args$X_test)),
                                  args$propensity_test)
    control <- predict_linear_bcf(fit, args$X_test, rep(0, nrow(args$X_test)),
                                  args$propensity_test)
    expect_equal(treated$y_hat - control$y_hat, prediction$tau_hat,
                 tolerance = 1e-12)
    expect_null(fit$theta_samples)
})

test_that("fixed orthogonalisation accepts endpoint scores without outcome standardisation", {
    args <- .fixed_ortho_args()
    args$general_params$standardize <- FALSE
    args$propensity_train <- rep(c(0, 1), length.out = nrow(args$X_train))
    args$propensity_test <- rep(c(1, 0), length.out = nrow(args$X_test))
    fit <- .fixed_ortho_fit(args)
    repeated <- .fixed_ortho_fit(args)
    other_args <- args
    other_args$general_params$random_seed <- 1710L
    other <- .fixed_ortho_fit(other_args)
    expect_identical(fit$orthogonal_coefficients, repeated$orthogonal_coefficients)
    expect_false(identical(fit$orthogonal_coefficients, other$orthogonal_coefficients))
    expect_identical(fit$propensity_train, args$propensity_train)
    expect_identical(fit$propensity_test, args$propensity_test)
    expect_equal(fit$model_params$outcome_scale, 1)
    prediction <- predict_linear_bcf(fit, args$X_test, args$Z_test,
                                     args$propensity_test)
    expect_equal(prediction$y_hat, fit$y_hat_test, tolerance = 1e-12)
    expect_equal(prediction$y_hat,
                 prediction$mu_hat + prediction$tau_hat *
                     (args$Z_test - args$propensity_test), tolerance = 1e-12)
    expect_true(all(is.finite(fit$sigma2_global_samples)))
    expect_true(all(fit$sigma2_global_samples > 0))
})
