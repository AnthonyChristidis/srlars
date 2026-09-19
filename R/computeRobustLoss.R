#' @title Compute Robust Loss (Internal)
#'
#' @description
#' Scores validation residuals with the requested loss (mse, trimmed, or huber). Extracted
#' from \code{computeCVError} so the same loss formulas can be reused by the outer
#' cross-validation loop in \code{cv.srlars}.
#'
#' @param r_val Numeric vector of validation residuals.
#' @param r_train Numeric vector of training residuals, used only to set the scale for the
#'   Huber loss. Ignored for \code{cv_loss = "mse"} or \code{"trimmed"}.
#' @param cv_loss Character. Loss function: "huber", "trimmed", or "mse".
#'
#' @return Numeric. The scalar loss value.
#'
#' @keywords internal
#'
#' @importFrom stats mad quantile sd
computeRobustLoss <- function(r_val, r_train, cv_loss) {

    if (cv_loss == "mse") {
        return(mean(r_val^2))

    } else if (cv_loss == "trimmed") {
        r2 <- r_val^2
        cutoff <- quantile(r2, 0.90, names = FALSE)
        return(mean(r2[r2 <= cutoff]))

    } else if (cv_loss == "huber") {
        k <- 1.345
        scale_est <- mad(r_train)
        if (scale_est < 1e-6) scale_est <- max(sd(r_train), 1e-6)

        u <- r_val / scale_est
        huber_loss <- ifelse(abs(u) <= k,
                             0.5 * u^2,
                             k * abs(u) - 0.5 * k^2)
        return(mean(huber_loss) * (scale_est^2))
    }
}
