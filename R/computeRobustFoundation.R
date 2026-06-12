#' @title Compute Robust Foundation (Internal)
#'
#' @description
#' Implements Stage 1 of the FSCRE algorithm: separate data cleaning for X and y,
#' and robust, PSD correlation estimation.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param x_preprocess Method for X cleaning ("ddc", "none").
#' @param y_preprocess Method for y cleaning ("wrap", "robust_z", "none").
#' @param cor_estimator Method for correlations ("wrap", "pearson").
#'
#' @return A list containing the imputed data (x.imp, y.imp), the correlation structures (Rx, ry),
#' and the DDC object for prediction.
#'
#' @keywords internal
#'
#' @importFrom stats cor median mad sd
#' @importFrom cellWise DDC wrap
#'
computeRobustFoundation <- function(x, y, x_preprocess, y_preprocess, cor_estimator) {

    ddc_object <- NULL # Default
    p <- ncol(x)

    if (is.null(colnames(x))) {
        colnames(x) <- paste0("V", 1:ncol(x))
    }

    # _______________________________________
    # 1. X Preprocessing (Independent of y)
    # _______________________________________

    # Prevents Target Leakage
    if (x_preprocess == "ddc") {
        # Run Fast DDC on predictors ONLY
        ddc_out <- cellWise::DDC(x, DDCpars = list(fastDDC = TRUE, silent = TRUE))
        ddc_object <- ddc_out
        x.imp <- ddc_out$Ximp
    } else {
        x.imp <- x
    }

    # ________________________________
    # 2. y Preprocessing (Univariate)
    # ________________________________

    if (y_preprocess == "wrap") {
        # cellWise::wrap expects a matrix
        yw_out <- cellWise::wrap(as.matrix(y))
        y.imp <- as.numeric(yw_out$Xw)
    } else if (y_preprocess == "robust_z") {
        med_y <- median(y)
        mad_y <- mad(y)
        if (mad_y == 0) mad_y <- sd(y)
        # Simple winsorization at 3 MADs
        y.imp <- pmin(pmax(y, med_y - 3 * mad_y), med_y + 3 * mad_y)
    } else {
        y.imp <- y
    }

    # ___________________________________________
    # 3. Correlation Estimation (PSD guaranteed)
    # ___________________________________________

    if (cor_estimator == "wrap") {
        # Raymaekers & Rousseeuw (2021): wrap then Pearson cor -> Guaranteed PSD
        # Note: We wrap the already cleaned data to ensure maximum robustness
        Xw <- cellWise::wrap(x.imp)$Xw
        yw <- cellWise::wrap(as.matrix(y.imp))$Xw
        
        Rx <- cor(Xw)
        ry <- as.numeric(cor(Xw, yw))
    } else {
        # Standard Pearson on imputed data
        Rx <- cor(x.imp)
        ry <- as.numeric(cor(x.imp, y.imp))
    }

    return(list(x.imp = x.imp,
                y.imp = y.imp,
                Rx = Rx,
                ry = ry,
                ddc.object = ddc_object))
}