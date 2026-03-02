#' @title Compute Robust Foundation (Internal)
#'
#' @description
#' Implements Stage 1 of the FSCRE algorithm: data cleaning via DDC and robust correlation estimation.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param robust Logical.
#'
#' @return A list containing the imputed data (x.imp, y.imp), the robust correlation structures (Rx, ry),
#' and the DDC object for prediction.
#'
#' @keywords internal
#'
#' @importFrom stats cor
#' @importFrom cellWise DDC
#'
computeRobustFoundation <- function(x, y, robust) {

    ddc_object <- NULL # Default

    # Ensure x has column names for DDC matching
    if (is.null(colnames(x))) {
        colnames(x) <- paste0("V", 1:ncol(x))
    }

    if (robust) {
        # 1. Form joint data matrix [X, y]
        xy_data <- cbind(x, y)
        p <- ncol(x)

        # Ensure y column has a specific name we can track
        colnames(xy_data)[p + 1] <- "response_y"

        # 2. Run Fast DDC
        # Uses fastDDC=TRUE as per Raymaekers & Rousseeuw (2021) for scalability
        ddc_out <- cellWise::DDC(xy_data, DDCpars = list(fastDDC = TRUE, silent = TRUE))

        # Save the object for use in predict()
        ddc_object <- ddc_out

        # 3. Extract Imputed Data
        xy_imp <- ddc_out$Ximp
        x.imp <- xy_imp[, 1:p, drop = FALSE]
        y.imp <- xy_imp[, p + 1]

        # 4. Compute Sample Correlations on Cleaned Data
        # "Detect-then-Analyze" paradigm
        R_total <- cor(xy_imp)
        Rx <- R_total[1:p, 1:p, drop = FALSE]
        ry <- R_total[1:p, p + 1]

    } else {
        # Non-Robust Case
        x.imp <- x
        y.imp <- y

        # Compute standard correlations
        Rx <- cor(x)
        ry <- cor(x, y)
        # Ensure ry is a vector
        ry <- as.numeric(ry)
    }

    return(list(x.imp = x.imp,
                y.imp = y.imp,
                Rx = Rx,
                ry = ry,
                ddc.object = ddc_object))
}
