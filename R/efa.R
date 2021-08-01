#' @import mice
#' @import psych
#' @import GPArotation
#' @import Cairo
#' @import data.table
#' @importFrom grDevices dev.off
#' @importFrom graphics abline grid hist par
#' @importFrom stats cor cov lm mahalanobis qchisq qqnorm rchisq rstudent symnum

# for multiple imputation library
library("mice")
# for running the efa analysis
library("psych")
library("GPArotation")

#' @export
filter_multi_df_outliers_missing = function(
        df,
        value_cols=NULL) {
    # first we eliminate the items with too many missing values
    if (is.null(value_cols)) {
        value_cols = colnames(df)
    }
    value_df = df[ , value_cols, drop=FALSE]

    percent_na = function (x) {
        sum(is.na(x)) / length(x) * 100}

    missing = apply(value_df, 1, percent_na)
    table(missing)

    missing_df = subset(df, missing > 5)
    df = subset(df, missing <= 5)

    df0 = mice::mice(df)
    df = complete(df0, 1)

    # outliers
    cutoff = stats::qchisq(1-.001, length(value_cols))
    mahal = stats::mahalanobis(
        value_df,
        colMeans(value_df),
        cov(value_df))
    # outliers
    outlier_df = subset(df, mahal >= cutoff)

    # exclude outliers
    df = subset(df, mahal < cutoff)

    return(list(
        "df"=df,
        "mahal_cutoff"=cutoff,
        "missing_df"=missing_df,
        "outlier_df"=outlier_df
    ))
}

#' @export
check_efa_assumptions = function(
        df,
        value_cols=NULL,
        figs_path=NULL) {
    # check efa assumptions
    if (is.null(value_cols)) {
        value_cols = colnames(df)
    }
    value_df = df[ , value_cols, drop=FALSE]
    cor_mat = stats::cor(value_df, use="pairwise.complete.obs")
    stats::symnum(cor_mat)

    # assumption set up
    random = stats::rchisq(nrow(value_df), 7)
    dummy_model = stats::lm(random~., data=value_df)
    standardized = stats::rstudent(dummy_model)
    fitted = base::scale(dummy_model[["fitted.values"]])

    plt_func = function() {
        par(mfrow=c(2, 2))
        # normality
        hist(standardized, main="residuals histogram")

        # linearity
        qqnorm(standardized, main="residuals qq plot")
        abline(0, 1)
        grid()

        # homogeneity
        plot(fitted, standardized, main="homogeneity")
        abline(0, 0)
        abline(v=0)
    }

    plt_file_name = paste0(figs_path, "dummy_model_check.png")

    if (!is.null(figs_path)) {
        Cairo::Cairo(file=plt_file_name, type="png")
        plt_func()
        dev.off()
    }

    # correlation adequacy Bartlett's test
    barlet_test = psych::cortest.bartlett(cor_mat, n=nrow(value_df))
    barlet_test

    # sampling adequacy KMO test
    # we want this to be larger than 0.7
    kaiser_sampling_adequacy = psych::KMO(cor_mat)
    kaiser_sampling_adequacy

    return(list(
        "cor_mat"=cor_mat,
        "barlet_test"=barlet_test,
        "kaiser_sampling_adequacy"=kaiser_sampling_adequacy,
        "plt_func"=plt_func
    ))

}

#' Determines how many factors are suitable
#' @export
get_factor_num = function(
        df,
        value_cols=NULL,
        figs_path=NULL) {
    # get factor number
    if (is.null(value_cols)) {
        value_cols = colnames(df)
    }
    value_df = df[ , value_cols, drop=FALSE]

    fa_parallel_model = psych::fa.parallel(value_df, fm="ml", fa="fa")
    kaiser_factor_num_old = sum(fa_parallel_model[["fa.values"]] >= 1.0) # old kaiser criterion
    kaiser_factor_num_new = sum(fa_parallel_model[["fa.values"]] >= .7) # new kaiser criterion

    plt_func = function() {
        fa_parallel_model = psych::fa.parallel(value_df, fm="ml", fa="fa")
    }

    plt_file_name = paste0(figs_path, "scree_plot.png")
    if (!is.null(figs_path)) {
        plt_func()
    }

    # Number of factors with eigen values > eigen values of random data
    factor_num = fa_parallel_model[["nfact"]]

    return(list(
        "fa_parallel_model"=fa_parallel_model,
        "kaiser_factor_num_old"=kaiser_factor_num_old,
        "kaiser_factor_num_new"=kaiser_factor_num_new,
        "factor_num"=factor_num,
        "plt_func"=plt_func))
}


#' @export
fit_fa_model = function(
        df,
        factor_num,
        value_cols=NULL,
        rotate="oblimin",
        fm="ml",
        loading_thresh=0.3) {
    # fit fa model and return loading_df rejected_items
    if (is.null(value_cols)) {
        value_cols = colnames(df)
    }
    value_df = df[ , value_cols, drop=FALSE]
    other_cols = setdiff(colnames(df), value_cols)

    fa_model = psych::fa(
        value_df,
        nfactors=factor_num,
        rotate=rotate,
        fm=fm)

    tuker_lewis_index = fa_model[["TLI"]]
    rmse = fa_model[["rms"]]
    rmsea = fa_model[["RMSEA"]]
    loadings = fa_model[["loadings"]]

    loadings_mat = round(data.frame((unclass(loadings))), 4)
    item = rownames(loadings_mat)
    loadings_df = cbind(item, as.data.frame(loadings_mat))
    loading_thresh_mat = (abs(loadings_mat) > loading_thresh)
    loadings_df[["associated_factor_num"]] = rowSums(loading_thresh_mat)

    # to each item assign the factor with which it is most associated
    loadings_df[["factor"]] = (
        apply(
            X=loadings_mat,
            FUN=function(x)which(abs(x) == max(abs(x)))[1],
            MARGIN=1))

    loadings_dt = data.table::data.table(loadings_df)

    # we reject any items which is not associated with exactly one factor
    rejected_items = loadings_dt[associated_factor_num != 1][["item"]]
    rejected_items = as.character(rejected_items)
    accepted_items = loadings_dt[associated_factor_num == 1][["item"]]
    accepted_items = as.character(accepted_items)
    is_accepted_bool = loadings_df[["item"]] %in% accepted_items
    loadings_df_filtered = loadings_df[is_accepted_bool, ]
    dim(loadings_df_filtered)
    factor_item_num = table(loadings_df_filtered[["factor"]])

    factors = unique(loadings_df_filtered[["factor"]])

    factor_item_list = list()
    for (f in factors) {
        f_ind = loadings_df_filtered[["factor"]] == f
        factor_item_list[[f]] = loadings_df_filtered[f_ind, "item"]
    }

    df_filtered = df[ , c(other_cols, accepted_items), drop=FALSE]
    min_item_num_per_factor = min(factor_item_num)

    stat = fa_model[["STATISTIC"]]
    dof = fa_model[["dof"]]
    null_chisq = fa_model[["null.chisq"]]
    null_dof = fa_model[["null.dof"]]
    cfi = 1 - (stat - dof) / (null_chisq - null_dof)

    return(list(
        "fa_model"=fa_model,
        "df_filtered"=df_filtered,
        "loadings"=loadings,
        "loadings_mat"=loadings_mat,
        "loadings_df_filtered"=loadings_df_filtered,
        "factor_item_list"=factor_item_list,
        "factor_item_num"=factor_item_num,
        "min_item_num_per_factor"=min_item_num_per_factor,
        "tuker_lewis_index"=tuker_lewis_index,
        "rmse"=rmse,
        "rmsea"=rmsea,
        "rejected_items"=rejected_items,
        "accepted_items"=accepted_items,
        "cfi"=cfi))
}
