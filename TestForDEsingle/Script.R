#生成ZINB模型
rzinbinom <- function(n, theta, size, prob = NULL, mu = NULL) {
    #n
    if (length(n) != 1)
        stop("n must be a single number.")
    if (!is.numeric(n))
        stop("Invalid class of n.")
    if (n <= 0)
        stop("The length must be above 0.")
    #theta
    if (length(theta) != 1)
        stop("theta must be a single number.")
    if (!is.numeric(theta))
        stop("invalid class of theta.")
    if (theta > 1.0 || theta < 0.0)
        stop("theta is invalid.")
    #size
    if (length(size) != 1)
        stop("size must be a single number.")
    if (!is.numeric(size))
        stop("Invalid class of size.")
    if (size <= 0.0)
        stop("size must be positive.")
    #prob and mu
    if (!(is.null(mu) || is.null(prob)))
        stop("mu and prob cannot be set at the same time.")
    if (is.null(mu) && is.null(prob))
        stop("mu or prob must be set.")
    if (length(mu) != 1 && length(prob) != 1)
        stop("mu or prob must be a single number.")
    if (!(is.numeric(prob) || is.numeric(mu)))
        stop("Invalid class of prob or mu.")
    #generate NB
    res <- NA
    if (is.null(mu))
        res <- rnbinom(n = n, size = size, prob = prob)
    if (is.null(prob))
        res <- rnbinom(n = n, size = size, mu = mu)
    if (any(is.na(res)))
        stop("NA detected.")
    #randomly set zeros
    res[sample(1:n, ceiling(theta * n), replace = FALSE)] <- 0
    return (res)
    }

#Get four types of the genes
DEsingle.res <- read.csv("./small_result.csv")
DEsingle.res.DEg <- DEsingle.res[which(DEsingle.res$Type == "DEg"),]
DEsingle.res.DEa <- DEsingle.res[which(DEsingle.res$Type == "DEa"),]
DEsingle.res.DEs <- DEsingle.res[which(DEsingle.res$Type == "DEs"),]
DEsingle.res.nDE <- DEsingle.res[which(is.na(DEsingle.res$Type)),]

set.seed(Sys.time())

#Generate genes with no different expression
Generate_nDE <- function(Cell_in_group, geneNum = 500) {
    if (length(Cell_in_group) != 1)
        stop("length of Cell_in_group is not 1.")
    if (!is.numeric(Cell_in_group))
        stop("Cell_in_group must be a number.")
    if (length(geneNum) != 1)
        stop("length of geneNum is not 1.")
    if (!is.numeric(geneNum))
        stop("geneNum must be a number.")
    if (geneNum <= 0)
        stop("geneNum must be positive.")

    res <- NULL
    TotalSample <- nrow(DEsingle.res.nDE)
    Sample <- DEsingle.res.nDE[sample(1:TotalSample, min(TotalSample, geneNum), replace = FALSE),]
    SampleNum <- nrow(Sample)
    for (i in 1:SampleNum) {
        res <- rbind(res, rzinbinom(Cell_in_group * 2, mean(Sample[i, "theta_1"], Sample[i, "theta_2"]), mean(Sample[i, "size_1"], Sample[i, "size_2"]), prob = mean(Sample[i, "prob_1"], Sample[i, "prob_2"])))
    }

    gc()
    return (res)
}

#start generation
Cell_in_Groups <- c(100, 200, 500, 800, 1000, 1500, 2000)

for (Sample_num in Cell_in_Groups) {

}
