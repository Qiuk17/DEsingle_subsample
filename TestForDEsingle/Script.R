options(warn = 1)
library(DEsingle)
#Generate zinb destribution
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
    record <- NULL
    TotalSample <- nrow(DEsingle.res.nDE)
    Sample <- DEsingle.res.nDE[sample(1:TotalSample, min(TotalSample, geneNum), replace = FALSE),]
    SampleNum <- nrow(Sample)
    for (i in 1: SampleNum) {
        res <- rbind(res, rzinbinom(Cell_in_group * 2, mean(Sample[i, "theta_1"], Sample[i, "theta_2"]), mean(Sample[i, "size_1"], Sample[i, "size_2"]), prob = mean(Sample[i, "prob_1"], Sample[i, "prob_2"])))
        record <- rbind(record, c(rep(mean(Sample[i, "theta_1"], Sample[i, "theta_2"]), 2), rep(mean(Sample[i, "size_1"], Sample[i, "size_2"]), 2), rep(mean(Sample[i, "prob_1"], Sample[i, "prob_2"]), 2), NA))
    }
    gc()
    #if sample is not enough, generate randomly.
    if (geneNum > SampleNum) {
        warning("Generating random data in nDE.")
        for (i in (SampleNum + 1): geneNum) {
            rand_theta <- runif(1)
            rand_size <- runif(1, min = 1, max = 1E5)
            rand_prob <- runif(1, min = 1E-4, max = 1)
            res <- rbind(res, rzinbinom(Cell_in_group * 2, rand_theta, rand_size, prob = rand_prob))
            record <- rbind(record, c(rand_theta, rand_theta, rand_size, rand_size, rand_prob, rand_prob, NA))
        }
    }
    rownames(res) <- paste0("gene_nDE_", 1:geneNum)
    rownames(record) <- paste0("gene_nDE_", 1:geneNum)
    write.csv(record, file = "./nDE_temp.csv")
    gc()
    return (res)
}

#Generate DEs genes
Generate_DEs <- function(Cell_in_group, geneNum = 500) {
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
    record <- NULL
    TotalSample <- nrow(DEsingle.res.DEs)
    Sample <- DEsingle.res.DEs[sample(1:TotalSample, min(TotalSample, geneNum), replace = FALSE),]
    SampleNum <- nrow(Sample)
    for (i in 1:SampleNum) {
        theta_1 <- Sample[i, "theta_1"]
        size_1 <- Sample[i, "size_1"]
        prob_1 <- Sample[i, "prob_1"]
        theta_2 <- Sample[i, "theta_2"]
        size_2 <- Sample[i, "size_2"]
        prob_2 <- Sample[i, "prob_2"]
        res <- rbind(res, c(rzinbinom(Cell_in_group, theta_1, size_1, prob = prob_1), rzinbinom(Cell_in_group, theta_2, size_2, prob = prob_2)))
        record <- rbind(record, c(theta_1, theta_2, size_1, size_2, prob_1, prob_2, "DEs"))
    }
    gc()
    #if sample is not enough, generate randomly.
    delta <- 0.3
    if (geneNum > SampleNum) {
        warning("Generating random data in DEs.")
        for (i in (SampleNum + 1):geneNum) {
            rand_theta_1 <- runif(1)
            rand_theta_2 <- runif(1)
            while (abs(rand_theta_1 - rand_theta_2) < delta) {
                rand_theta_1 <- runif(1)
                rand_theta_2 <- runif(1)
            }
            rand_size <- runif(1, min = 1, max = 1E5)
            rand_prob <- runif(1, min = 1E-4, max = 1)
            res <- rbind(res, c(rzinbinom(Cell_in_group, rand_theta_1, rand_size, prob = rand_prob), rzinbinom(Cell_in_group, rand_theta_2, rand_size, prob = rand_prob)))
            record <- rbind(record, c(rand_theta_1, rand_theta_2, rand_size, rand_size, rand_prob, rand_prob, "DEs"))
        }
    }
    rownames(res) <- paste0("gene_DEs_", 1:geneNum)
    rownames(record) <- paste0("gene_DEs_", 1:geneNum)
    write.csv(record, file = "./DEs_temp.csv")
    gc()
    return(res)
}

#Generate DEg genes
Generate_DEg <- function(Cell_in_group, geneNum = 500) {
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
    record <- NULL
    TotalSample <- nrow(DEsingle.res.DEg)
    Sample <- DEsingle.res.DEg[sample(1:TotalSample, min(TotalSample, geneNum), replace = FALSE),]
    SampleNum <- nrow(Sample)
    for (i in 1:SampleNum) {
        theta_1 <- Sample[i, "theta_1"]
        size_1 <- Sample[i, "size_1"]
        prob_1 <- Sample[i, "prob_1"]
        theta_2 <- Sample[i, "theta_2"]
        size_2 <- Sample[i, "size_2"]
        prob_2 <- Sample[i, "prob_2"]
        res <- rbind(res, c(rzinbinom(Cell_in_group, theta_1, size_1, prob = prob_1), rzinbinom(Cell_in_group, theta_2, size_2, prob = prob_2)))
        record <- rbind(record, c(theta_1, theta_2, size_1, size_2, prob_1, prob_2, "DEg"))
    }
    gc()
    #if sample is not enough, generate randomly.
    delta <- 0.3
    delta_size <- 1000
    if (geneNum > SampleNum) {
        warning("Generating random data in DEg.")
        for (i in (SampleNum + 1):geneNum) {
            rand_theta_1 <- runif(1)
            rand_theta_2 <- runif(1)
            while (abs(rand_theta_1 - rand_theta_2) < delta) {
                rand_theta_1 <- runif(1)
                rand_theta_2 <- runif(1)
            }
            rand_size_1 <- runif(1, min = 1, max = 1E5)
            rand_size_2 <- runif(1, min = 1, max = 1E5)
            while (abs(rand_size_1 - rand_size_2) < delta_size) {
                rand_size_1 <- runif(1, min = 1, max = 1E5)
                rand_size_2 <- runif(1, min = 1, max = 1E5)
            }
            rand_prob_1 <- runif(1, min = 1E-4, max = 1)
            rand_prob_2 <- runif(1, min = 1E-4, max = 1)
            while (abs(rand_prob_1 - rand_prob_2) < delta) {
                rand_prob_1 <- runif(1, min = 1E-4, max = 1)
                rand_prob_2 <- runif(1, min = 1E-4, max = 1)
            }
            res <- rbind(res, c(rzinbinom(Cell_in_group, rand_theta_1, rand_size_1, prob = rand_prob_1), rzinbinom(Cell_in_group, rand_theta_2, rand_size_2, prob = rand_prob_2)))
            record <- rbind(record, c(rand_theta_1, rand_theta_2, rand_size_1, rand_size_2, rand_prob_1, rand_prob_2, "DEg"))
        }
    }
    rownames(res) <- paste0("gene_DEg_", 1:geneNum)
    rownames(record) <- paste0("gene_DEg_", 1:geneNum)
    write.csv(record, file = "./DEg_temp.csv")
    gc()
    return(res)
}

#Generate DEa genes
Generate_DEa <- function(Cell_in_group, geneNum = 500) {
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
    record <- NULL
    TotalSample <- nrow(DEsingle.res.DEa)
    Sample <- DEsingle.res.DEa[sample(1:TotalSample, min(TotalSample, geneNum), replace = FALSE),]
    SampleNum <- nrow(Sample)
    for (i in 1:SampleNum) {
        theta_1 <- Sample[i, "theta_1"]
        size_1 <- Sample[i, "size_1"]
        prob_1 <- Sample[i, "prob_1"]
        theta_2 <- Sample[i, "theta_2"]
        size_2 <- Sample[i, "size_2"]
        prob_2 <- Sample[i, "prob_2"]
        res <- rbind(res, c(rzinbinom(Cell_in_group, theta_1, size_1, prob = prob_1), rzinbinom(Cell_in_group, theta_2, size_2, prob = prob_2)))
        record <- rbind(record, c(theta_1, theta_2, size_1, size_2, prob_1, prob_2, "DEa"))
    }
    gc()
    #if sample is not enough, generate randomly.
    delta <- 0.2
    delta_size <- 1E3
    if (geneNum > SampleNum) {
        warning("Generating random data in DEa.")
        for (i in (SampleNum + 1):geneNum) {
            rand_theta <- runif(1)
            rand_size_1 <- runif(1, min = 1, max = 1E5)
            rand_size_2 <- runif(1, min = 1, max = 1E5)
            while (abs(rand_size_1 - rand_size_2) < delta_size) {
                rand_size_1 <- runif(1, min = 1, max = 1E5)
                rand_size_2 <- runif(1, min = 1, max = 1E5)
            }
            rand_prob_1 <- runif(1, min = 0.5, max = 1)
            rand_prob_2 <- runif(1, min = 0.5, max = 1)
            while (abs(rand_prob_1 - rand_prob_2) < delta) {
                rand_prob_1 <- runif(1, min = 0.5, max = 1)
                rand_prob_2 <- runif(1, min = 0.5, max = 1)
            }
            res <- rbind(res, c(rzinbinom(Cell_in_group, rand_theta, rand_size_1, prob = rand_prob_1), rzinbinom(Cell_in_group, rand_theta, rand_size_2, prob = rand_prob_2)))
            record <- rbind(record, c(rand_theta, rand_theta, rand_size_1, rand_size_2, rand_prob_1, rand_prob_2, "DEa"))
        }
    }
    rownames(res) <- paste0("gene_DEa_", 1:geneNum)
    rownames(record) <- paste0("gene_DEa_", 1:geneNum)
    write.csv(record, file = "./DEa_temp.csv")
    gc()
    return(res)
}



#start generation
#set these nums
Cell_in_one_group <- 500
nDE_num <- 500
DEg_num <- 500
DEs_num <- 500
DEa_num <- 500
gene_sum <- sum(nDE_num, DEg_num, DEa_num, DEs_num)
counts <- Generate_DEs(Cell_in_one_group, DEs_num)
counts <- rbind(counts, Generate_DEg(Cell_in_one_group, DEg_num))
counts <- rbind(counts, Generate_DEa(Cell_in_one_group, DEa_num))
counts <- rbind(counts, Generate_nDE(Cell_in_one_group, nDE_num))

colnames(counts) <- paste0("cell_", 1:(Cell_in_one_group * 2))
record <- rbind(read.csv("./nDE_temp.csv"), read.csv("./DEs_temp.csv"), read.csv("./DEg_temp.csv"), read.csv("./DEa_temp.csv"))
record <- cbind(record, rep(NA, gene_sum))
colnames(record) <- c("gene_name", "theta_1", "theta_2", "size_1", "size_2", "prob_1", "prob_2", "Type", "DEsingle_Type")
groups <- as.factor(c(rep(TRUE, Cell_in_one_group), rep(FALSE, Cell_in_one_group)))

res <- DEtype(DEsingle(counts, groups), 0.05)
for (i in 1:gene_sum) {
    try(record[i, "DEsingle_Type"] <- res[which(rownames(res) == as.character(record[i, "gene_name"])), "Type"])
}

write.csv(res, file = "./DEsingle_result.csv")
write.csv(record, file = "./origin_data.csv")
message("Script done.")
