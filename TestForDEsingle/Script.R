# 产生符合ZINB(θ, μ, size, prob)的一组数字，长度为n 
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

