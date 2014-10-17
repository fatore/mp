library(svd)

lamp = function(X, sample.indices = NULL, Ys = NULL) {
    if (is.null(sample.indices)) {
        n = nrow(X)
        sample.indices = sample(1:n, sqrt(n))
    }

    Xs = X[sample.indices, ]

    if (is.null(Ys)) {
        Ys = forceScheme(as.matrix(dist(Xs)))
    }

    Y = t(apply(X, 1, function(point) {
        alphas = apply(Xs, 1, function(sample.point) sum((sample.point - point)^2))
        alphas.sum = sum(alphas)
        alphas.sqrt = sqrt(alphas)

        X.til = colSums(Xs * alphas) / alphas.sum
        Y.til = colSums(Ys * alphas) / alphas.sum
        X.hat = sweep(Xs, 2, X.til, '-')
        Y.hat = sweep(Ys, 2, Y.til, '-')

        A = X.hat * alphas.sqrt
        B = Y.hat * alphas.sqrt
        AtB = t(A) %*% B
        s = propack.svd(AtB, neig = 2) # We just need the first 2
        M = s$u %*% s$v

        y = (point - X.til) %*% M + Y.til
        y
    }))
}
