
cv.svd <- function(x, K, L, max.rank = floor(min(nrow(x)/K, ncol(x)/L)))
{
    call <- match.call()
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    x <- as.matrix( x )
    M <- nrow(x)
    N <- ncol(x)
    
    if (M < 2)
        stop ("x should have at least two rows")
    if (N < 2)
        stop ("x should have at least two columns")
    if ((K > M) || (K <= 1)) 
        stop("K outside allowable range")
    if ((L > N) || (L <= 1))
        stop("L outside allowable range")
    if (max.rank < 0)
        stop("max.rank should be non-negative")
    
    K.o <- K; K <- round.fold(M, K);
    L.o <- L; L <- round.fold(N, L);
    
    if (K != K.o) 
        warning("K has been set to ", K)
    if (L != L.o)
        warning("L has been set to ", L)
    
    s.r <- choose.sets(M,K)
    s.c <- choose.sets(N,L)

    m.r <- min( table( s.r ) )
    n.c <- min( table( s.c ) )
    max.rank.o <- max.rank
    max.rank <- min(m.r, n.c, round( max.rank.o) )
    
    if (!missing(max.rank) && max.rank != max.rank.o)
        warning("max.rank has been set to ", max.rank)
    
    cv <- .Call("driver_svd", x, K, L, max.rank, s.r, s.c)
    list(call=call, K=K, L=L, max.rank=max.rank, cv=cv, 
         row.sets=s.r, col.sets=s.c, seed=seed)
}

round.fold <- function (n, K)
{
    K <- round(K)
    kvals <- unique(round(n/(1:floor(n/2))))
    temp <- abs(kvals - K)
    if (!any(temp == 0)) 
        K <- kvals[temp == min(temp)][1]
    K
}

choose.sets <- function (n, K)
{
    f <- ceiling(n/K)
    s <- sample(rep(1:K, f), n)
    n.s <- table(s)
    
    if ( length( n.s ) != K )
        choose.sets(n,K)
    else
        s
}
