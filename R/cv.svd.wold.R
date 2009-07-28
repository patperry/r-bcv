
cv.svd.wold.check <- function( cv.svd ) {
    function( x, k, maxrank=min(n,p), tol=max(n,p)*1e-10, maxiter=100 ) {
        x <- as.matrix( x )
        n <- nrow( x )
        p <- ncol( x )
    
        if (is.complex( x ))
            stop ("x cannot be complex")
        if (n*p < 2)
            stop ("x should have at least two elements")
        if (k > n*p || k <= 1)
            stop("k outside allowable range")
        if (maxrank < 0 || maxrank > min(n,p))
            stop("maxrank outside allowable range")
    
        storage.mode( x ) <- "double"
    
        k.o <- k; k <- round.fold( n*p, k );
        if (k != k.o) 
            warning("k has been set to ", k)
    
        sets <- choose.sets( n*p, k )
        cv   <- cv.svd ( x, k, maxrank, tol, maxiter, sets )
        list( k=k, maxrank=maxrank, cv=cv, sets=sets )
    }
}

cv.svd.wold.R.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    cv <- matrix( NA, maxrank+1, k )
    
    for( j in seq_len( k ) ) {
        train  <- sets != j
        test   <- !train
        
        xtrain <- x
        xtrain[ test ] <- NA
        
        for( rank in seq( 0, maxrank ) ) {
            xhat <- suppressWarnings(
                        svd.impute( xtrain, rank, tol, maxiter )$x )
                        
            err  <- sum( ( xhat[ test ] - x[ test ] )^2 )
            cv[ rank+1, j ] <- err
        }
    }
    
    cv
}
cv.svd.wold.R <- cv.svd.wold.check( cv.svd.wold.R.unchecked )

cv.svd.wold.C.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    n    <- nrow( x )
    p    <- ncol( x )
    sets <- choose.sets( n*p, k )
    
    res <- .Call( "R_bcv_svd_wold", x, k, maxrank, tol, maxiter, 
                  as.integer( sets-1 ) )    
    res
}
cv.svd.wold.C <- cv.svd.wold.check( cv.svd.wold.C.unchecked )

cv.svd.wold <- cv.svd.wold.C
