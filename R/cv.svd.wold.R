
cv.svd.wold.check <- function( cv.svd ) {
    function( x, k=5, maxrank=min(n,p), tol=max(n,p)*1e-10, maxiter=100 ) {
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
    
        sets  <- choose.sets( n*p, k )
        
        press <- cv.svd( x, k, maxrank, tol, maxiter, sets )
        colnames( press ) <- 0:maxrank
        
        res   <- list( call=match.call(), k=k, maxrank=maxrank, 
                       press=press, sets=sets )
        class( res ) <- c("cvsvd_wold", "cvsvd")
        res
    }
}

cv.svd.wold.R.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    cv <- matrix( NA, k, maxrank+1 )
    
    for( j in seq_len( k ) ) {
        train  <- sets != j
        test   <- !train
        
        xtrain <- x
        xtrain[ test ] <- NA
        
        for( rank in seq( 0, maxrank ) ) {
            xhat <- suppressWarnings(
                        impute.svd( xtrain, rank, tol, maxiter )$x )
                        
            err  <- sum( ( xhat[ test ] - x[ test ] )^2 )
            cv[ j, rank+1 ] <- err
        }
    }
    
    cv
}
cv.svd.wold.R <- cv.svd.wold.check( cv.svd.wold.R.unchecked )

cv.svd.wold.C.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    presst <- .Call( "R_cv_svd_wold", x, k, maxrank, tol, maxiter, 
                     as.integer( sets-1 ) )    
    press  <- t( presst )
    press
}
cv.svd.wold.C <- cv.svd.wold.check( cv.svd.wold.C.unchecked )

cv.svd.wold <- cv.svd.wold.C
