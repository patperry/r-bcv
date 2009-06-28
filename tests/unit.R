
get.x11 <- function( cv, x, i, j ) x[ cv$row != i, cv$col != j ]
get.x12 <- function( cv, x, i, j ) x[ cv$row != i, cv$col == j ]
get.x21 <- function( cv, x, i, j ) x[ cv$row == i, cv$col != j ]
get.x22 <- function( cv, x, i, j ) x[ cv$row == i, cv$col == j ]

check <- function( cv, x )
{
    expected <- matrix( 0, cv$max.rank + 1, cv$K * cv$L )
    for( i in 1:cv$K )
    {
        for( j in 1:cv$L )
        {
            x11 <- get.x11( cv, x, i, j );
            x11.svd <- svd( x11, nu=cv$max.rank, nv=cv$max.rank )
            u1 <- x11.svd$u
            v1 <- x11.svd$v
            d  <- x11.svd$d
            
            x12 <- get.x12( cv, x, i, j )
            x21 <- get.x21( cv, x, i, j )
            x22 <- get.x22( cv, x, i, j )
            
            col <- i + ( j-1 )*cv$K 
            expected[ 1, col ] <- mean( x22^2 )
            for( k in seq_len( cv$max.rank ) )
            {
                x22.hat <- (     ( x21 %*% v1[,1:k,drop=FALSE] )
                             %*% diag( 1.0 / d[1:k], k, k )
                             %*% ( t( u1[,1:k,drop=FALSE] ) %*% x12 ) )
                expected[ k + 1, col ] <- mean( ( x22 - x22.hat )^2 )                
            }
        }
    }
    
    actual <- cv$cv
    
    if( !all( abs( ( expected - actual )/expected < 1e-8 ) ) )
    {
        cat( "Expected:\n")
        print( expected )
        
        cat( "Actual:\n")
        print( actual )
        TRUE
    }
    else
        FALSE
}

rmatrix <- function( size )
{
    M <- floor( runif( 1, 4, max( 5, size ) ) )
    N <- floor( runif( 1, 4, max( 5, size ) ) )
    x <- matrix( rnorm( M*N ), M, N )
    x
}

rfold <- function( N )
{
    floor( runif( 1, 2, max( 3, N ) ) )
}

rmaxrank <- function( M, N, K, L )
{
    max.rank <- min( M*( 1 - 1/L ) , N*( 1 - 1/K ) )
    floor( runif( 1, 0, max.rank + 1 ) )
}

rtest <- function( size )
{
    x        <- rmatrix( size )
    K        <- rfold( nrow( x ) )
    L        <- rfold( ncol( x ) )
    max.rank <- rmaxrank(  nrow( x ), ncol( x ), K, L )

    list( x=x, K=K, L=L, max.rank=max.rank )
}

sizes <- function( ntests, each=2 )
{
    if( ntests > 0 )
        rep(4 + sqrt( seq( 0, ntests/each, length=ceiling( ntests/each ) ) )
                    , each=each )[ 1:ntests ]
    else
        c()
}

require( bcv )
ntests <- 1000
s <- sizes( ntests )
set.seed( 0 )

nsuc <- 0
for (size in s)
{
    cat( '.' )
    test <- rtest( size )
    
    cv <- suppressWarnings( cv.svd( test$x, test$K, test$L, test$max.rank ) )
    if( !check( cv, test$x ) )
        nsuc <- nsuc + 1
}

if( nsuc == ntests ) {
    cat("Passed", nsuc, "tests!\n")
} else {
    cat("Not all tests passed.\n")
}

