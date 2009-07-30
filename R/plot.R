
plot.cvsvd <- function( x, errorbars=TRUE,
                        col="blue", col.errorbars="gray50", 
                        add=FALSE, xlab="Rank", ylab="Prediction Error",
                         ... ) {
    press   <- x$press
    maxrank <- x$maxrank
    
    K           <- nrow( press )
    rank        <- seq( from=0, to=maxrank, by=1 )
    press.mean  <- apply( press, 2, mean )
    press.se    <- apply( press, 2, sd ) / sqrt( K )
    
    if( !add ) {
        if( errorbars ) {
            plot( c(rank-0.2,rank+0.2), press.mean+c(-press.se, press.se), 
                  t='n', xlab=xlab, ylab=ylab, ... )
        } else {
            plot( rank, press.mean, t='n', xlab=xlab, ylab=ylab, ... )
        }
    }

    lines( rank, press.mean, col=col )
    
    if( errorbars ) {
        segments( rank-0.2, press.mean-press.se, 
                  rank+0.2, press.mean-press.se,
                  col=col.errorbars )

        segments( rank, press.mean-press.se, 
                  rank, press.mean+press.se,
                  col=col.errorbars )

        segments( rank-0.2, press.mean+press.se, 
                  rank+0.2, press.mean+press.se,
                  col=col.errorbars )
    }

    points( rank, press.mean, col=col, pch=16, cex=0.6 )
    
    invisible( list( k=rank, press=press.mean ) )
}
