
print.cvsvd <- function( cv, digits=max(3, getOption("digits") - 3), ... ) {
    cat( "\nCall:\n", deparse( cv$call ), "\n\n", sep = "" )
    
    press   <- cv$press
    maxrank <- cv$maxrank
    
    K           <- nrow( press )
    rank        <- seq( from=0, to=maxrank, by=1 )
    press.mean  <- apply( press, 1, mean )
    press.se    <- apply( press, 1, sd ) / sqrt( K )
    min.rank    <- which.min( press.mean ) - 1
    min.rank.se <- min( which( press.mean 
                               <= press.mean[ min.rank+1 ] 
                                  + press.mean[ min.rank+1 ] ) ) - 1
    
    rank.fmt   <- format( rank )
    rank.width <- max( nchar( rank.fmt[ 1 ] ), nchar( "Rank" ) + 1 )
    
    mean.fmt   <- format( press.mean, digits=digits )
    mean.width <- max( nchar( mean.fmt[ 1 ] ), nchar( "PRESS" ) + 1 )
    
    se.fmt     <- format( press.se, digits=digits )
    se.width   <- max( nchar( se.fmt[ 1 ] ), nchar( "SE" ) + 1 )
    
    fmt <- paste( " %", rank.width, "s", 
                  "  %", mean.width, "s",
                  "  %", se.width, "s", sep="" )
    cat( sprintf( fmt, "Rank", "PRESS", "SE" ), "\n", sep="" )
    cat( " " )
    cat( do.call( paste, 
         c( as.list( rep("-", 
                rank.width + 2 + mean.width + 2 + se.width) ), 
            sep='' ) ) )
    cat( "\n" )
    
    for( i in seq_len( maxrank+1 ) ) {
        rank <- i - 1
        mean <- press.mean[ i ]
        se   <- press.se[ i ]
        cat( sprintf( fmt, rank.fmt[ i ], mean.fmt[ i ], se.fmt[ i ] ) )
        if ( rank == min.rank && rank == min.rank.se ) {
            cat( " * " )
        } else if ( rank == min.rank ) {
            cat( " * " )
        } else if ( rank == min.rank.se ) {
            cat( " + " )
        }
        cat( "\n" )
    }
    
    cat( "\n" )
    
    invisible( cv )
}
