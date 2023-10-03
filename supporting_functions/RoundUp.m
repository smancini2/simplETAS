
%% CODE TO BE TRANSLATED

RoundUp <- function(x, digits=0) {
    factor <- 10^digits
    as.numeric(trunc(x*factor + 0.5)/factor)
}
