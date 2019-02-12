#Writing functions

#nameoffunction <- function(arg1, arg2) {statements
#return(something)}

df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
# df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
# df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
# df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
# df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d)) 

#intent: df$a <- rescale(df$a)
x <- function(x, y){
  xmin <- min(x)
  x <- (x - xmin)/(max(x)-xmin)
}



rescale <- function(x, na.rm=TRUE, plot=FALSE) {
  rng <-range(x, na.rm=na.rm)
  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see it in ...")
  return(answer)
}