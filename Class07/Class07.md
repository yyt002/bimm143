Class 7 functions and packages
================
Yvette Tan
1/29/2019

Functions revisit
-----------------

``` r
source("http://tinyurl.com/rescale-R")
```

Lets try the rescale() function out

``` r
rescale( c(1, 5, 10))
```

    ## [1] 0.0000000 0.4444444 1.0000000

Let's try **rescale2()** with the **stop()** function catch for non-numeric input

**which()** : which element is true using

``` r
x <- c(3, 7, NA, 4, 8, NA)
which(is.na(x))
```

    ## [1] 3 6

``` r
#Let's define an example x and y
x <- c(1,2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
sum(is.na(x))
```

    ## [1] 2

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
#putting together!
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)
y3 <- c(1, NA, NA, NA, NA)
y4 <- c(1, 1, NA)
both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na(x, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
both_na(y1, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na(y4, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 1

``` r
#when length isn't equal, the shorter will rerun the vector from left to right until it matches the length of the longer vector
```

``` r
both_na3(x, y1)
```

    ## Found 2 NA's at position(s):2, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 2 3
