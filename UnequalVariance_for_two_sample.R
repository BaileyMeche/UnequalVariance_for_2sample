# -*- coding: utf-8 -*-
"""
# Preliminaries
"""

install.packages('openxlsx')
library(openxlsx)

library ( MASS )
library(tidyverse)


install.packages('Partiallyoverlapping')
library(Partiallyoverlapping)

install.packages('lme4')
library(lme4)

"""# ---- Generate data

"""

data <- function (n , n1 , n2 , mu1 , mu2 , sigma_x,sigma_y, cor , diff ,
              dist =c ( " Normal " ,
              " Symmetric ␣ heavy ␣ tailed " ,
              " Skewed ␣ heavy ␣ tailed " ,
              " Skewed ␣ light ␣ tailed " ))
{
  mu <- c ( mu1 , mu2 )
  # Correlation matrix
  normCor <- matrix ( c (sigma_x^2 , cor*sigma_x*sigma_y , cor*sigma_x*sigma_y , sigma_y^2) , nrow = 2)

  # Generate data from Normal distribution
  data <- mvrnorm ( n = n+ n1 + n2 , mu = mu , Sigma = normCor )
  if ( dist == " Normal " ){
    data = data
    } else ( print ( " Specify ␣ distribution " ))

    # Delete missing data from variables
    if (( n1 != 0) && ( n2 != 0)){
    data [1: n2 ,1] <- NA
    data [( n + n2 +1):( n + n1 + n2 ) ,2] <- NA

    # Rename data to paired and unpaired
    xp <- data [( n2 +1):( n + n2 ) ,1]
    yp <- data [( n2 +1):( n + n2 ) ,2]
    xu <- data [( n + n2 +1):( n + n1 + n2 ) ,1]
    yu <- data [1: n2 ,2]
    }
    else {
      xp <- data [ ,1]
      yp <- data [ ,2]
      xu <- NA
      yu <- NA
    }
    # Sample correlation
  r <- cor ( xp , yp ) #Pearson correlation
  return ( list ("xp" = xp , "yp" = yp , "xu" =xu , "yu"= yu , "r"=r ,
  "n" =n , "n1" =n1 , "n2" = n2 ))
}

paired <- function ( x1 , x2 )
  { n = length ( x1 )
  y1 = numeric ( n )
  y2 = numeric ( n )
  for ( i in 1: n)
  { u = runif (1)
  if ((u -0.5) <=0) {
    y1 [i ]= x1 [ i]
    y2 [i ]= x2 [ i]
  } else {
    if (( u -0.5) >0) {
      y1 [i ]= x2 [ i]
      y2 [i ]= x1 [ i]
  }
  }
  }
  return ( list ( "y1" = y1 ,"y2" = y2 ))
  }
  unpaired <- function ( x1 , x2 )
  {
    m1 = length ( x1 )
    m2 = length ( x2 )
    c = c (x1 , x2 )
    n = length ( c )
    y = sample (c , replace = F )
    y1 =y [1: m1 ]
    y2 =y [( m1 +1): n ]
  return ( list ("y1" = y1 , "y2" = y2 ))
  }

"""# Statistics

## Bhoj
"""

Bhoj_Zb <- function(xp,xu,yp,yu,n,n1,n2){
 xbar1 <- mean ( xp )
    xbar2 <- mean ( xu )
    ybar1 <- mean ( yp )
    ybar2 <- mean ( yu )

    a11 <- sum (( xp - xbar1 )^2)
    a22 <- sum (( yp - ybar1 )^2)
    a12 <- sum (( xp - xbar1 ) * ( yp - ybar1 ))
    b1 <- sum (( xu - xbar2 )^2)
    b2 <- sum (( yu - ybar2 )^2)

    u <- (2 * a12 ) / ( a11 + a22 )
    w <- ( n1 * ( n +((1+ u )* n2 ))) / (( n * ( n1 + n2 ))+(2 * (1+ u ) * n1 * n2 ))

    s <- (1+ u ) /2
    f1 <- n -1
    f2 <- n1 + n2 -2
    f3 <- n + n1 + n2 -3
    d2 <- 2 * xbar2 - xbar1 - ybar1
    d3 <- xbar1 + ybar1 -2 * ybar2
    d <- w * d2 + (1 - w ) * d3

    t1 <- (( xbar1 - ybar1 ) * sqrt ( n )) /sqrt (( a11 + a22 -2 * a12 ) / (n -1))
    t3 <- d / sqrt (((4 *s * ( b1 + b2 )+ a11 + a22 +2 * a12 ) /( n + n1 + n2 -3)) * (( w ^2 / (s * n1) )+(((1 - w )^2) /(s * n2) )+(((1 -2 * w )^2) /n ))) #corrected
    F1 <- 1+((2 * t1 ^2) / f1 )+(2 * t1 / sqrt ( f1 )) *sqrt (1+(( t1 ^2) / f1 ))
    F3 <- 1+((2 * t3 ^2) / f3 )+(2 * t3 / sqrt ( f3 )) *sqrt (1+(( t3 ^2) / f3 ))

    U1 <- ((1 -(2 / (9 * f1 ))) * ( F1 ^(1 / 3) -1)) /sqrt ((2 / (9 * f1 )) * ( F1 ^(2 / 3)+1))
    U3 <- ((1 -(2 / (9 * f3 ))) * ( F3 ^(1 / 3) -1)) /sqrt ((2 / (9 * f3 )) * ( F3 ^(2 / 3)+1))

    lb <- (1+ sqrt (( n1 * n2 * (1 - u )) /(2 *n * n2 * ( w ^2)+ 2 * n* n1 * ((1 - w )^2)+ n1 * n2 * ((1 -2 * w )^2) * (1+ u ))))^(-1)
    Zb <- ( lb * U1 +(1 - lb ) * U3 ) / sqrt (( lb ^2)+(1 - lb )^2)

    #pvalue
    if (( n1 == 0) && ( n2 == 0)) { pvalue.zb <- t.test ( xp ,yp ,alternative = "greater" , paired = TRUE ,   var.equal = TRUE )$p.value
    }
    else {
      pvalue.zb <- pnorm ( Zb, lower.tail = FALSE ) # ONE SIDED
    }

  return(c(Zb, pvalue.zb))
}

Bhoj_T <- function(xp,xu,yp,yu,n,n1,n2){
  m1 <- n/(n+n1)
  m2  <- n/(n+n2)
  xbar1 <- mean ( xp )
  xbar2 <- mean ( xu )
  ybar1 <- mean ( yp )
  ybar2 <- mean ( yu )
  xbar <- sum (c(xp,xu))/(n+n1)
  ybar <- sum (c(yp,yu))/(n+n2)
  a11 <- sum (( xp - xbar1 )^2)
  a22 <- sum (( yp - ybar1 )^2)
  a12 <- sum (( xp - xbar1 ) * ( yp - ybar1 ))
  b1 <- sum (( xu - xbar2 )^2)
  b2 <- sum (( yu - ybar2 )^2)
  u <- (2*a12)/(a11 + a22)

  K1 <- (m1^2+ m2^2-2*m1*m2*u)/((1-m1)^2)
  K2 <- (m1^2+ m2^2-2*m1*m2*u)/((1-m2)^2)

  #checked
  denom <-(((m1^2 * a11 + m2^2 * a22 -2*m1*m2*a12 +
    K1*b1*(1-m1)^2 + K2*b2*(1-m2)^2)/(n+n1+n2-3)) *
      (1/n + 1/(K1*n1) + 1/(K2*n2)))^(1/2)
  num   <-  m1 * xbar1 - m2*ybar1 + (1-m1)*xbar2 - (1-m2)*ybar2

  Tb <- num/denom
  df <- n+n1+n2-3
  p_value <- pt(Tb, df, lower.tail = F) # ONE SIDED

  return(c(Tb,p_value))
}

LS_Z1s <- function(xp,xu,yp,yu,n,n1,n2){
  xbar1 <- mean ( xp )
  xbar2 <- mean ( xu )
  ybar1 <- mean ( yp )
  ybar2 <- mean ( yu )

  a11 <- sum (( xp - xbar1 )^2)
  a22 <- sum (( yp - ybar1 )^2)
  a12 <- sum (( xp - xbar1 ) * ( yp - ybar1 ))

  r <- a12/(sqrt(a11*a22))

  f1 <- n-1

  g <- n*(n+n2+ (n1*a12)/a11 )/( (n+n1)*(n+n2) - n1*n2*r^2) # ~ a_hat
  h <- n*(n+n1+ (n2*a12)/a22 )/( (n+n1)*(n+n2) - n1*n2*r^2) # ~ b_hat

  V1 <- ( (g^2/n + (1-g)^2/n1)*a11 + (h^2/n + (1-h)^2/n2)*a22 - (2*h*g*a12)/n)/f1 #checked

  Z1s <- (g*xbar1 + (1-g)*xbar2 - h*ybar1 - (1-h)*ybar2 )/(sqrt(V1))

  p_value <- pt(Z1s, n, lower.tail = F) # ONE SIDED

  return(c(Z1s, p_value))
}

"""## Derrick"""

Tnew1_derrick <- function(xp,xu,yp,yu,n,r){
  x <- c(xp,xu)
  y <- c(yp, yu)

  n_1 <- length(x)
  n_2 <- length(y)

  n_a <- length(xu)
  n_b <- length(yu)

  X_1 <- mean(x)
  X_2 <- mean(y)

  S_py <- sqrt(((n_1-1)*var(x) + (n_2-1)*var(y))/((n_1-1)+(n_2-1)))

  denom_Tnew1 <- S_py * sqrt(1/(n_1)+1/(n_2)-2*r*(n/(n_1*n_2)) )

  #compute statistic
  Tnew1 <- (X_1 - X_2)/denom_Tnew1

  #Reference distribution is T distribution with df= v_1
  v1 <- (n-1)+ ((n_a+n_b+n-1)/(n_a+n_b+2*n))*(n_a+n_b)

  #p_value
  #p_value <- 2*pt(-abs(Tnew1), v1) # TWO SIDED
  p_value <- pt(Tnew1, v1, lower.tail = F) # ONE SIDED

  return(c(Tnew1, p_value))
}

Tnew2_derrick <- function(xp,xu,yp,yu,n,r){
  x <- c(xp,xu)
  y <- c(yp, yu)

  n_1 <- length(x)
  n_2 <- length(y)

  n_a <- length(xu)
  n_b <- length(yu)

  X_1 <- mean(x)
  X_2 <- mean(y)

  denom_Tnew1 <- sqrt(var(x)/(n_1)+var(y)/(n_2)-2*r*((sqrt(var(x)*var(y))*n)/(n_1*n_2)) )

  #compute statistic
  Tnew2 <- (X_1 - X_2)/denom_Tnew1

  #Reference distribution is T distribution with df= v_1
  gamma_df <- (var(x)/n_1 + var(y)/n_2)^2/((var(x)/n_1)^2/(n_1-1) + (var(y)/n_2)^2/(n_2-1))
  v2 <- (n-1)+ ((gamma_df -n+1)/(n_a +n_b + 2*n))*(n_a+n_b)

  #get p_value from this reference distribution
  #p_value <- 2*pt(-abs(Tnew2), v2) #TWO SIDED
  p_value <- pt(Tnew2, v2, lower.tail = F)  # ONE SIDED

  return(c(Tnew2, p_value))
}

T_adj <- function(Tnew2_stat,xp,xu,yp,yu,n ){
  x <- c(xp,xu)
  y <- c(yp, yu)

  n_1 <- length(c(xp,xu))
  n_2 <- length(c(yp,yu))

  #satterthwaite df
  satterth_df <- (var(x)/n_1 + var(y)/n_2)^2/((var(x)/n_1)^2/(n_1-1) + (var(y)/n_2)^2/(n_2-1))
  paired_df <- n-1

  v2 <- satterth_df+paired_df

  #get p_value from this reference distribution
  p_value <- pt(Tnew2_stat, v2, lower.tail = F) # ONE SIDED

  return(p_value)
}

"""## Initialization"""

np <- 566
m <- 5000

#sigma_x is parameterized for gamma
sigma_y = 1

conf = 0.05
dist = " Normal "
col_names <- c(' Tnew1 ', ' Tnew2 ', ' T_Partover ',' T_adj ', " Zb ", ' Z1s ', " Tb ", " M", " R "," Rw ", " T test "," S ")

#Set the desired levels to be tested
delta_levels <- c(0, 0.5,1)
gamma_levels <- c(1,2,4)
cor_levels <- c(0.1, 0.3, 0.5, 0.9)

#Initialize power levels used in simulation loop
D = rep ( NA ,m )
Dw = rep ( NA , m )
M = rep ( NA ,m )
Robs = rep ( NA , m )
Rwobs = rep ( NA , m)
W = rep ( NA ,m )
t.stat = rep ( NA , m )
Tobs = rep ( NA , m )
Zb = rep ( NA , m )
pwr.t = 0
pwr.Z1s =0
pwr.r = 0
pwr.rw = 0
pwr.zb = 0
pwr.Tb = 0
pwr.m = 0
pwr.ttest = 0
pwr.tvar1 = 0
pwr.tvar2 = 0
pwr.d = 0
pwr.dw = 0
pwr.w = 0
pwr.package_test = 0
pwr.T_adj = 0

"""
# Outer data loop
"""

outer_data_loop <- function(samp, sampx, sampy, mu1 , mu2, sigma_x, cor){
for ( j in 1: m ) {
  dt <- data( samp , sampx , sampy , mu1 , mu2, sigma_x, sigma_y , cor , diff , dist )
  xp <- dt $ xp
  yp <- dt $ yp
  xu <- dt $ xu
  yu <- dt $ yu
  n <- dt $ n
  n1 <- dt $ n1
  n2 <- dt $ n2
  r <- dt $ r

  # Dubnicka
  if (( n1 != 0) && ( n2 != 0)){
    N = n1 + n2
    wr =2 / (( n * N )+(2 * n1 * n2 ))

    Sobs = wilcox.test( xp , yp , alternative = "greater" ,paired = TRUE )$statistic
    Uobs = wilcox.test( xu , yu , alternative = "greater" ,paired = FALSE )$ statistic
    Robs [j ]= Sobs + Uobs
    Rwobs [ j ]=(( N / ( n +1)) * wr * Sobs )+( wr * Uobs )
    }
  else {
      wr =1
      Sobs = wilcox.test( xp , yp , alternative = "greater" ,paired = TRUE )$statistic
      Robs [j ]= Sobs
      Rwobs [ j ]= wr * Sobs
    }

  # Magel
  mus <- n * ( n +1) / 4
  vars <- n * ( n +1) * (2 *n +1) / 24
  muu <- n1 * n2 / 2
  varu <- n1 * n2 * ( n1 + n2 +1) / 12
  if (( n1!= 0) && ( n2 != 0)) {
    M [ j ] = (1 / sqrt (2)) * ((( Sobs - mus ) / sqrt ( vars ))+
      (( Uobs - muu ) / sqrt ( varu )))
  } else M [ j ] = (( Sobs - mus )/ sqrt ( vars ))

  # Dubnicka asymptotic R, Rw
  D [ j ] = ( Robs [ j ] -( mus + muu )) / sqrt ( vars + varu )
  Sz = ( Sobs - mus ) / sqrt ( vars )
  if (( n1 != 0) && ( n2 != 0)) {
    Uz = ( Uobs - muu ) / sqrt ( varu )
    Dw [j ] = ( Rwobs [ j ] -(1 / 2)) /
    sqrt ((((( N / (n +1)) * wr )^2) * vars )+((( wr )^2) * varu ))
  } else Dw [ j] = D [ j]

  # Bhoj
  pvalue.zb <- Bhoj_Zb(xp,xu,yp,yu,n,n1,n2)[2]

  #Lin and Stivers Z1s
  LS_Z1s_pvalue <- LS_Z1s(xp,xu,yp,yu,n,n1,n2)[2]

  #Bhoj T statistic
  Bhoj_T_pvalue <- Bhoj_T(xp,xu,yp,yu,n,n1,n2)[2]

  # Paired t - test
  ttest <- t.test (xp , yp , alternative = "greater" ,paired = TRUE , var.equal =  (sigma_x == sigma_y)  ) #one sided
  t.stat [ j ] <- ttest $ statistic

  # Derrick statistics
  tvar1_pvalue <- Tnew1_derrick(xp,xu,yp,yu,n,r)[2]
  tvar2_pvalue <- Tnew2_derrick(xp,xu,yp,yu,n,r)[2]

  #Derrick test from package - returns p value based on Boolean variance equality
  package_test <- Partover.test(xu,yu,xp,yp,var.equal= (sigma_x == sigma_y) )$p.value

  #Proposed Tnew2 improvement
    #UNBIASED ESTIMATOR OF CORRELATION
  xbar1 <- mean ( xp )
  ybar1 <- mean ( yp )
  a11 <- sum (( xp - xbar1 )^2)
  a22 <- sum (( yp - ybar1 )^2)
  a12 <- sum (( xp - xbar1 ) * ( yp - ybar1 ))
  cov <- a12/(sqrt(a11*a22))
  rho_est <- cov*(1+ (1-cov^2)/(2*(n-3)))

  Tnew2_stat <- Tnew2_derrick(xp,xu,yp,yu,n, rho_est)[1]
  T_adj_pvalue <- T_adj(Tnew2_stat,xp,xu,yp,yu,n )

  # Signed Rank test
  W.test <- wilcox.test ( xp ,yp ,alternative ="greater" , paired = TRUE ) # one sided
  W [ j ] <- ( W.test $ statistic - mus ) / sqrt ( vars )

  # Inner Permutation loop
  t10 <- rep ( NA , np )
  R <- rep ( NA , np )
  Rw <- rep (NA , np )
  pt = 0
  pr = 0
  prw = 0
  wt =(1 / n2 +1 / n1 ) / ((2 -2 * r ) /n +(1 / n2 +1 / n1 ))


  # P VALUES
  pvalue.t <- pt / np
  pvalue.r <- pr / np
  pvalue.rw <- prw / np
  pvalue.m <- pnorm ( M[ j ] , lower.tail = FALSE )
  pvalue.ttest <- ttest $ p.value
  pvalue.d <- pnorm ( D[ j ] , lower.tail = FALSE )
  pvalue.dw <- pnorm ( Dw [ j ] , lower.tail = FALSE )
  pvalue.w <- pnorm ( W[ j ] , lower.tail = FALSE )

  #TYPE 1 ERROR / POWER TESTING
  if ( LS_Z1s_pvalue < conf) pwr.Z1s = pwr.Z1s +1
  if ( pvalue.zb < conf) pwr.zb = pwr.zb +1
  if ( Bhoj_T_pvalue < conf) pwr.Tb = pwr.Tb +1
  if ( pvalue.t < conf) pwr.t = pwr.t +1
  if ( pvalue.r < conf) pwr.r = pwr.r +1                # r perm
  if ( pvalue.rw < conf) pwr.rw = pwr.rw +1             # r w perm
  if ( pvalue.m < conf) pwr.m = pwr.m +1
  if ( pvalue.ttest < conf) pwr.ttest = pwr.ttest +1
  if ( tvar1_pvalue <  conf) pwr.tvar1 = pwr.tvar1 +1   #Derrick Tvar1
  if ( tvar2_pvalue <  conf) pwr.tvar2 = pwr.tvar2 +1   #Derrick Tvar2
  if ( package_test <  conf) pwr.package_test = pwr.package_test +1   #Derrick Tvar package
  if ( T_adj_pvalue <  conf) pwr.T_adj = pwr.T_adj +1                   #Derrick Tvar2 improvement
  if ( pvalue.d < conf) pwr.d = pwr.d +1
  if ( pvalue.dw < conf) pwr.dw = pwr.dw +1
  if ( pvalue.w < conf) pwr.w = pwr.w +1
  }


  return(list (
  'Delta' =  mu1,
  'Gamma' =sigma_x,
  'Rho' = cor,
  ' Tnew1 ' = pwr.tvar1 / m,
  ' Tnew2 ' = pwr.tvar2 / m,
  ' T_Partover ' = pwr.package_test / m,
  ' T_adj ' = pwr.T_adj / m,
  " Zb " = pwr.zb /m ,
  ' Z1s ' =  pwr.Z1s/m,
  " Tb " = pwr.Tb /m ,
  " M" = pwr.m/m ,
  " R " = pwr.d /m ,
  " Rw " = pwr.dw /m ,
  " T test " = pwr.ttest /m ,
  " S " = pwr.w / m))
}

"""# execute

## Execute loop
"""

run_simulation <- function(samp,sampx,sampy){  #Can also initialize:  delta_levels, gamma_levels,cor_levels
  results.compile <- data.frame()

    #executing loop
  for(mu1 in delta_levels){
    for(sigma_x in gamma_levels){

        results <- data.frame(matrix(ncol = length(col_names) + 3, nrow = length(cor_levels)))  # Initialize one-run  dataframe, cases x (tests + 3 labels)

        for(cor in  cor_levels){
          set.seed(2018)
          results[match(cor, cor_levels), ] <- outer_data_loop(samp, sampx, sampy, mu1, 0, sigma_x, cor)
          }
      colnames(results) <- c( 'Delta','Gamma', 'Rho', col_names)
      results.compile <- rbind(results.compile, results)
    }
  }
  return(results.compile)
}

#The final cases used for the associated paper were as follows:
final_cases <- list(c(5,5,10),
c(5,10,5),
c(5,7,8) ,
c(5,8,7),
c(10,5,5),
c(10,7,3) ,
c(10,3,7) ,
c(15,2,3),
c(15,3,2) ,
c(15,15,30),
c(15,30,15) ,
c(15,22,23) ,
c(30,15,15) ,
c(30,21,9),
c(30,9,21) ,
c(45,7,8) ,
c(45,8,7) )

"""**Run a single example**"""

results <- run_simulation(15,15,30)
results

"""**Run a list of sample size cases and save to xlsx file**"""

#Name the exporting file
path <- 'runagain.xlsx'

#Initialize cases list
selected_cases = final_cases
#selected_cases = list(c() ) #format: list(c(1,1,1), c(2,2,2), ... )

#Create workbook
wb <- createWorkbook()

sheet_counter <- 1
for(set in selected_cases){
  #Rests the loop for 15 seconds every 4 runs

  #Splits the sample size tuple into parts
  set <- unlist(set)
  nc<- set[1]
  nux<- set[2]
  nuy<- set[3]

  name <-  paste('Sheet',sheet_counter,'-', as.character(nc),as.character(nux),as.character(nuy))


  addWorksheet(wb, name)
  results <- run_simulation(nc,nux,nuy) #Execute the simulation for the sample settings
  writeData(wb, name, results)          #Write the data to the Excel sheet

  #Print how long the execution took

  #Iterate the next index
  sheet_counter <- sheet_counter + 1
}

#Result will be saved as a workbook with path name
saveWorkbook(wb, path)
