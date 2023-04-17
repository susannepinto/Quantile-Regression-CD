Heterogeneous associations of gut microbiota with Crohn’s disease
activity
================
2023-04-14

## Introduction

R code belonging to the paper: “Heterogeneous associations of gut
microbiota with Crohn’s disease activity” by Susanne Pinto, Elisa
Benincà, Gianluca Galazzo, Daisy Jonkers, John Penders and Johannes
Bogaards.

We investigated the relationship between bacterial relative abundances
and disease activity in a longitudinal cohort of CD patients (n = 57)
and healthy controls (n = 15). We applied quantile regression, a
statistical technique that allows investigation of possible
relationships between bacterial relative abundance and Crohn’s disease
activity outside the mean response. For this example script we have used

## Load the packages

Load the packages in R.

``` r
library( tidyverse )
library( here )
library( lqmm )
library( nlme )
library( sjPlot )
```

# The data

The data need to be in a long dataframe format. Here an example dataset
for a single family is given. Prior to the analyses, relative abundances
were multiplied with 1000 and log-transformed with the natural log
function assuming a lower detection limit of 100 reads (which is 1/4th
of the lowest measurable value in the real dataset). Also the age
variable was transformed (centered around mean age) else you will run
into convergence errors in the model.

``` r
IBD.data.Family <- data.frame( individual = c( do.call( paste0, expand.grid( c( 'HC' ), 1:15 )),
                                               do.call( paste0, expand.grid( c( 'IBD' ), 1:57 )),
                                               do.call( paste0, expand.grid( c( 'HC' ), 1:15 )),
                                               do.call( paste0, expand.grid( c( 'IBD' ), 1:57 ))),
                               visit = c( rep( 1, 72 ), rep( 2, 72 )),
                               Group_01 <- c( rep( 0, 15 ), rep( 1, 35 ), rep( 2, 22 ),
                                             rep( 0, 15 ), rep( 1, 35 ), rep( 2, 22 )),
                               Status_01 <- c( rep( 0, 15 ), rep( 1, 35 ), rep( 1, 22 ),
                                             rep( 0, 15 ), rep( 1, 35 ), rep( 2, 22 )),
                               age <- runif( 72, 20, 70 ),
                               Gender <- rep( sample( c( rep( c( "male", "female" ), 36 ))), 2 ),
                               Smoking <- rep( sample( c( rep( c( "ex", "never", "current" ), 24 ))), 2 ),
                               density <- abs( rnorm( 144, mean = 0, sd = 1 )))

colnames( IBD.data.Family ) <- c( "individual", "visit", "Group_01", "Status_01", "age", "Gender", "Smoking", "density" )
IBD.data.Family$visit <- as.factor( IBD.data.Family$visit )
IBD.data.Family$Group_01 <- as.factor( IBD.data.Family$Group_01 )
IBD.data.Family$Status_01 <- as.factor( IBD.data.Family$Status_01 )

# Transform age
mean.age <- IBD.data.Family$age %>% mean()
IBD.data.Family$age2 <- as.numeric( IBD.data.Family$age - mean.age )

# Transform density
IBD.data.Family$density.1000 <- IBD.data.Family$density * 1000
IBD.data.Family$density.log1000 <- ifelse( IBD.data.Family$density.1000 > 0,
                              log( IBD.data.Family$density.1000 ), log( 100 ))

str(IBD.data.Family)
```

    ## 'data.frame':    144 obs. of  11 variables:
    ##  $ individual     : chr  "HC1" "HC2" "HC3" "HC4" ...
    ##  $ visit          : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Group_01       : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Status_01      : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ age            : num  32.5 66.6 37.2 49 45.8 ...
    ##  $ Gender         : chr  "male" "female" "male" "male" ...
    ##  $ Smoking        : chr  "never" "current" "current" "ex" ...
    ##  $ density        : num  0.8182 0.0946 1.0476 0.286 0.222 ...
    ##  $ age2           : num  -12.834 21.323 -8.104 3.655 0.487 ...
    ##  $ density.1000   : num  818.2 94.6 1047.6 286 222 ...
    ##  $ density.log1000: num  6.71 4.55 6.95 5.66 5.4 ...

# The smallest model

For the first model, we investigated whether the relative abundances of
the bacterial families could be explained by the group to which each
individual belongs (i.e., healthy control (HC), remission-remission
(RR), or remission-exacerbation (RE)). We added the interaction with
visit number, to allow for different temporal changes in bacterial
relative abundance over time between healthy controls, CD patients who
experienced an exacerbation at the second visit, and those who remained
in remission.

The models contain two timepoints per individual. Therefore, we used a
random intercept per patient as well as a random effect for the variable
‘visit number’, because temporal changes in bacterial family’s relative
abundance may differ within patients, even when accounting for the fixed
effect of disease trajectory (e.g., experiencing an exacerbation at the
second visit).

``` r
# The minimal model
# Model selection will only be run on the 50% quantile
fit.model1 <- lqmm( density.log1000 ~  Group_01 * visit, random = ~visit, 
                    # group = individual, because two time points
                    group = individual, tau = 0.5, covariance = "pdSymm",
                    data = IBD.data.Family, na.action = na.omit,
                    control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
fit.model1
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.5, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.5 
    ## 
    ## Fixed effects:
    ##      (Intercept)         Group_011         Group_012            visit2  
    ##           6.0360            0.6004            0.4393            0.8291  
    ## Group_011:visit2  Group_012:visit2  
    ##          -1.3430           -0.7513  
    ## 
    ## Covariance matrix of the random effects:
    ##             (Intercept) visit2   
    ## (Intercept)  0.119606   -0.008205
    ## visit2      -0.008205    0.120040
    ## 
    ## Residual scale parameter: 0.3172 (standard deviation 0.8973)
    ## Log-likelihood: -197.9 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

``` r
sum.fit.model1 <- summary( fit.model1 )
sum.fit.model1
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.5, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.5 
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.035961   0.229982    5.573795      6.4981 < 2.2e-16 ***
    ## Group_011         0.600442   0.289254    0.019165      1.1817   0.04318 *  
    ## Group_012         0.439298   0.341990   -0.247956      1.1266   0.20499    
    ## visit2            0.829122   0.176449    0.474535      1.1837 2.154e-05 ***
    ## Group_011:visit2 -1.343015   0.296848   -1.939554     -0.7465 3.868e-05 ***
    ## Group_012:visit2 -0.751264   0.361222   -1.477167     -0.0254   0.04280 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## AIC:
    ## [1] 415.7 (df = 10)

# Model selection

Variable selection was performed by running all possible models and then
selecting the model with the lowest Bayesian Information Criterion (BIC)
in the 50% quantile.

``` r
# Covariates for selection procedure (Gender, smoking and transformed age)
vars = c( "Gender", "Smoking" ,"age2")

# Create a NULL vector called model so we have something to add our layers to
model.group = NULL
# First add the minimal model without covariates
model.group[[1]] <- fit.model1

for( i in 1:length( vars )) {
  # The combn function will run every different combination of variables and then run the model
  # Below function makes every possible combination of the variables
  xx = combn( vars, i )
  for( j in 1:dim( xx )[2] ){
    fla = paste( "density.log1000 ~  Group_01 * visit +", paste( xx[ 1:dim( xx )[ 1 ], j ], collapse = "+" ))
    model.group[[ length( model.group ) + 1 ]] = lqmm( as.formula( fla ),
                                                       random = ~visit, group = individual, tau = 0.5, 
                                                       covariance = "pdSymm", data = IBD.data.Family, na.action = na.omit,
                                                       control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
    # print( j )
  }
  # print( i )
}

# Make a list of all AICs and BICs
AICs = NULL
BICs = NULL
for(k in 1:length( model.group )){
  AICs[[k]] = AIC( model.group[[k]] ) 
  BICs[[k]] = BIC( model.group[[k]] ) 
}

AICs <- AICs %>% unlist()
AICs
```

    ## [1] 415.7214 415.8723 415.4665 424.5069 417.9106 428.1015 422.6746 424.4939

``` r
# See which model was chosen based on the lowest AIC
AIC <- which( AICs == min( AICs ))
AIC
```

    ## [1] 3

``` r
BICs <- BICs %>% unlist()
BICs
```

    ## [1] 445.4195 448.5402 451.1043 457.1748 456.5182 463.7393 461.2821 466.0713

``` r
# See which model was chosen based on the lowest AIC
BIC <- which( BICs == min( BICs ))
BIC
```

    ## [1] 1

``` r
# We use the model with the lowest AIC
selected.model.group <- model.group[[ BIC ]] 
selected.model.group
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.5, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.5 
    ## 
    ## Fixed effects:
    ##      (Intercept)         Group_011         Group_012            visit2  
    ##           6.0360            0.6004            0.4393            0.8291  
    ## Group_011:visit2  Group_012:visit2  
    ##          -1.3430           -0.7513  
    ## 
    ## Covariance matrix of the random effects:
    ##             (Intercept) visit2   
    ## (Intercept)  0.119606   -0.008205
    ## visit2      -0.008205    0.120040
    ## 
    ## Residual scale parameter: 0.3172 (standard deviation 0.8973)
    ## Log-likelihood: -197.9 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

# The selected model

``` r
# run the model with for more quantiles
qs <- seq(0.1, 0.9, by = 0.1) # seq(0.05, 0.95, by = 0.05)

# run the model with the selected covariates (if present)
fit.model1.2 <- lqmm( density.log1000 ~  Group_01 * visit, random = ~visit, 
                      # group = individual, because two time points
                      group = individual, tau = qs, covariance = "pdSymm",
                      data = IBD.data.Family, na.action = na.omit,
                      control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
fit.model1.2
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = qs, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Fixed effects:
    ##                   tau = 0.1  tau = 0.2  tau = 0.3  tau = 0.4  tau = 0.5
    ## (Intercept)        5.1957     5.4325     5.6694     5.9348     6.0360  
    ## Group_011          0.1033     0.3503     0.5397     0.4748     0.6004  
    ## Group_012         -0.1210     0.1329     0.2864     0.3333     0.4393  
    ## visit2             0.9565     0.8955     0.9232     0.8114     0.8291  
    ## Group_011:visit2  -1.3453    -1.5518    -1.3474    -1.3825    -1.3430  
    ## Group_012:visit2  -0.7377    -0.6628    -0.7499    -0.7308    -0.7513  
    ##                   tau = 0.6  tau = 0.7  tau = 0.8  tau = 0.9
    ## (Intercept)        6.1767     6.3782     6.5084     6.7214  
    ## Group_011          0.6691     0.6637     0.6343     0.5928  
    ## Group_012          0.4706     0.5140     0.6158     0.5918  
    ## visit2             0.8177     0.8110     0.7382     0.6584  
    ## Group_011:visit2  -1.2725    -1.2382    -1.0699    -1.0905  
    ## Group_012:visit2  -0.7409    -0.7443    -0.6744    -0.6826  
    ## 
    ## Covariance matrix of the random effects:
    ## tau = 0.1
    ##             (Intercept) visit2  
    ## (Intercept)  0.96146    -0.04853
    ## visit2      -0.04853     0.15077
    ## tau = 0.2
    ##             (Intercept) visit2 
    ## (Intercept)  0.7681     -0.2308
    ## visit2      -0.2308      0.5122
    ## tau = 0.3
    ##             (Intercept) visit2   
    ## (Intercept)  0.187000   -0.004086
    ## visit2      -0.004086    0.187151
    ## tau = 0.4
    ##             (Intercept) visit2  
    ## (Intercept)  0.24412    -0.05755
    ## visit2      -0.05755     0.25759
    ## tau = 0.5
    ##             (Intercept) visit2   
    ## (Intercept)  0.119606   -0.008205
    ## visit2      -0.008205    0.120040
    ## tau = 0.6
    ##             (Intercept) visit2 
    ## (Intercept) 0.08650     0.03178
    ## visit2      0.03178     0.01642
    ## tau = 0.7
    ##             (Intercept) visit2   
    ## (Intercept)  0.136103   -0.009873
    ## visit2      -0.009873    0.012189
    ## tau = 0.8
    ##             (Intercept) visit2  
    ## (Intercept)  0.14345    -0.05172
    ## visit2      -0.05172     0.07617
    ## tau = 0.9
    ##             (Intercept) visit2 
    ## (Intercept)  0.1875     -0.0607
    ## visit2      -0.0607      0.1538
    ## 
    ## Residual scale parameter: 0.09859 (tau = 0.1)  0.16881 (tau = 0.2)  0.28397 (tau = 0.3)  0.29224 (tau = 0.4)  0.31725 (tau = 0.5)  0.30581 (tau = 0.6)  0.25685 (tau = 0.7)  0.18773 (tau = 0.8)  0.09277 (tau = 0.9) 
    ## Log-likelihood: -228.7 (tau = 0.1)  -223.2 (tau = 0.2)  -215.3 (tau = 0.3)  -204.4 (tau = 0.4)  -197.9 (tau = 0.5)  -193.9 (tau = 0.6)  -190.6 (tau = 0.7)  -189.7 (tau = 0.8)  -189.3 (tau = 0.9) 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

``` r
sum.fit.model1.2 <- summary( fit.model1.2 )
sum.fit.model1.2 
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = qs, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## tau = 0.1
    ## 
    ## Fixed effects:
    ##                     Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       5.19571    0.22911     4.73529      5.6561 < 2.2e-16 ***
    ## Group_011         0.10325    0.34304    -0.58612      0.7926 0.7646988    
    ## Group_012        -0.12099    0.41204    -0.94901      0.7070 0.7702813    
    ## visit2            0.95650    0.24137     0.47146      1.4415 0.0002402 ***
    ## Group_011:visit2 -1.34530    0.36523    -2.07925     -0.6113 0.0005741 ***
    ## Group_012:visit2 -0.73765    0.49139    -1.72513      0.2498 0.1397321    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.2
    ## 
    ## Fixed effects:
    ##                     Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       5.43252    0.20171     5.02717      5.8379 < 2.2e-16 ***
    ## Group_011         0.35032    0.27971    -0.21178      0.9124 0.2163584    
    ## Group_012         0.13290    0.32343    -0.51705      0.7829 0.6829199    
    ## visit2            0.89554    0.24330     0.40661      1.3845 0.0005788 ***
    ## Group_011:visit2 -1.55178    0.33828    -2.23159     -0.8720 3.134e-05 ***
    ## Group_012:visit2 -0.66278    0.36827    -1.40285      0.0773 0.0780624 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.3
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       5.669394   0.180978    5.305706      6.0331 < 2.2e-16 ***
    ## Group_011         0.539655   0.247554    0.042177      1.0371 0.0340947 *  
    ## Group_012         0.286433   0.282270   -0.280809      0.8537 0.3152103    
    ## visit2            0.923157   0.203465    0.514279      1.3320 3.704e-05 ***
    ## Group_011:visit2 -1.347416   0.344533   -2.039781     -0.6551 0.0002831 ***
    ## Group_012:visit2 -0.749851   0.317753   -1.388400     -0.1113 0.0223095 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.4
    ## 
    ## Fixed effects:
    ##                        Value  Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       5.93484946  0.17612896  5.58090506      6.2888 < 2.2e-16 ***
    ## Group_011         0.47484007  0.23616351  0.00025174      0.9494   0.04988 *  
    ## Group_012         0.33328624  0.25680453 -0.18278178      0.8494   0.20043    
    ## visit2            0.81143574  0.18523868  0.43918468      1.1837 6.226e-05 ***
    ## Group_011:visit2 -1.38254824  0.30401934 -1.99349797     -0.7716 3.578e-05 ***
    ## Group_012:visit2 -0.73079315  0.29351461 -1.32063284     -0.1410   0.01622 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.5
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.035961   0.154277    5.725929      6.3460 < 2.2e-16 ***
    ## Group_011         0.600442   0.217429    0.163503      1.0374  0.008075 ** 
    ## Group_012         0.439298   0.263437   -0.090099      0.9687  0.101782    
    ## visit2            0.829122   0.183615    0.460134      1.1981 3.981e-05 ***
    ## Group_011:visit2 -1.343015   0.299302   -1.944484     -0.7415 4.375e-05 ***
    ## Group_012:visit2 -0.751264   0.279709   -1.313361     -0.1892  0.009847 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.6
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.176718   0.153422    5.868406      6.4850 < 2.2e-16 ***
    ## Group_011         0.669103   0.209069    0.248964      1.0892  0.002409 ** 
    ## Group_012         0.470633   0.272834   -0.077648      1.0189  0.090834 .  
    ## visit2            0.817726   0.188443    0.439036      1.1964 7.127e-05 ***
    ## Group_011:visit2 -1.272520   0.286847   -1.848960     -0.6961 5.180e-05 ***
    ## Group_012:visit2 -0.740941   0.286241   -1.316164     -0.1657  0.012652 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.7
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.378186   0.183659    6.009110      6.7473 < 2.2e-16 ***
    ## Group_011         0.663708   0.218483    0.224650      1.1028 0.0038139 ** 
    ## Group_012         0.514011   0.263762   -0.016039      1.0441 0.0570612 .  
    ## visit2            0.810998   0.192620    0.423913      1.1981 0.0001086 ***
    ## Group_011:visit2 -1.238198   0.261088   -1.762874     -0.7135  1.86e-05 ***
    ## Group_012:visit2 -0.744311   0.284888   -1.316814     -0.1718 0.0118959 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.8
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.508390   0.212410    6.081536      6.9352 < 2.2e-16 ***
    ## Group_011         0.634283   0.237365    0.157280      1.1113 0.0102038 *  
    ## Group_012         0.615761   0.289068    0.034857      1.1967 0.0382046 *  
    ## visit2            0.738242   0.197551    0.341249      1.1352 0.0004869 ***
    ## Group_011:visit2 -1.069915   0.227361   -1.526814     -0.6130 2.105e-05 ***
    ## Group_012:visit2 -0.674357   0.293675   -1.264520     -0.0842 0.0259764 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## tau = 0.9
    ## 
    ## Fixed effects:
    ##                     Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.72135    0.25159     6.21576      7.2269 < 2.2e-16 ***
    ## Group_011         0.59280    0.23262     0.12532      1.0603  0.014010 *  
    ## Group_012         0.59176    0.28500     0.01903      1.1645  0.043129 *  
    ## visit2            0.65837    0.19498     0.26654      1.0502  0.001445 ** 
    ## Group_011:visit2 -1.09055    0.20828    -1.50911     -0.6720 3.429e-06 ***
    ## Group_012:visit2 -0.68261    0.31901    -1.32369     -0.0415  0.037382 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## AIC:
    ## [1] 477.4 (df = 10) 466.5 (df = 10) 450.6 (df = 10) 428.9 (df = 10)
    ## [5] 415.7 (df = 10) 407.7 (df = 10) 401.2 (df = 10) 399.5 (df = 10)
    ## [9] 398.5 (df = 10)

# Run the model 10 times

We run the selected model 10 times and use the median values, because
the p-values slightly change every time. The lqmm package uses the
optimization algorithm “Broyden-Fletcher-Goldfarb-Shanno” (BFGS) to
estimate the parameters of the linear quantile mixed model. The BFGS
algorithm is an iterative method that relies on an initial starting
point and a set of convergence criteria to estimate the parameters.

Because the BFGS algorithm is an iterative method, the estimates
obtained can depend on the starting values and the convergence criteria
used. If the starting values are very different, or if the convergence
criteria are not tight enough, then the algorithm may converge to
different estimates on different runs.

``` r
allsumms.model1 = NULL
for(l in 1:10){
  allsumms.model1[[l]] = summary( fit.model1.2 )
  print(l)
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10

``` r
# Subtract the table with the results
all.ttables.model1 = NULL
for(m in 1:length(allsumms.model1)){
  all.ttables.model1[[m]] = allsumms.model1[[m]][["tTable"]]
  #print(m)
}

# Subtract the values per quantile
all.values.model1 = NULL
for(n in 1:length(all.ttables.model1)){
  all.values.model1[["0.1"]][[n]] = all.ttables.model1[[n]][["0.1"]]
  all.values.model1[["0.2"]][[n]] = all.ttables.model1[[n]][["0.2"]]
  all.values.model1[["0.3"]][[n]] = all.ttables.model1[[n]][["0.3"]]
  all.values.model1[["0.4"]][[n]] = all.ttables.model1[[n]][["0.4"]]
  all.values.model1[["0.5"]][[n]] = all.ttables.model1[[n]][["0.5"]]
  all.values.model1[["0.6"]][[n]] = all.ttables.model1[[n]][["0.6"]]
  all.values.model1[["0.7"]][[n]] = all.ttables.model1[[n]][["0.7"]]
  all.values.model1[["0.8"]][[n]] = all.ttables.model1[[n]][["0.8"]]
  all.values.model1[["0.9"]][[n]] = all.ttables.model1[[n]][["0.9"]]
  #print(n)
}

# Make a list with the different results
nrows <- all.values.model1[["0.1"]][[1]] %>% nrow()
median.table.model1 = list()

for( p in 1:nrows ){
  median.table.model1[["0.1"]][[p]] <- list()
  median.table.model1[["0.2"]][[p]] <- list()
  median.table.model1[["0.3"]][[p]] <- list()
  median.table.model1[["0.4"]][[p]] <- list()
  median.table.model1[["0.5"]][[p]] <- list()
  median.table.model1[["0.6"]][[p]] <- list()
  median.table.model1[["0.7"]][[p]] <- list()
  median.table.model1[["0.8"]][[p]] <- list()
  median.table.model1[["0.9"]][[p]] <- list()
  
  for(o in 1:length(all.ttables.model1)){
    median.table.model1[["0.1"]][[p]][[o]] = all.values.model1[["0.1"]][[o]][p,] %>% as.matrix() 
    median.table.model1[["0.2"]][[p]][[o]] = all.values.model1[["0.2"]][[o]][p,] %>% as.matrix() 
    median.table.model1[["0.3"]][[p]][[o]] = all.values.model1[["0.3"]][[o]][p,] %>% as.matrix()
    median.table.model1[["0.4"]][[p]][[o]] = all.values.model1[["0.4"]][[o]][p,] %>% as.matrix() 
    median.table.model1[["0.5"]][[p]][[o]] = all.values.model1[["0.5"]][[o]][p,] %>% as.matrix()
    median.table.model1[["0.6"]][[p]][[o]] = all.values.model1[["0.6"]][[o]][p,] %>% as.matrix()
    median.table.model1[["0.7"]][[p]][[o]] = all.values.model1[["0.7"]][[o]][p,] %>% as.matrix()
    median.table.model1[["0.8"]][[p]][[o]] = all.values.model1[["0.8"]][[o]][p,] %>% as.matrix()
    median.table.model1[["0.9"]][[p]][[o]] = all.values.model1[["0.9"]][[o]][p,] %>% as.matrix()
  } 
}

# Get the median values per quantile
median.values.model1 = NULL
for( q in 1:nrows ){
  median.values.model1[["0.1"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.1"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.2"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.2"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.3"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.3"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.4"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.4"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.5"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.5"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.6"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.6"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.7"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.7"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.8"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.8"]][[q]]), 1, median, na.rm=T)
  median.values.model1[["0.9"]][[q]] <- apply(do.call(cbind, median.table.model1[["0.9"]][[q]]), 1, median, na.rm=T)
}

# Make new results tables with the median values
rownames.results <- all.values.model1[["0.1"]][[1]] %>% rownames()
results.model1 = NULL
results.model1[["0.1"]] <- do.call(rbind, median.values.model1[["0.1"]])
rownames(results.model1[["0.1"]]) <- rownames.results

results.model1[["0.2"]] <- do.call(rbind, median.values.model1[["0.2"]])
rownames(results.model1[["0.2"]]) <- rownames.results

results.model1[["0.3"]] <- do.call(rbind, median.values.model1[["0.3"]])
rownames(results.model1[["0.3"]]) <- rownames.results

results.model1[["0.4"]] <- do.call(rbind, median.values.model1[["0.4"]])
rownames(results.model1[["0.4"]]) <- rownames.results

results.model1[["0.5"]] <- do.call(rbind, median.values.model1[["0.5"]])
rownames(results.model1[["0.5"]]) <- rownames.results

results.model1[["0.6"]] <- do.call(rbind, median.values.model1[["0.6"]])
rownames(results.model1[["0.6"]]) <- rownames.results

results.model1[["0.7"]] <- do.call(rbind, median.values.model1[["0.7"]])
rownames(results.model1[["0.7"]]) <- rownames.results

results.model1[["0.8"]] <- do.call(rbind, median.values.model1[["0.8"]])
rownames(results.model1[["0.8"]]) <- rownames.results

results.model1[["0.9"]] <- do.call(rbind, median.values.model1[["0.9"]])
rownames(results.model1[["0.9"]]) <- rownames.results
```

# Plot p-values model 1

Differences at baseline (visit 1) are visualized in light green (RR) and
red (RE), the interaction with visit number (in dark green and orange)
displays changes over time. Significant variables are indicated with a
closed circle.

``` r
fam.name <- "Example family"
values.model1 <- matrix(ncol = 7, nrow =  0) %>% as.data.frame()
colnames(values.model1) <- c("Family", "Variable", "Quantile", "Estimate", "lower", "upper", "Pvalue")

for( s in 1:length(results.model1 )){
  t.table.qs <- results.model1[[s]]
  values.model1[ nrow(values.model1) + 1 , ] <- NA
  
  values.model1[nrow(values.model1), 1] <- fam.name
  values.model1[nrow(values.model1), 2] <- "RR"
  values.model1[nrow(values.model1), 3] <- qs[s] * 100
  values.model1[nrow(values.model1), 4] <- t.table.qs["Group_011",1]
  values.model1[nrow(values.model1), 5] <- t.table.qs["Group_011",3]
  values.model1[nrow(values.model1), 6] <- t.table.qs["Group_011",4]
  values.model1[nrow(values.model1), 7] <- t.table.qs["Group_011",5]
}

for( s in 1:length(results.model1 )){
  t.table.qs <- results.model1[[s]]
  values.model1[ nrow(values.model1) + 1 , ] <- NA
  
  values.model1[nrow(values.model1), 1] <- fam.name
  values.model1[nrow(values.model1), 2] <- "RA"
  values.model1[nrow(values.model1), 3] <- qs[s] * 100
  values.model1[nrow(values.model1), 4] <- t.table.qs["Group_012",1]
  values.model1[nrow(values.model1), 5] <- t.table.qs["Group_012",3]
  values.model1[nrow(values.model1), 6] <- t.table.qs["Group_012",4]
  values.model1[nrow(values.model1), 7] <- t.table.qs["Group_012",5]
}

for( s in 1:length(results.model1 )){
  t.table.qs <- results.model1[[s]]
  values.model1[ nrow(values.model1) + 1 , ] <- NA
  
  values.model1[nrow(values.model1), 1] <- fam.name
  values.model1[nrow(values.model1), 2] <- "RR:visit"
  values.model1[nrow(values.model1), 3] <- qs[s] * 100
  values.model1[nrow(values.model1), 4] <- t.table.qs["Group_011:visit2",1]
  values.model1[nrow(values.model1), 5] <- t.table.qs["Group_011:visit2",3]
  values.model1[nrow(values.model1), 6] <- t.table.qs["Group_011:visit2",4]
  values.model1[nrow(values.model1), 7] <- t.table.qs["Group_011:visit2",5]
}

for( s in 1:length(results.model1 )){
  t.table.qs <- results.model1[[s]]
  values.model1[ nrow(values.model1) + 1 , ] <- NA
  
  values.model1[nrow(values.model1), 1] <- fam.name
  values.model1[nrow(values.model1), 2] <- "RA:visit"
  values.model1[nrow(values.model1), 3] <- qs[s] * 100
  values.model1[nrow(values.model1), 4] <- t.table.qs["Group_012:visit2",1]
  values.model1[nrow(values.model1), 5] <- t.table.qs["Group_012:visit2",3]
  values.model1[nrow(values.model1), 6] <- t.table.qs["Group_012:visit2",4]
  values.model1[nrow(values.model1), 7] <- t.table.qs["Group_012:visit2",5]
}

values.model1 %>% head()
```

    ##           Family Variable Quantile  Estimate       lower     upper     Pvalue
    ## 1 Example family       RR       10 0.1032508 -0.68449099 0.8909925 0.79334773
    ## 2 Example family       RR       20 0.3503192 -0.32860945 1.0292479 0.30478908
    ## 3 Example family       RR       30 0.5396547 -0.07705765 1.1563671 0.08497246
    ## 4 Example family       RR       40 0.4748401 -0.07727520 1.0269553 0.09082648
    ## 5 Example family       RR       50 0.6004424  0.08418358 1.1167013 0.02461453
    ## 6 Example family       RR       60 0.6691028  0.15082094 1.1873847 0.01298380

``` r
values.plot.model1 <- ggplot( values.model1, 
                              aes( x = Quantile, y = Estimate, 
                                   group = Variable, ymin = lower, ymax = upper)) +  
  theme_bw() +
  geom_line( aes( alpha = 0.4), size = 1.2 )+
  geom_point( aes(shape = Pvalue < 0.05, color = Variable, size = Pvalue < 0.05), size = 8, stroke = 2  ) +
  scale_shape_manual(name = 'P value < 0.05', values = setNames( c( 19, 1 ), c( T, F ))) +
  scale_color_manual(values = c("RR" = "Green", "RA" = "Red", 
                                "RR:visit" = "Dark green", "RA:visit" = "Orange")) +
  theme( legend.position = "none",
         axis.text=element_text( size = 12 ),
         axis.title.x = element_text(size = 14),
         axis.title.y = element_text(size = 14)#, 
         #strip.text.x = element_blank()
  ) + 
  labs(title = "Effect estimates Group") + 
  geom_hline( yintercept = 0, linetype = "dashed", color = "Black", size = 1 )   +
  ggtitle( fam.name )  
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

``` r
values.plot.model1
```

![](README_files/figure-gfm/Plot%20p-values-1.png)<!-- -->

# Ordinary linear mixed effect model

We compared our results from the lqmm models with the results obtained
from an ordinary linear mixed effect model (with similar variables as
used in the lqmm models) by using the lme function from the nlme package
(version 3.1).

``` r
model1 <- lme( density.log1000 ~ Group_01*visit, random = ~1+visit|individual, data = IBD.data.Family,
               na.action = na.omit )
summary( model1 ) 
```

    ## Linear mixed-effects model fit by REML
    ##   Data: IBD.data.Family 
    ##        AIC      BIC    logLik
    ##   444.8429 474.1155 -212.4215
    ## 
    ## Random effects:
    ##  Formula: ~1 + visit | individual
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##             StdDev    Corr  
    ## (Intercept) 0.9585886 (Intr)
    ## visit2      1.2560100 -0.678
    ## Residual    0.4738403       
    ## 
    ## Fixed effects:  density.log1000 ~ Group_01 * visit 
    ##                      Value Std.Error DF   t-value p-value
    ## (Intercept)       6.001941 0.2760938 69 21.738773  0.0000
    ## Group_011         0.447089 0.3299952 69  1.354834  0.1799
    ## Group_012         0.323136 0.3580521 69  0.902483  0.3699
    ## visit2            0.867274 0.3675695 69  2.359483  0.0211
    ## Group_011:visit2 -1.320037 0.4393296 69 -3.004661  0.0037
    ## Group_012:visit2 -0.834522 0.4766823 69 -1.750687  0.0844
    ##  Correlation: 
    ##                  (Intr) Gr_011 Gr_012 visit2 G_011:
    ## Group_011        -0.837                            
    ## Group_012        -0.771  0.645                     
    ## visit2           -0.684  0.572  0.527              
    ## Group_011:visit2  0.572 -0.684 -0.441 -0.837       
    ## Group_012:visit2  0.527 -0.441 -0.684 -0.771  0.645
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.67633889 -0.23364300  0.04367702  0.29648852  0.77560160 
    ## 
    ## Number of Observations: 144
    ## Number of Groups: 72

# Plot models to compare LQMM and LMM

``` r
# model using only the first timepoint (so no patients are in exacerbation)
fit.model1.02 <- lqmm( density.log1000 ~  Group_01 * visit, random = ~visit, 
                    # group = individual, because two time points
                    group = individual, tau = 0.2, covariance = "pdSymm",
                    data = IBD.data.Family, na.action = na.omit,
                    control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
fit.model1.02
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.2, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.2 
    ## 
    ## Fixed effects:
    ##      (Intercept)         Group_011         Group_012            visit2  
    ##           5.4325            0.3503            0.1329            0.8955  
    ## Group_011:visit2  Group_012:visit2  
    ##          -1.5518           -0.6628  
    ## 
    ## Covariance matrix of the random effects:
    ##             (Intercept) visit2 
    ## (Intercept)  0.7681     -0.2308
    ## visit2      -0.2308      0.5122
    ## 
    ## Residual scale parameter: 0.1688 (standard deviation 0.87)
    ## Log-likelihood: -223.2 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

``` r
sum.fit.model1.02 <- summary( fit.model1.02 )
sum.fit.model1.02
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.2, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.2 
    ## 
    ## Fixed effects:
    ##                     Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       5.43252    0.22184     4.98671      5.8783 < 2.2e-16 ***
    ## Group_011         0.35032    0.25314    -0.15839      0.8590    0.1727    
    ## Group_012         0.13290    0.43345    -0.73815      1.0040    0.7604    
    ## visit2            0.89554    0.20882     0.47591      1.3152 8.414e-05 ***
    ## Group_011:visit2 -1.55178    0.29144    -2.13745     -0.9661 2.519e-06 ***
    ## Group_012:visit2 -0.66278    0.45866    -1.58449      0.2589    0.1548    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## AIC:
    ## [1] 466.5 (df = 10)

``` r
# model using only the first timepoint (so no patients are in exacerbation)
fit.model1.05 <- lqmm( density.log1000 ~  Group_01 * visit, random = ~visit, 
                    # group = individual, because two time points
                    group = individual, tau = 0.5, covariance = "pdSymm",
                    data = IBD.data.Family, na.action = na.omit,
                    control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
fit.model1.05
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.5, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.5 
    ## 
    ## Fixed effects:
    ##      (Intercept)         Group_011         Group_012            visit2  
    ##           6.0360            0.6004            0.4393            0.8291  
    ## Group_011:visit2  Group_012:visit2  
    ##          -1.3430           -0.7513  
    ## 
    ## Covariance matrix of the random effects:
    ##             (Intercept) visit2   
    ## (Intercept)  0.119606   -0.008205
    ## visit2      -0.008205    0.120040
    ## 
    ## Residual scale parameter: 0.3172 (standard deviation 0.8973)
    ## Log-likelihood: -197.9 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

``` r
sum.fit.model1.05 <- summary( fit.model1.05 )
sum.fit.model1.05
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.5, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.5 
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.035961   0.183814    5.666573      6.4053 < 2.2e-16 ***
    ## Group_011         0.600442   0.266459    0.064973      1.1359  0.028742 *  
    ## Group_012         0.439298   0.302475   -0.168549      1.0471  0.152780    
    ## visit2            0.829122   0.148326    0.531049      1.1272 9.961e-07 ***
    ## Group_011:visit2 -1.343015   0.311212   -1.968419     -0.7176 7.709e-05 ***
    ## Group_012:visit2 -0.751264   0.245942   -1.245502     -0.2570  0.003639 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## AIC:
    ## [1] 415.7 (df = 10)

``` r
# model using only the first timepoint (so no patients are in exacerbation)
fit.model1.08 <- lqmm( density.log1000 ~  Group_01 * visit, random = ~visit, 
                    # group = individual, because two time points
                    group = individual, tau = 0.8, covariance = "pdSymm",
                    data = IBD.data.Family, na.action = na.omit,
                    control = list( verbose = FALSE, max_iter = 1000, LP_tol_ll = 1e-04 ))
fit.model1.08
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.8, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.8 
    ## 
    ## Fixed effects:
    ##      (Intercept)         Group_011         Group_012            visit2  
    ##           6.5084            0.6343            0.6158            0.7382  
    ## Group_011:visit2  Group_012:visit2  
    ##          -1.0699           -0.6744  
    ## 
    ## Covariance matrix of the random effects:
    ##             (Intercept) visit2  
    ## (Intercept)  0.14345    -0.05172
    ## visit2      -0.05172     0.07617
    ## 
    ## Residual scale parameter: 0.1877 (standard deviation 0.9675)
    ## Log-likelihood: -189.7 
    ## 
    ## Number of observations: 144 
    ## Number of groups: 72

``` r
sum.fit.model1.08 <- summary( fit.model1.08 )
sum.fit.model1.08
```

    ## Call: lqmm(fixed = density.log1000 ~ Group_01 * visit, random = ~visit, 
    ##     group = individual, covariance = "pdSymm", tau = 0.8, data = IBD.data.Family, 
    ##     na.action = na.omit, control = list(verbose = FALSE, max_iter = 1000, 
    ##         LP_tol_ll = 1e-04))
    ## 
    ## Quantile 0.8 
    ## 
    ## Fixed effects:
    ##                      Value Std. Error lower bound upper bound  Pr(>|t|)    
    ## (Intercept)       6.508390   0.198976    6.108533      6.9082 < 2.2e-16 ***
    ## Group_011         0.634283   0.224031    0.184077      1.0845  0.006708 ** 
    ## Group_012         0.615761   0.272635    0.067882      1.1636  0.028397 *  
    ## visit2            0.738242   0.163157    0.410367      1.0661 3.861e-05 ***
    ## Group_011:visit2 -1.069915   0.252161   -1.576651     -0.5632 9.766e-05 ***
    ## Group_012:visit2 -0.674357   0.354645   -1.387043      0.0383  0.063125 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## AIC:
    ## [1] 399.5 (df = 10)

``` r
# Model 1: Normal regression
# Model 1.02: Quantile regression, quantile 20%
# Model 1.05: Quantile regression, quantile 50%
# Model 1.08: Quantile regression, quantile 80%

plot_models(fit.model1.02, model1, fit.model1.05 ,fit.model1.08, show.p = TRUE, title = paste(fam.name), 
            vline.color = "black", line.size = 1.5, dot.size = 5, spacing = 0.7) + 
  theme_bw() +  
  theme(axis.text = element_text(size=14), legend.position="none") + 
  scale_color_manual(values=c( "darkblue", "Blue", "Red", "#6699FF" ))
```

![](README_files/figure-gfm/Plot%20to%20compare-1.png)<!-- -->
