###############################################################################
##  Simulations to test GIV and various other methods. This code simulates   ##
##  markers and individuals in a GWAS setting and then estimates heritability##
##  using polygenic scores                                                   ##
##                                                                           ##
##  See 'Genetic Instrumental Variable (GIV) regression:                     ##
##  Explaining socioeconomic and health outcomes in non-experimental data'   ##
##  for more details.                                                        ##
##                                                                           ##
##  This script particularly evaluates the ability to estimate heritabilities##
###############################################################################
cat('Simulations started at: ',strftime(Sys.time()))

### Import libraries ###
library(MASS)
library(foreach)
library(doParallel)

### Set paralell pool ###
no_cores <- detectCores() - 1

cat('\nno cores available: ',detectCores(), '\n')
cat('no cores used: ',no_cores, '\n')

cl <- makeCluster(no_cores, type="FORK") # FORK does not work on windows
registerDoParallel(cl) # Register paralell pool so GWAS utilise multiple cores

### Command Line input ###
args = commandArgs(TRUE) # take command line input
n_pop = as.integer(args[1]) * 1000 # total GWAS sample
n_markers = as.integer(args[2]) * 1000 # number of markers
h2 = as.numeric(args[3]) # heritability of the trait


### Other parameters, have not changed for this set of simulations ###
set.seed(1) # note: changed seed
name = 'heritability'
n_replic = 10000; # size of replication sample
rep = 20  # number of replications
spl = round(n_pop/2); # where to split the subsamples

### prep ###
output=sprintf('%s_N_%i_M_%i_h2_%1.2f.RData',name,n_pop,n_markers,h2)

### Define a few functions ###
# Fast OLS
fols = function (y, x) {
     XtX = crossprod(x)
     Xty = crossprod(x, y)
     solve(XtX, Xty)
}

# calculate variance-covariance matrix of OLS
olscov = function(y,X,bhat,n) {
    k = dim(bhat)[1]
    e = y - X %*% bhat
    s2 = crossprod(e)/(n-k)
    return(solve(crossprod(X)) * s2[1,1]) # first element of 1 by 1 matrix
}

# calculate variance-covariance matrix of 2SLS
ivcov = function(y,X,Xpred,bhat,n) {
    e = y - X %*% bhat
    dim(X)
    dim(bhat)
    s2 = crossprod(e)/n
    return(solve(crossprod(Xpred)) * s2[1,1]) # first element of 1 by 1 matrix
}


### Allocating Matrices for Results ###
est_IV = matrix(0L,nrow=rep,ncol=2)
est_h2 = matrix(0L,nrow=rep,ncol=1)

cov_IV=array(0L, c(2, 2, rep))
var_h2 = matrix(0L,nrow=rep,ncol=1)

ones1 = matrix(1L,nrow=spl,ncol=1)
ones2 = matrix(1L,nrow=n_pop-spl,ncol=1)
onesr = matrix(1L,nrow=n_replic,ncol=1)


for(i_rep in 1:rep){
    ### Generate Data ###
    markers = matrix(rbinom(n_markers*n_pop,2,0.5),ncol=n_markers)
        # Drawing markers from a binomial distribution, all have MAF 0.5
    zeta = mvrnorm(n_markers, mu=0, Sigma = 1) # Coefficients for the markers
    errors = mvrnorm(n_pop, mu=0, Sigma = 1) # error term in for y

    # Setting the correct weight for the markers and error terms
    m_weight = sqrt(h2)/sqrt(0.5*n_markers)
    e_weight = sqrt(1-h2)

    # Generate y
    y = markers %*% zeta * m_weight + errors * e_weight

    ### Run GWAS ###
    # Runs paralell over multiple cores, each loop is passed one column of the
    # markers matrix.
    res <- foreach(marker=iter(markers,by='col'), .combine=rbind) %dopar% {
        # Split the sample in two: #
        X1 = cbind(ones1,marker[1:spl])
        X2 = cbind(ones2,marker[(spl+1):n_pop])

        ## Normal GWAS ##
        beta1 = fols(y[1:spl],X1) # OLS on half of the sample
        beta2 = fols(y[(spl+1):n_pop],X2) # Other half

        return(c(beta1[2],beta2[2]))
    }
    # Drop some variables for memory constraints
    rm(markers,y,errors)
    gc()

    ## Generate Replication Sample ##
    markers_R = matrix(rbinom(n_markers*n_replic,2,0.5),ncol=n_markers)
    errors_R = mvrnorm(n_replic, mu=0, Sigma=1)

    y_R = markers_R %*% zeta * m_weight + errors_R * e_weight
    scores = markers_R %*% res

    # Drop some variables for memory constraints
    rm(markers_R)
    gc()

    # standardize the scores
    S1 = scores[,1]/sqrt(var(scores[,1]))
    S2 = scores[,2]/sqrt(var(scores[,2]))

    ## GIV ##
    X = cbind(onesr, S1) # regressors
    Z = cbind(onesr, S2) # instruments
    A = fols(X,Z) # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A # Predicted values
    est = fols(y_R,pred_X) # OLS (Second stage 2sls)
    est_IV[i_rep,]= est
    cov_IV[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)

    # estimate the heritiability
    est_h2[i_rep] = est_IV[i_rep,2]^2 * cor(S1,S2)

    # estimate the variance of the estimate using the delta method
    var_h2[i_rep] = cov_IV[2,2,i_rep] * (2 *  cor(S1,S2) * est_IV[i_rep,2])^2

    cat(sprintf('\nRepitition %i completed at:',i_rep),strftime(Sys.time()))
    save(list = ls(all.names = TRUE),file=paste(output));

}
cat('\n\nMean h2:',mean(est_h2))
cat('\nMean se:',mean(sqrt(var_h2)))
