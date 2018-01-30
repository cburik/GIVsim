###############################################################################
##  Simulations to test GIV and various other methods. This code simulates   ##
##  markers and individuals in a GWAS setting and then estimates the various ##
##  methods in an independent sample.                                        ##
##                                                                           ##
##  See 'Genetic Instrumental Variable (GIV) regression:                     ##
##  Explaining socioeconomic and health outcomes in non-experimental data'   ##
##  for more details.                                                        ##
##                                                                           ##
##  This script particularly evaluates a model with pleiotropy and           ##
##  genetic-unrelated endogeneity with partial control for the endogeneity   ##
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
delta = as.numeric(args[1]) # coefficient for T
rho = as.numeric(args[2]) # correlation in case of pleitropy
rho_e = as.numeric(args[3])  # correlation for endogeneity
s = as.numeric(args[4]) # share of variance of measurable endogeneity


### Other parameters, have not changed for this set of simulations ###
set.seed(1)
name = 'partial_control'

n_pop = 100000  # population size
spl = round(n_pop/2)  # where to split the subsamples
n_replic = 10000  # size of replication sample
n_markers = 10000  # number of markers
h2 = 0.5 # heritability for both traits
rep = 20  # number of replications

output=sprintf('%s_s_%1.1f_rho_%1.1f_rhoe_%1.1f_delta_%i.RData',name,s,rho,
    rho_e,delta)

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
GIV_U = matrix(0L,nrow=rep,ncol=4)
GIV_C = matrix(0L,nrow=rep,ncol=4)
OLS_T = matrix(0L,nrow=rep,ncol=3)
OLS_T_Sy =  matrix(0L,nrow=rep,ncol=4)
OLS_T_Sy_cond= matrix(0L,nrow=rep,ncol=4)
MR = matrix(0L,nrow=rep,ncol=3)
EMR = matrix(0L,nrow=rep,ncol=4)
EMR2 = matrix(0L,nrow=rep,ncol=4)

cov_GIV_U = array(0L, c(4, 4, rep))
cov_GIV_C = array(0L, c(4, 4, rep))
cov_OLS_T = array(0L, c(3, 3, rep))
cov_OLS_T_Sy = array(0L, c(4, 4, rep))
cov_OLS_T_Sy_cond = array(0L, c(4, 4, rep))
cov_MR =array(0L, c(3, 3, rep))
cov_EMR = array(0L, c(4, 4, rep))
cov_EMR2 = array(0L, c(4, 4, rep))

### Setting parameters to draw markers and error terms ###
mub = c(0,0)
sigmab = matrix(c(1,rho,rho,1), ncol=2)

mue = c(0,0,0)

# calculalculating variance for the two parts of the error terms
e1w = sqrt(s) # weight of error part that will get controlled for
e2w = sqrt(1-s) # other part
cor_tmp = rho_e * sqrt(e1w^2 + e2w^2)/(e1w + e2w)
sigmae = matrix(c(1,0,cor_tmp,
                0,1,cor_tmp,
                cor_tmp,cor_tmp,1), ncol=3) # 3x3 (2 for y, 1 for T)

ones1 = matrix(1L,nrow=spl,ncol=1)
ones2 = matrix(1L,nrow=n_pop-spl,ncol=1)
onesr = matrix(1L,nrow=n_replic,ncol=1)

for(i_rep in 1:rep){
    ### Generate Data ###
    markers = matrix(rbinom(n_markers*n_pop,2,0.5),ncol=n_markers)
            # Drawing markers from a binomial distribution, all have MAF 0.5
    zeta = mvrnorm(n_markers, mu=mub, Sigma = sigmab) # Coefficients for the
                                        # markers, one column of y and one for T
    errors = mvrnorm(n_pop, mu=mue, Sigma = sigmae) # two columns of errors for y
                                                   # and one for T

   # Setting the correct weight for the markers and error terms for T
   m_weight = sqrt(h2)/sqrt(0.5*n_markers) #ignoring pleiotropy and endogeneity
   e_weight = sqrt(1-h2)

   # Generating y and T
    T = markers %*% zeta[,2] * m_weight  + errors[,3] * e_weight
    y =  markers %*% zeta[,1]  * m_weight + delta*T + (e1w * errors[,1] + e2w * errors[,2]) * e_weight

    ### Run GWAS ###
    # Runs paralell over multiple cores, each loop is passed one column of the
    # markers matrix.
    res <- foreach(marker=iter(markers,by='col'), .combine=rbind) %dopar% {

        # Split the sample in two: #
        X1 = cbind(ones1,marker[1:spl]) # Regressor matrix for first half
        X2 = cbind(ones2,marker[(spl+1):n_pop]) # Regressor matrix for second

        ## Normal GWAS ##
        beta1_y = fols(y[1:spl],X1)  # OLS on half of the sample
        beta2_y = fols(y[(spl+1):n_pop],X2)  # Other half
        beta_y_full = fols(y,rbind(X1,X2)) # OLS on the full sample
        beta_T = fols(T[(spl+1):n_pop],X2) # GWAS for T on the second sample

        ## Conditional GWAS ##
        beta1_y_cond =  fols(y[1:spl],cbind(X1,T[1:spl]))  # OLS conditioning on T
        beta2_y_cond =  fols(y[(spl+1):n_pop],cbind(X2,T[(spl+1):n_pop]))
        beta_y_cond_full = fols(y,cbind(rbind(X1,X2),T))

        # only return the coefficients for the markers, not the constant
        return(c(beta1_y[2],beta2_y[2],beta_T[2],beta1_y_cond[2],
            beta2_y_cond[2], beta_y_full[2], beta_y_cond_full[2]))
    }

    # Drop some variables for memory constraints
    rm(markers,T,y,errors)
    gc()

    ## Generate Replication Sample ##
    markers_R = matrix(rbinom(n_markers*n_replic,2,0.5),ncol=n_markers)
    errors_R = mvrnorm(n_replic, mu=mue, Sigma = sigmae)
    T_R = markers_R %*% zeta[,2] * m_weight + errors_R[,3]* e_weight
    y_R = T_R*delta + markers_R %*% zeta[,1] * m_weight + (e1w * errors_R[,1] + e2w * errors_R[,2]) * e_weight

    ## create scores ##
    scores = markers_R %*% res[,1:3]  # N by 3, two scores for y, one for T
    cond_scores = markers_R %*% res[,4:5]  # N by 2, two scores for y
    full_score = markers_R %*% res[,6]  # one score for y using the full sample
    cond_full_score = markers_R %*% res[,7] # one score for y using the full sample

    # Drop some variables for memory constraints
    rm(markers_R)
    gc()

    ### The different estimation methods ###
    # note: here we plug in the partial control
    ## OLS on only T ##
    X = cbind(onesr, T_R, errors_R[,1])   # regressors
    est = fols(y_R,X)  # OLS
    OLS_T[i_rep,]= est
    cov_OLS_T[,,i_rep] = olscov(y_R,X,est,n_replic)

    ## OLS on T and S(y) ##
    X = cbind(onesr, full_score, T_R, errors_R[,1])   # regressors
    est = fols(y_R,X)  # OLS
    OLS_T_Sy[i_rep,]= est
    cov_OLS_T_Sy[,,i_rep] = olscov(y_R,X,est,n_replic)

    ## OLS on T and S(y|T) ##
    X = cbind(onesr, cond_full_score, T_R, errors_R[,1])   # regressors
    est = fols(y_R,X)  # OLS
    OLS_T_Sy_cond[i_rep,]= est
    cov_OLS_T_Sy_cond[,,i_rep] = olscov(y_R,X,est,n_replic)

    ## MR ##
    X = cbind(onesr, T_R, errors_R[,1])   # regressors
    Z = cbind(onesr, scores[,3], errors_R[,1])  # instruments
    A = fols(X,Z)  # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A  # Predicted values
    est = fols(y_R,pred_X)  # OLS (Second stage 2sls)
    MR[i_rep,] = est
    cov_MR[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)

    ## GIV conditional ##
    X = cbind(onesr, cond_scores[,1], T_R, errors_R[,1])   # regressors
    Z = cbind(onesr, scores[,2],  T_R, errors_R[,1])   # instruments are the S(y2) and T itself.
    A = fols(X,Z)  # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A  # Predicted values
    est = fols(y_R,pred_X)  # OLS (Second stage 2sls)
    GIV_C[i_rep,] = est
    cov_GIV_C[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)

    ## GIV unconditional ##
    X = cbind(onesr, scores[,1], T_R, errors_R[,1])   # regressors
    Z = cbind(onesr, scores[,2],  T_R, errors_R[,1])   # instruments are the S(y2) and T itself.
    A = fols(X,Z)  # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A  # Predicted values
    est = fols(y_R,pred_X)  # OLS (Second stage 2sls)
    GIV_U[i_rep,] = est # OLS (Second stage 2sls)
    cov_GIV_U[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)

    ## EMR ##
    X = cbind(onesr, scores[,1], T_R, errors_R[,1])   # regressors
    Z = cbind(onesr, scores[,1], scores[,3], errors_R[,1])  # instruments
    A = fols(X,Z)  # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A  # Predicted values
    est = fols(y_R,pred_X)  # OLS (Second stage 2sls)
    EMR[i_rep,] = est
    cov_EMR[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)

    ## EMR2 ##
    X = cbind(onesr, cond_scores[,1], T_R, errors_R[,1])   # regressors
    Z = cbind(onesr, scores[,2:3], errors_R[,1])   # instruments
    A = fols(X,Z)  # Multivariate OLS (first stage 2sls)
    pred_X = Z %*% A  # Predicted values
    est = fols(y_R,pred_X)  # OLS (Second stage 2sls)
    EMR2[i_rep,]= est
    cov_EMR2[,,i_rep] = ivcov(y_R,X,pred_X,est,n_replic)


    cat(sprintf('\nRepitition %i completed at:',i_rep),strftime(Sys.time()))
    save(list = ls(all.names = TRUE),file=paste(output)) # save all variables
                                                        # including the results
}
