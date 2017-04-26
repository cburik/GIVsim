%% Simulations for Genetic IV regression
%   Simulation of all 7 methods described in the main text
%
%   Data Generating Model: 
%       y = a + d T + g sy0 + e
%       T = sT0 + u
%
%       Sy0 and ST0 are correlated, e and u correlated
%   
%   Casper Burik
%      06-Mar-2017

clear;
mkdir('Simulation Output')
outputname = 'Simulation Output/standard.endogenous.mat';

%% Model specifics
h2y = 0.2;              % heritability y
h2T = 0.55;             % heritability T

rho   = 0.15;           % correlation sT0 and sy0
rhoer = 0.4;            % correlation between error in T and y for additional endogeneity

delta = 0.15;           % effect size T on y

sp  = 0.5;              % sample split of GWAS for IV 
sp2 = [1/3,1/3,1/3];    % for when it's split three ways            

my  = 300000;           % number of SNPs y
mT  = 300000;           % number of SNPs T

ngy = 200000:50000:2000000; % total gwas sample size for y
ngT = 200000:50000:2000000; % total gwas sample size for T
    % at the moment sizes need to be equal! 
    
nrs = 8600;             % Replication sample size (Approximately HRS sample size)
s   = 1;             % Number of repitions

%% Caculate Model Parameters 
beta  = sqrt(h2T);
gamma = -delta*beta*rho + sqrt(delta^2*beta^2*rho^2 - delta^2*beta^2+h2y) ;
v_eta = 1 - h2T;
v_eps = (-delta*rhoer*sqrt(v_eta) + sqrt(delta^2*rhoer^2*v_eta - delta^2*v_eta + 1 - h2y))^2; 

%% Allocate Matrix Sizes
[~,h]   = size(ngy);

ghat_1  = zeros(s,h);
ghat_2  = zeros(s,h);
ghat_3  = zeros(s,h);
ghat_4  = zeros(s,h);
ghat_5  = zeros(s,h);
ghat_6  = zeros(s,h);
ghat_7  = zeros(s,h);

dhat_1  = zeros(s,h);
dhat_2  = zeros(s,h);
dhat_3  = zeros(s,h);
dhat_4  = zeros(s,h);
dhat_5  = zeros(s,h);
dhat_6  = zeros(s,h);
dhat_7  = zeros(s,h);

%% Simulation loops
 
for ih = 1:h % loop over GWAS sample sizes
    % variance of measurement errors
    v_ey1    = 1/h2y * (my - h2y)/ ngy(ih);          % full GWAS sample
    v_ey2a   = 1/h2y * (my - h2y)/ (ngy(ih)*sp);     % GWAS sample split in two
    v_ey2b   = 1/h2y * (my - h2y)/ (ngy(ih)*(1-sp));
    v_ey3a   = 1/h2y * (my - h2y)/ (ngy(ih)*sp2(1)); % GWAS sample split three ways
    v_ey3b   = 1/h2y * (my - h2y)/ (ngy(ih)*sp2(2));
    v_ey3c   = 1/h2y * (my - h2y)/ (ngy(ih)*sp2(3));
    v_eT     = 1/h2T * (mT - h2T)/ (ngT(ih));
    
    for is = 1:s % loop over simulated samples
        % generate data
        % true scores
        mu = [0,0]; 
        sigma = [1,rho;rho,1];

        scores = mvnrnd(mu,sigma,nrs); % True scores are draw from multi variate normal
        sT0    = scores(:,1);
        sy0    = scores(:,2);
        
        % measured scores
        % measurement errors
        ey1     = sqrt(v_ey1) * randn(nrs,1);
        ey2a    = sqrt(v_ey2a) * randn(nrs,1);
        ey2b    = sqrt(v_ey2b) * randn(nrs,1);
        ey3a    = sqrt(v_ey3a) * randn(nrs,1);
        ey3b    = sqrt(v_ey3b) * randn(nrs,1);
        ey3c    = sqrt(v_ey3c) * randn(nrs,1);
        eT      = sqrt(v_eT) * randn(nrs,1);
        
        % measured scores
        sy1     = sy0 + ey1;
        sy2a    = sy0 + ey2a;
        sy2b    = sy0 + ey2b;
        sy3a    = sy0 + ey3a;
        sy3b    = sy0 + ey3b;
        sy3c    = sy0 + ey3c;
        sT      = sT0 + eT;
        
        % T and y
        muer  = [0,0];
        cover = rhoer * sqrt(v_eta * v_eps);
        siger = [v_eta, cover;cover, v_eps];
        errors= mvnrnd(muer,siger,nrs); % draw from multivariate normal

        T = beta * sT0 + errors(:,1);
        y = delta * T + gamma * sy0 + errors(:,2);
        
        % calculations
        % OLS (1)
        xmat = [ones(nrs,1), T,sy1]; % matrix of regressors
        bhat = xmat\y;               % OLS
        dhat_1(is,ih) = bhat(2);     % Save coefficients
        ghat_1(is,ih) = bhat(3);
        
        % Only MR (2)
        xmat = [ones(nrs,1), T];  % matrix of regressors
        zmat = [ones(nrs,1), sT]; % matrix of instruments
        A    = zmat\xmat;         % First stage of 2SLS
        xfit = zmat*A;            % Fitted values
        bhat = xfit\y;            % Second stage of 2SLS
        
        dhat_2(is,ih) = bhat(2);  % Save coefficients
        ghat_2(is,ih) = 0;        % Sy1 is left out of regression
                
                
        % EMR with sT (3)
        xmat = [ones(nrs,1), T, sy1];  % matrix of regressors
        zmat = [ones(nrs,1), sT, sy1]; % matrix of instruments
        A    = zmat\xmat;              % First stage of 2SLS
        xfit = zmat*A;                 % Fitted values
        bhat = xfit\y;                 % Second stage of 2SLS
                
        dhat_3(is,ih) = bhat(2);  % Save coefficients
        ghat_3(is,ih) = bhat(3);
        
        % EMR with sy2 (4)
        xmat = [ones(nrs,1), T, sy2a];    % matrix of regressors
        zmat = [ones(nrs,1), sy2b, sy2a]; % matrix of instruments
        A    = zmat\xmat;                 % First stage of 2SLS
        xfit = zmat*A;                    % Fitted values
        bhat = xfit\y;                    % Second stage of 2SLS
                
        dhat_4(is,ih) = bhat(2);  % Save coefficients
        ghat_4(is,ih) = bhat(3);
        
               
        % GIV with sy2 and sy3 as instruments (5)
        xmat = [ones(nrs,1), T, sy3a];     % matrix of regressors
        zmat = [ones(nrs,1), sy3b, sy3c];  % matrix of instruments
        A    = zmat\xmat;                  % First stage of 2SLS
        xfit = zmat*A;                     % Fitted values
        bhat = xfit\y;                     % Second stage of 2SLS
        
        dhat_5(is,ih) = bhat(2);  % Save coefficients
        ghat_5(is,ih) = bhat(3);
              

        % GIV with sT and sy2 as instruments (6)
        xmat = [ones(nrs,1), T, sy2a];    % matrix of regressors
        zmat = [ones(nrs,1), sT, sy2b];   % matrix of instruments
        A    = zmat\xmat;                 % First stage of 2SLS
        xfit = zmat*A;                    % Fitted values
        bhat = xfit\y;                    % Second stage of 2SLS
        
        dhat_6(is,ih) = bhat(2);  % Save coefficients
        ghat_6(is,ih) = bhat(3);
                
        % GIV with sT, sy2 and sy3 as instruments (7)
        xmat = [ones(nrs,1), T, sy3a];         % matrix of regressors
        zmat = [ones(nrs,1), sT, sy3b, sy3c];  % matrix of instruments
        A    = zmat\xmat;                      % First stage of 2SLS
        xfit = zmat*A;                         % Fitted values
        bhat = xfit\y;                         % Second stage of 2SLS
        
        dhat_7(is,ih) = bhat(2);  % Save coefficients
        ghat_7(is,ih) = bhat(3);
        
    end
end
 
save(outputname);
