## A stan model using mono, bi and trispecific time series to fit demographic parameters on Tet, Col and Ble.

rm(list=ls())

# load packages ##############################################################
library(rstan)
library(deSolve)
library(coda)

dir.create("out", showWarnings = FALSE) ## create a directory to store outputs

## Declare Stan model ########################################################

stanmodelcode = 'functions{
real[] ode_monospecific(real t, real[] N, real[] p, real[] x_r, int[] x_i){  // function to simulate monospecific time series
// p[1]=r1, p[2]=r2, p[3]=r3
// p[4]=a11, p[5]=a21, p[6]=a31
// p[7]=a12, p[8]=a22, p[9]=a32
// p[10]=a13, p[11]=a23, p[12]=a33
// aij is the effect of species i on j
real dNdt[3];
if(N[1]<1e-6) dNdt[1]=0.0; else dNdt[1] = (p[1] - p[4] *N[1]) *N[1];
if(N[2]<1e-6) dNdt[2]=0.0; else dNdt[2] = (p[2] - p[8] *N[2]) *N[2]; 
if(N[3]<1e-6) dNdt[3]=0.0; else dNdt[3] = (p[3] - p[12]*N[3]) *N[3]; 
return dNdt;
}

real[] ode_Ble_Col(real t, real[] N, real[] p, real[] x_r, int[] x_i){  // simulates bispecific (Ble + Col) time series
// p[1]=r1, p[2]=r2, p[3]=r3
// p[4]=a11, p[5]=a21, p[6]=a31
// p[7]=a12, p[8]=a22, p[9]=a32
// p[10]=a13, p[11]=a23, p[12]=a33
// aij is the effect of species i on j
real dNdt[2];
if(N[1]<1e-6) dNdt[1]=0.0; else dNdt[1] = (p[1] - p[4] *N[1] - p[5] *N[2]) *N[1];
if(N[2]<1e-6) dNdt[2]=0.0; else dNdt[2] = (p[2] - p[7] *N[1] - p[8] *N[2]) *N[2];
return dNdt;
}

real[] ode_Ble_Tet(real t, real[] N, real[] p, real[] x_r, int[] x_i){  // simulates bispecific (Ble + Tet) time series
// p[1]=r1, p[2]=r2, p[3]=r3
// p[4]=a11, p[5]=a21, p[6]=a31
// p[7]=a12, p[8]=a22, p[9]=a32
// p[10]=a13, p[11]=a23, p[12]=a33
// aij is the effect of species i on j
real dNdt[2];
if(N[1]<1e-6) dNdt[1]=0.0; else dNdt[1] = (p[1] - p[4] *N[1] - p[6] *N[2]) *N[1];
if(N[2]<1e-6) dNdt[2]=0.0; else dNdt[2] = (p[3] - p[10]*N[1] - p[12]*N[2]) *N[2];
return dNdt;
}

real[] ode_trispecific(real t, real[] N, real[] p, real[] x_r, int[] x_i){  // simulates trispecific time series
// p[1]=r1, p[2]=r2, p[3]=r3
// p[4]=a11, p[5]=a21, p[6]=a31
// p[7]=a12, p[8]=a22, p[9]=a32
// p[10]=a13, p[11]=a23, p[12]=a33
// aij is the effect of species i on j
real dNdt[3];
if(N[1]<1e-6) dNdt[1]=0.0; else dNdt[1] = (p[1] - p[4] *N[1] - p[5] *N[2] - p[6] *N[3]) *N[1];
if(N[2]<1e-6) dNdt[2]=0.0; else dNdt[2] = (p[2] - p[7] *N[1] - p[8] *N[2] - p[9] *N[3]) *N[2];
if(N[3]<1e-6) dNdt[3]=0.0; else dNdt[3] = (p[3] - p[10]*N[1] - p[11]*N[2] - p[12]*N[3]) *N[3];
return dNdt;
}
}



data{
int n; // observations
int m; // replicates

real t[n]; //vector of times

int Ble[m,n]; // monospecific densities
int Col[m,n];
int Tet[m,n];


int Ble_w_col[m,n]; // bispecific of Ble (with col)
int Ble_w_tet[m,n]; // bispecific of Ble (with tet)

int Col_w_ble[m,n]; // bispecific of Col (with ble)
int Tet_w_ble[m,n]; // bispecific of Tet (with ble)


int Ble_tri[m,n]; // trispecific densities
int Col_tri[m,n];
int Tet_tri[m,n];
}


parameters{
real<lower=0, upper=3> r[3]; // growth rates (1 Ble, 2 Col, 3 Tet)


real<lower=0> a[3,3]; // competition rates


real<lower=0> N0_sim_Ble[m]; // individual initial values (monospecific)
real<lower=0> N0_sim_Col[m];
real<lower=0> N0_sim_Tet[m];

real<lower=0> N0_sim_Ble_w_col[m]; // individual initial values (bispecific)
real<lower=0> N0_sim_Ble_w_tet[m];
real<lower=0> N0_sim_Col_w_ble[m];
real<lower=0> N0_sim_Tet_w_ble[m];

real<lower=0> N0_sim_Ble_tri[m]; // individual initial values (trispecific)
real<lower=0> N0_sim_Col_tri[m];
real<lower=0> N0_sim_Tet_tri[m];


real<lower=0> sdev[3];
}


model{
real p[12]; // parameters for integrator

real Nsim_mono[n-1,3]; // simulated values, monospecific cultures
real Nsim_Ble_Col[n-1,2]; // bispecific
real Nsim_Ble_Tet[n-1,2];
real Nsim_tri[n-1,3]; // trispecific

//priors
r[1] ~ lognormal(-2,1);
r[2] ~ lognormal(-2,1);
r[3] ~ lognormal(-2,1);

a[1,1] ~ gamma(2,1);
a[1,2] ~ gamma(2,1);
a[1,3] ~ gamma(2,1);

a[2,1] ~ gamma(2,1);
a[2,2] ~ gamma(2,1);
a[2,3] ~ gamma(2,1);

a[3,1] ~ gamma(2,1);
a[3,2] ~ gamma(2,1);
a[3,3] ~ gamma(2,1);

sdev ~ gamma(2,1);

N0_sim_Ble ~ normal(0,10);
N0_sim_Col ~ normal(0,100);
N0_sim_Tet ~ normal(0,1000);

N0_sim_Ble_w_col ~ normal(0,10);
N0_sim_Ble_w_tet ~ normal(0,10);
N0_sim_Col_w_ble ~ normal(0,100);
N0_sim_Tet_w_ble ~ normal(0,1000);

N0_sim_Ble_tri ~ normal(0,10);
N0_sim_Col_tri ~ normal(0,100);
N0_sim_Tet_tri ~ normal(0,1000);

// parameters for integrator
p[1] = r[1];
p[2] = r[2];
p[3] = r[3];

p[4] = a[1,1];
p[5] = a[2,1];
p[6] = a[3,1];

p[7] = a[1,2];
p[8] = a[2,2];
p[9] = a[3,2];

p[10] = a[1,3];
p[11] = a[2,3];
p[12] = a[3,3];


for (j in 1:m){
// integrate ODE
Nsim_mono = integrate_ode_rk45(ode_monospecific,{N0_sim_Ble[j],N0_sim_Col[j],N0_sim_Tet[j]},t[1],t[2:n],p,rep_array(0.0,0),rep_array(0,0));

Nsim_Ble_Col = integrate_ode_rk45(ode_Ble_Col,{N0_sim_Ble_w_col[j],N0_sim_Col_w_ble[j]},t[1],t[2:n],p,rep_array(0.0,0),rep_array(0,0));

Nsim_Ble_Tet = integrate_ode_rk45(ode_Ble_Tet,{N0_sim_Ble_w_tet[j],N0_sim_Tet_w_ble[j]},t[1],t[2:n],p,rep_array(0.0,0),rep_array(0,0));

Nsim_tri = integrate_ode_rk45(ode_trispecific,{N0_sim_Ble_tri[j],N0_sim_Col_tri[j],N0_sim_Tet_tri[j]},t[1],t[2:n],p,rep_array(0.0,0),rep_array(0,0));


// likelihood
Ble[j,1] ~ neg_binomial_2(N0_sim_Ble[j],1/sdev[1]);
Col[j,1] ~ neg_binomial_2(N0_sim_Col[j],1/sdev[2]);
Tet[j,1] ~ neg_binomial_2(N0_sim_Tet[j],1/sdev[3]);

Ble_w_col[j,1] ~ neg_binomial_2(N0_sim_Ble_w_col[j],1/sdev[1]);
Ble_w_tet[j,1] ~ neg_binomial_2(N0_sim_Ble_w_tet[j],1/sdev[1]);
Col_w_ble[j,1] ~ neg_binomial_2(N0_sim_Col_w_ble[j],1/sdev[2]);
Tet_w_ble[j,1] ~ neg_binomial_2(N0_sim_Tet_w_ble[j],1/sdev[3]);

Ble_tri[j,1] ~ neg_binomial_2(N0_sim_Ble_tri[j],1/sdev[1]);
Col_tri[j,1] ~ neg_binomial_2(N0_sim_Col_tri[j],1/sdev[2]);
Tet_tri[j,1] ~ neg_binomial_2(N0_sim_Tet_tri[j],1/sdev[3]);

for(i in 2:n){
Ble[j,i] ~ neg_binomial_2(Nsim_mono[i-1,1],1/sdev[1]);
Col[j,i] ~ neg_binomial_2(Nsim_mono[i-1,2],1/sdev[2]);
Tet[j,i] ~ neg_binomial_2(Nsim_mono[i-1,3],1/sdev[3]);

Ble_w_col[j,i] ~ neg_binomial_2(Nsim_Ble_Col[i-1,1],1/sdev[1]);
Ble_w_tet[j,i] ~ neg_binomial_2(Nsim_Ble_Tet[i-1,1],1/sdev[1]);
Col_w_ble[j,i] ~ neg_binomial_2(Nsim_Ble_Col[i-1,2],1/sdev[2]);
Tet_w_ble[j,i] ~ neg_binomial_2(Nsim_Ble_Tet[i-1,2],1/sdev[3]);

Ble_tri[j,i] ~ neg_binomial_2(Nsim_tri[i-1,1],1/sdev[1]);
Col_tri[j,i] ~ neg_binomial_2(Nsim_tri[i-1,2],1/sdev[2]);
Tet_tri[j,i] ~ neg_binomial_2(Nsim_tri[i-1,3],1/sdev[3]);
}
}
}
'

s_model = stan_model(model_code=stanmodelcode) ## compiling the model


density_data = read.table("./Data/density_data.csv", header = T, sep = ',') ## reading data

sampling_times = unique(density_data$hours)

n_reps = 3
Ble = data.frame(matrix(nrow = 3, ncol = 13))
Col = data.frame(matrix(nrow = 3, ncol = 13))
Tet = data.frame(matrix(nrow = 3, ncol = 13))

Ble_w_col = data.frame(matrix(nrow = 3, ncol = 13))
Ble_w_tet = data.frame(matrix(nrow = 3, ncol = 13))
Col_w_ble = data.frame(matrix(nrow = 3, ncol = 13))
Tet_w_ble = data.frame(matrix(nrow = 3, ncol = 13))

Ble_tri = data.frame(matrix(nrow = 3, ncol = 13))
Col_tri = data.frame(matrix(nrow = 3, ncol = 13))
Tet_tri = data.frame(matrix(nrow = 3, ncol = 13))

for (rep in 1:n_reps){
  Ble[rep,] = round(density_data$Ble[density_data$replicate == rep & density_data$community == "Ble"])
  Col[rep,] = round(density_data$Col[density_data$replicate == rep & density_data$community == "Col"])
  Tet[rep,] = round(density_data$Tet[density_data$replicate == rep & density_data$community == "Tet"])
  
  Ble_w_col[rep,] = round(density_data$Ble[density_data$replicate == rep & density_data$community == "Ble_Col"])
  Ble_w_tet[rep,] = round(density_data$Ble[density_data$replicate == rep & density_data$community == "Ble_Tet"])
  Col_w_ble[rep,] = round(density_data$Col[density_data$replicate == rep & density_data$community == "Ble_Col"])
  Tet_w_ble[rep,] = round(density_data$Tet[density_data$replicate == rep & density_data$community == "Ble_Tet"])
  
  Ble_tri[rep,] = round(density_data$Ble[density_data$replicate == rep & density_data$community == "Tet_Ble_Col"])
  Col_tri[rep,] = round(density_data$Col[density_data$replicate == rep & density_data$community == "Tet_Ble_Col"])
  Tet_tri[rep,] = round(density_data$Tet[density_data$replicate == rep & density_data$community == "Tet_Ble_Col"])
}



data = list(n  = length(sampling_times),
            m  = nrow(Ble),
            
            Ble = Ble,
            Tet = Tet,
            Col = Col,
            
            Ble_w_col = Ble_w_col,
            Ble_w_tet = Ble_w_tet,
            Col_w_ble = Col_w_ble,
            Tet_w_ble = Tet_w_ble,
            
            Ble_tri = Ble_tri,
            Col_tri = Col_tri,
            Tet_tri = Tet_tri,
            
            t  = sampling_times)

# stan options
chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = chains)

iter   =  10000
warmup =  2000
thin   =     1

# initial values for sampling 
init=rep(list(list(r=c(0.03, 0.05, 0.1),
                   a=matrix(c(0.01,0.001,0.001, 0.0001,0.001,0.0001, 0.0001,0.0001,0.001 ),3,3,byrow=TRUE),
                   N0_sim_Ble=rep(2,data$m),
                   N0_sim_Col=rep(90,data$m),
                   N0_sim_Tet=rep(500,data$m),
                   
                   N0_sim_Ble_w_col=rep(2,data$m),
                   N0_sim_Ble_w_tet=rep(2,data$m),
                   N0_sim_Col_w_ble=rep(90,data$m),
                   N0_sim_Tet_w_ble=rep(500,data$m),
                   
                   N0_sim_Ble_tri=rep(2,data$m),
                   N0_sim_Col_tri=rep(90,data$m),
                   N0_sim_Tet_tri=rep(500,data$m),
                   sdev=c(0.1,0.1,0.1)))
         ,chains)

# run model and print result
fit_obs = sampling(s_model,
                   data=data,
                   iter=iter,
                   warmup=warmup,
                  thin=thin,
                  chains=chains,
                  init=init,
                  control = list(adapt_delta = 0.9, max_treedepth=12),
                  refresh=10
)
 
save(fit_obs, file="./out/3_species_fit_posterior.RData")

#load("./out/3_species_fit_posterior.RData")

print(fit_obs)

samples=As.mcmc.list(fit_obs)


pdf("./out/3_species_fit_chains.pdf")
plot(samples)
dev.off()

# pairs(fit_obs, pars=c("r[1]","r[2]", "r[3]",
#                       "a[1,1]","a[1,2]", "a[1,3]",
#                       "a[2,1]","a[2,2]", "a[2,3]",
#                       "a[3,1]","a[3,2]", "a[3,3]"))
# 
# check_hmc_diagnostics(fit_obs)
# 
# results_obs = summary(fit_obs)$summary[c("r[1]","r[2]", "r[3]",
#                                          "a[1,1]","a[1,2]", "a[1,3]",
#                                          "a[2,1]","a[2,2]", "a[2,3]",
#                                          "a[3,1]","a[3,2]", "a[3,3]",
#                                          "sdev[1]","sdev[2]","sdev[3]"),]
# 
# results_obs = cbind(results_obs, "divergent"=rep(get_num_divergent(fit_obs),nrow(results_obs)))
# 
# save(results_obs, file = "./out/3_species_fit.RData")


## declaring function to check the fit ###############################################
ode.model.mono = function(t,N,p){
  with(as.list(p),{
    dNdt = c(0, 0, 0)
    dNdt[1] = N[1] * (p[1] - p[4]*N[1])
    dNdt[2] = N[2] * (p[2] - p[8]*N[2])
    dNdt[3] = N[3] * (p[3] - p[12]*N[3])
    return(list(dNdt))
  })
}

ode.model.ble.col = function(t,N,p){
  with(as.list(p),{
    dNdt = c(0, 0)
    dNdt[1] = N[1] * (p[1] - p[4]*N[1] - p[5]*N[2])
    dNdt[2] = N[2] * (p[2] - p[7]*N[1] - p[8]*N[2])
    return(list(dNdt))
  })
}

ode.model.ble.tet = function(t,N,p){
  with(as.list(p),{
    dNdt = c(0, 0)
    dNdt[1] = N[1] * (p[1] - p[4]*N[1]  - p[6]*N[2])
    dNdt[2] = N[2] * (p[3] - p[10]*N[1] - p[12]*N[2])
    return(list(dNdt))
  })
}

ode.model.tri = function(t,N,p){
  with(as.list(p),{
    dNdt = c(0, 0, 0)
    dNdt[1] = N[1] * (p[1] - p[4]*N[1] - p[5]*N[2] - p[6]*N[3])
    dNdt[2] = N[2] * (p[2] - p[7]*N[1] - p[8]*N[2] - p[9]*N[3])
    dNdt[3] = N[3] * (p[3] - p[10]*N[1] - p[11]*N[2] - p[12]*N[3])
    return(list(dNdt))
  })
}

# posterior as matrix for predictions
post = as.matrix(fit_obs)
str(post)

n.post = nrow(post)
# n.post = 500 # uncomment to speed-up: lower number of samples used for predictions

# set time interval for prediction
timessim <- seq(min(sampling_times),max(sampling_times),len=100)

# matrices to save predictions
y1 = matrix(0, nrow=n.post, ncol=length(timessim))
y2 = y1
y3 = y1

# function for extracting CI and median of prediction
get_quantiles <- function(y){
  quantiles = matrix(0, nrow=3, ncol=ncol(y))
  quantiles[1, ] = apply(y, 2, function(x) quantile(x, probs=0.05))
  quantiles[2, ] = apply(y, 2, function(x) quantile(x, probs=0.50))
  quantiles[3, ] = apply(y, 2, function(x) quantile(x, probs=0.95))
  return(quantiles)
}


### Plots ###################################################################################

pdf("./out/3_species_fit.pdf")
## one species fit
par(mfrow = c(3,1))
for (k in 1:3){
  
  # loop over posterior samples
  for(l in 1:n.post){
    
    # ode simulation
    Nsim = as.data.frame(lsoda(y=c(post[l,paste0("N0_sim_Ble[", k, "]")],
                                   post[l,paste0("N0_sim_Col[", k, "]")],
                                   post[l,paste0("N0_sim_Tet[", k, "]")]),
                               times=timessim,
                               func=ode.model.mono,
                               parms=c(post[l,"r[1]"],
                                       post[l,"r[2]"],
                                       post[l,"r[3]"],
                                       post[l,"a[1,1]"],
                                       post[l,"a[2,1]"],
                                       post[l,"a[3,1]"],
                                       
                                       post[l,"a[1,2]"],
                                       post[l,"a[2,2]"],
                                       post[l,"a[3,2]"],
                                       
                                       post[l,"a[1,3]"],
                                       post[l,"a[2,3]"],
                                       post[l,"a[3,3]"])))[,c(2,3,4)]
    
    # save predictions for 3 species
    y1[l, ] = Nsim[,1]
    y2[l, ] = Nsim[,2]
    y3[l, ] = Nsim[,3]
    
  }
  
  # calculate median and CIs
  y1q = get_quantiles(y1) 
  y2q = get_quantiles(y2) 
  y3q = get_quantiles(y3) 
  
  plot(sampling_times, Ble[k,], col = "blue", xlab = "Time (hours)", ylab = "Blepharisma (mono)", ylim=c(0,225))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y1q[1, ], rev(y1q[3, ])),
          border=F,col=adjustcolor("blue",alpha.f=0.2))
  lines(timessim, y1q[2, ], col = "blue")
  
  plot(sampling_times, Col[k,], col = "red", xlab = "Time (hours)", ylab = "Colpidium (mono)", ylim=c(0,2500))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y2q[1, ], rev(y2q[3, ])),
          border=F,col=adjustcolor("red",alpha.f=0.2))
  lines(timessim, y2q[2, ], col = "red")
  
  plot(sampling_times, Tet[k,], col = "green", xlab = "Time (hours)", ylab = "Tetrahymena (mono)", ylim=c(0,8000))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y3q[1, ], rev(y3q[3, ])),
          border=F,col=adjustcolor("green",alpha.f=0.2))
  lines(timessim, y3q[2, ], col = "green")
  
}

## two species - Ble Col
par(mfrow = c(3,1))
for (k in 1:3){
  
  # loop over posterior samples
  for(l in 1:n.post){
    
    # ode simulation
    Nsim = as.data.frame(lsoda(y=c(post[l,paste0("N0_sim_Ble_w_col[", k, "]")],
                                   post[l,paste0("N0_sim_Col_w_ble[", k, "]")]),
                               times=timessim,
                               func=ode.model.ble.col,
                               parms=c(post[l,"r[1]"],
                                       post[l,"r[2]"],
                                       post[l,"r[3]"],
                                       post[l,"a[1,1]"],
                                       post[l,"a[2,1]"],
                                       post[l,"a[3,1]"],
                                       
                                       post[l,"a[1,2]"],
                                       post[l,"a[2,2]"],
                                       post[l,"a[3,2]"],
                                       
                                       post[l,"a[1,3]"],
                                       post[l,"a[2,3]"],
                                       post[l,"a[3,3]"])))[,c(2,3)]
    
    # save predictions for 2 species
    y1[l, ] = Nsim[,1]
    y2[l, ] = Nsim[,2]
  }
  
  # calculate median and CIs
  y1q = get_quantiles(y1) 
  y2q = get_quantiles(y2)
  
  plot(sampling_times, Ble_w_col[k,], col = "blue", xlab = "Time (hours)", ylab = "Blepharisma", ylim=c(0,100))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y1q[1, ], rev(y1q[3, ])),
          border=F,col=adjustcolor("blue",alpha.f=0.2))
  lines(timessim, y1q[2, ], col = "blue")
  
  plot(sampling_times, Col_w_ble[k,], col = "red", xlab = "Time (hours)", ylab = "Colpidium", ylim=c(0,2700))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y2q[1, ], rev(y2q[3, ])),
          border=F,col=adjustcolor("red",alpha.f=0.2))
  lines(timessim, y2q[2, ], col = "red")
  
  plot.new()
}

## two species - Ble Tet
par(mfrow = c(3,1))
for (k in 1:3){
  
  # loop over posterior samples
  for(l in 1:n.post){
    
    # ode simulation
    Nsim = as.data.frame(lsoda(y=c(post[l,paste0("N0_sim_Ble_w_tet[", k, "]")],
                                   post[l,paste0("N0_sim_Tet_w_ble[", k, "]")]),
                                   times=timessim,
                                   func=ode.model.ble.tet,
                                   parms=c(post[l,"r[1]"],
                                           post[l,"r[2]"],
                                           post[l,"r[3]"],
                                           post[l,"a[1,1]"],
                                           post[l,"a[2,1]"],
                                           post[l,"a[3,1]"],
                                           
                                           post[l,"a[1,2]"],
                                           post[l,"a[2,2]"],
                                           post[l,"a[3,2]"],
                                           
                                           post[l,"a[1,3]"],
                                           post[l,"a[2,3]"],
                                           post[l,"a[3,3]"])))[,c(2,3)]
                         
                         # save predictions for 2 species
                         y1[l, ] = Nsim[,1]
                         y2[l, ] = Nsim[,2]
  }
  
  # calculate median and CIs
  y1q = get_quantiles(y1) 
  y2q = get_quantiles(y2)
  
  plot(sampling_times, Ble_w_tet[k,], col = "blue", xlab = "Time (hours)", ylab = "Blepharisma", ylim=c(0,225))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y1q[1, ], rev(y1q[3, ])),
          border=F,col=adjustcolor("blue",alpha.f=0.2))
  lines(timessim, y1q[2, ], col = "blue")
  
  plot.new()
  
  plot(sampling_times, Tet_w_ble[k,], col = "green", xlab = "Time (hours)", ylab = "Tetrahymena", ylim=c(0,6750))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y2q[1, ], rev(y2q[3, ])),
          border=F,col=adjustcolor("green",alpha.f=0.2))
  lines(timessim, y2q[2, ], col = "green")
}


## three species fit
par(mfrow = c(3,1))
for (k in 1:3){
  
  # loop over posterior samples
  for(l in 1:n.post){
    
    # ode simulation
    Nsim = as.data.frame(lsoda(y=c(post[l,paste0("N0_sim_Ble_tri[", k, "]")],
                                   post[l,paste0("N0_sim_Col_tri[", k, "]")],
                                   post[l,paste0("N0_sim_Tet_tri[", k, "]")]),
                               times=timessim,
                               func=ode.model.tri,
                               parms=c(post[l,"r[1]"],
                                       post[l,"r[2]"],
                                       post[l,"r[3]"],
                                       post[l,"a[1,1]"],
                                       post[l,"a[2,1]"],
                                       post[l,"a[3,1]"],
                                       
                                       post[l,"a[1,2]"],
                                       post[l,"a[2,2]"],
                                       post[l,"a[3,2]"],
                                       
                                       post[l,"a[1,3]"],
                                       post[l,"a[2,3]"],
                                       post[l,"a[3,3]"])))[,c(2,3,4)]
    
    # save predictions for 3 species
    y1[l, ] = Nsim[,1]
    y2[l, ] = Nsim[,2]
    y3[l, ] = Nsim[,3]
    
  }
  
  # calculate median and CIs
  y1q = get_quantiles(y1) 
  y2q = get_quantiles(y2) 
  y3q = get_quantiles(y3) 
  
  plot(sampling_times, Ble_tri[k,], col = "blue", xlab = "Time (hours)", ylab = "Blepharisma", ylim=c(0,135))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y1q[1, ], rev(y1q[3, ])),
          border=F,col=adjustcolor("blue",alpha.f=0.2))
  lines(timessim, y1q[2, ], col = "blue")
  
  plot(sampling_times, Col_tri[k,], col = "red", xlab = "Time (hours)", ylab = "Colpidium", ylim=c(0,2200))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y2q[1, ], rev(y2q[3, ])),
          border=F,col=adjustcolor("red",alpha.f=0.2))
  lines(timessim, y2q[2, ], col = "red")
  
  plot(sampling_times, Tet_tri[k,], col = "green", xlab = "Time (hours)", ylab = "Tetrahymena", ylim=c(0,5500))
  polygon(x = c(timessim,rev(timessim)), 
          y= c(y3q[1, ], rev(y3q[3, ])),
          border=F,col=adjustcolor("green",alpha.f=0.2))
  lines(timessim, y3q[2, ], col = "green")
  
}

dev.off()
