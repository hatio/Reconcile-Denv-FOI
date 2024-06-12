# Parameters from Salje, 2018
Sigma = 0.49    # assay SD
Omega = 1.33        # permanent titer: mean
Omega_var = 1.85    # permanent titer: variance
Gamma = 5.39        # temporary titer rise: mean
Gamma_var = 4.54    # temporary titer rise: variance
Delta = 0.017 * 365.25  # decay rate of temporary titer per year

Omega_rate  = Omega/Omega_var
Omega_shape = Omega_rate * Omega


# Functions to work with these parameters
adjustTiter = function(x){
    log(x/10, 2) + 1
}
adjustTiter.inv = function(x){
    2^x * 10 / 2
}