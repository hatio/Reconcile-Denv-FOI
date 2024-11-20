functions {
    vector pow_vector(vector x, real y){
        vector[num_elements(x)] out;
        for(i in 1:num_elements(out)){
            out[i] = x[i] ^ y;
        }
        return out;
    }

    // function to expand binned variables based on given index
    row_vector expandBinnedByIndex( row_vector binned, int[] index){
        row_vector[size(index)] out;
        for(i in 1:size(index)){
            out[i] = binned[index[i]];
        }            
        return out;
    }
}
data {  
    int<lower=1> num_epoch; // number of epoches to estimate FOI
    // length of binned parameters
    int len_binned_kappa;
    int len_binned_phi_a;
    int len_binned_phi_t;    

    // ... case data
    int<lower=1> len_cohorts;
    int<lower=0> len_times_case;
    int<lower=1> epoch_time[len_cohorts]; // which epoch does each time index map to?
    int<lower=0> len_ages;
    int<lower=0> len_Ag;
    int obs_case[ len_Ag , len_times_case ];    // observed case counts
    real pop[ len_cohorts , len_times_case ];   // cohort population size
    matrix[ len_Ag , len_ages ] mapAge;         // mapping between Ag and age

    // index of binned parameters
    int<lower=1> index_kappa[len_ages];
    int<lower=1> index_phi_a[len_ages];
    int<lower=1> index_phi_t[len_times_case];
    
    // NA value
    int naVal;


}
transformed data {
    int<lower=0> len_pretimes = len_cohorts - len_times_case;
}

parameters {
    // infection parameters
    row_vector<lower=0,upper=1>[num_epoch] lambda;
    row_vector<lower=0.0>[len_binned_kappa - 1] binned_kappa_rel;
    
    // case reporting
    row_vector<lower=0.0>[len_binned_phi_a - 1] binned_phi_a_rel;
    row_vector<lower=0.0,upper=1.0>[len_binned_phi_t] binned_phi_t; // time-varying reporting rate
    vector<lower=0.0,upper=1.0>[3] p_dhf_rel; // proportion of i-th infection that is DHF (severe) relative to 2nd infection

    // over dispersion of cases
    real<lower=0> psi;

}
transformed parameters {
    // proportion of i-th infection that is DHF (severe)
    vector[4] p_dhf = [p_dhf_rel[1], 1.0, p_dhf_rel[2], p_dhf_rel[3]]';

    // expand binned parameters
    row_vector[len_cohorts] lambda_cohort = expandBinnedByIndex( lambda, epoch_time );
    row_vector[len_ages] kappa = expandBinnedByIndex(append_col(1.0, binned_kappa_rel), index_kappa);
    row_vector[len_ages] phi_a = expandBinnedByIndex(append_col(1.0, binned_phi_a_rel), index_phi_a);
    row_vector[len_times_case] phi_t = expandBinnedByIndex(binned_phi_t, index_phi_t);

}

model {    
    // .... case data
    // expected counts for each age
    matrix[ len_ages , len_times_case ] ageCountExpected = rep_matrix(0.0, len_ages , len_times_case);
    // expected counts for each age group
    matrix[ len_Ag , len_times_case ] countExpected;

    // for each birth cohort
    for (h in 1:len_cohorts){
        // maximum time index to evaluate for this birth cohort
        int Tmax = ((h - 1 + len_ages) > len_cohorts ? len_cohorts : h - 1 + len_ages);
        // susceptibility state: 0,1,2,3,4
        matrix[Tmax - h + 1 + 1, 5] S_cohort = rep_matrix(0.0, Tmax - h + 1 + 1, 5);
        // infection
        matrix[Tmax - h + 1, 4] I;
        // probability of this cohort escaping each serotype at each age
        row_vector[Tmax - h + 1] p_esc = exp(-cumulative_sum(lambda_cohort[h:Tmax] .* kappa[1:(Tmax - h + 1)]));
        
        S_cohort[1, 1] = 1.0; // everyone born susceptible
        for (a in 1:num_elements(p_esc)){
            for(i in 0:4){
                S_cohort[a+1, i+1] = (p_esc[a]^(4-i)) * (1-p_esc[a])^(i) * choose(4,i);
            }
        }

        // fraction with susceptibility > i
        for (a in 1:rows(S_cohort)){
            S_cohort[a, 1:4] = 1 - cumulative_sum(S_cohort[a, 1:4]) ;
        }
        
        // infections that happened = fraction susceptible > i that increased from last year
        I = S_cohort[ 2 : , 1:4] - S_cohort[ :(rows(S_cohort)-1) , 1:4];
        for (a in 1:rows(I)){ // age
            // time = h + a - 1
            int iTime = h + a - 1 - len_pretimes; // with(stanInput, len_cohorts - len_times_case)
            if(iTime < 1){ continue; }
            // compute expected cases from this age,year
            ageCountExpected[a,iTime] = pop[h, iTime] * I[a,] * p_dhf * phi_t[iTime] * phi_a[a];
        }
    }
    countExpected = mapAge * ageCountExpected;
    
    // likelihood: case data
    for(T in 1:len_times_case){
        for(ag in 1:len_Ag){
            if(obs_case[ag,T]==naVal){ continue; }
            obs_case[ag,T] ~ neg_binomial_2(countExpected[ag,T], psi) ;
        }
    }

    // priors: infection parameters
    lambda ~ exponential(2); // Very non-informative prior slightly favoring smaller values
    binned_kappa_rel ~ lognormal(0.0, 0.1); // E(x) close to 1
    // priors: case reporting and pathogenesis
    binned_phi_a_rel ~ lognormal(0.0, 0.1); // E(x) close to 1
    binned_phi_t ~ beta(2,2);  // E(x) = 0.5
    p_dhf_rel[1] ~ beta(1, 9);  // E(x) = 0.1 with zero being possible
    p_dhf_rel[2:3] ~ beta(1, 19);  // E(x) = 0.05 with zero being possible

    // priors: hyperparameter
    psi ~ normal(0, 50);

}
generated quantities {

    // expected values of observed data
    matrix[len_Ag , len_times_case] expect_case;    // case counts

    if(1){

        // .... case data
        // expected counts for each age
        matrix[ len_ages , len_times_case ] ageCountExpected = rep_matrix(0.0, len_ages , len_times_case);
        // expected counts for each age group
        matrix[ len_Ag , len_times_case ] countExpected;        

        // for each birth cohort
        for (h in 1:len_cohorts){
            // maximum time index to evaluate for this birth cohort
            int Tmax = ((h - 1 + len_ages) > len_cohorts ? len_cohorts : h - 1 + len_ages);
            // susceptibility state: 0,1,2,3,4
            matrix[Tmax - h + 1 + 1, 5] S_cohort = rep_matrix(0.0, Tmax - h + 1 + 1, 5);
            // infection
            matrix[Tmax - h + 1, 4] I;
            // probability of this cohort escaping each serotype at each age
            row_vector[Tmax - h + 1] p_esc = exp(-cumulative_sum(lambda_cohort[h:Tmax] .* kappa[1:(Tmax - h + 1)]));
            
            S_cohort[1, 1] = 1.0; // everyone born susceptible
            for (a in 1:num_elements(p_esc)){
                for(i in 0:4){
                    S_cohort[a+1, i+1] = (p_esc[a]^(4-i)) * (1-p_esc[a])^(i) * choose(4,i);
                }
            }

            // fraction with susceptibility > i
            for (a in 1:rows(S_cohort)){
                S_cohort[a, 1:4] = 1 - cumulative_sum(S_cohort[a, 1:4]) ;
            }
            
            // infections that happened = fraction susceptible > i that increased from last year
            I = S_cohort[ 2 : , 1:4] - S_cohort[ :(rows(S_cohort)-1) , 1:4];
            for (a in 1:rows(I)){ // age
                // time = h + a - 1
                int iTime = h + a - 1 - len_pretimes; // with(stanInput, len_cohorts - len_times_case)
                if(iTime < 1){ continue; }
                // compute expected cases from this age,year
                ageCountExpected[a,iTime] = pop[h, iTime] * I[a,] * p_dhf * phi_t[iTime] * phi_a[a];
            }
        }
        expect_case = mapAge * ageCountExpected;

    }


}
