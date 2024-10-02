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
    
    // function to generate kappa * presence row_vector from kappa of epoches and presence matrix
    row_vector kappaPresenceRowvector( matrix binned, row_vector mult, matrix Presence, int[] epochBreaks ){
        row_vector[cols(Presence)] out = rep_row_vector(0.0, cols(Presence));
        row_vector[cols(binned)] m = rep_row_vector(1.0, cols(binned));
        if(num_elements(mult)){ m = mult; }
        for(i in 1:rows(binned)){
            out[epochBreaks[i]:epochBreaks[i+1]] = (binned[i, ] .* m) * Presence[ , epochBreaks[i]:epochBreaks[i+1]];
        }
        return out;
    }
    
}
data {      
    int<lower=1> num_epoch; // number of epoches to estimate FOI
    int<lower=1> num_epoch_kappa; // number of epoches of Kappa(a) to estimate
    int<lower=0> num_thresh; // number of titer thresholds for seroprevalence
    int<lower=0> num_study;  // number of studies in this dataset
    // length of binned parameters
    int len_binned_kappa;
    int len_binned_phi_a;
    int len_binned_phi_t;    
    // whether to allow cohort studies to have risks that differs from the general population
    int<lower=0,upper=1> cohortKappa;
    // whether to have study-specific assay parameters
    int<lower=0,upper=1> cohortSigma;

    // ... seroincidence
    int<lower=0> num_intervals;
    int<lower=0> num_intervals_study[num_study];  // number of interval in each study
    int<lower=0> include_interval[num_intervals, num_thresh];
    matrix<lower=0>[len_binned_kappa, num_epoch] Presence[num_intervals]; // Time spent facing each FOI
    matrix<lower=0>[len_binned_kappa, num_epoch] prePresence[num_intervals]; // Time spent facing each FOI
    int<lower=0,upper=1> obs_inc[num_intervals, num_thresh]; // observed seropositivity

    // ... seroprevalence
    int<lower=0> num_indiv; // number of individuals
    int<lower=0> num_indiv_study[num_study]; // number of individuals by study
    matrix<lower=0>[len_binned_kappa, num_epoch] Presence_prev[num_indiv]; // Time spent facing each FOI
    int<lower=0,upper=1> obs_prev[num_indiv, num_thresh]; // observed seropositivity

    // ... case data
    int<lower=0> len_cohorts;
    int<lower=0> len_times_case;
    int<lower=1> epoch_time[len_cohorts]; // which epoch does each time index map to?
    int<lower=0> len_ages;
    int<lower=0> len_Ag;
    int obs_case[ len_Ag , len_times_case ];    // observed case counts
    real pop[ len_cohorts , len_times_case ];   // cohort population size
    matrix[ len_Ag , len_ages ] mapAge;         // mapping between Ag and age
    matrix<lower=0>[len_binned_kappa, len_cohorts] Presence_birthcohort[len_cohorts]; // Presence of birth cohort: Kappa age bin x time


    // index/breaks of binned parameters
    int<lower=0> break_epoch_kappa[2, num_epoch_kappa + 1]; // row 1: for cases, row 2: for serology
    int<lower=0> index_phi_a[len_ages];
    int<lower=0> index_phi_t[len_times_case];

    
    matrix<lower=0>[1,num_thresh] titer_thresh;     // titer threshold (log-scaled)
    real<lower=0> Omega_long_shape_prior;      // priors: permanent titer upon 1st infection
    real<lower=0> Omega_long_rate_prior;
    real<lower=0> Omega_short_shape_prior;      // priors: temporary titer upon 1st infection
    real<lower=0> Omega_short_rate_prior;
    real<lower=0> Sigma_prior;      // assay SD in exposed individuals

    // NA value
    int naVal;

}

transformed data {
    int<lower=0> len_pretimes = len_cohorts - len_times_case;
}

parameters {
    // infection parameters
    row_vector<lower=0,upper=1>[num_epoch] lambda;
    matrix<lower=0.0>[num_epoch_kappa, len_binned_kappa - 1] binned_kappa_rel; // two epoches of K(a) with K(a,epoch2) = K_rel(a,epoch1) * K_rel(a,epoch2)
    row_vector<lower=0.0>[(cohortKappa ? len_binned_kappa : 0)] kappa_study[num_study];  // risk multiplier compared to first study (reference class)
    
    // expected monotypic titer
    real<lower=0> Omega_long; // Long-lived Omega
    real<lower=0> Omega_short_ref; // Short-lived Omega of the first study (reference class)
    vector<lower=0>[(num_study > 1 ? num_study - 1 : 0)] Omega_short_mult;
    // expected cross-reactive titer in naive individuals
    real<lower=0,upper=1> Omega0_rel; // average CXR titer relative to long-lived DENV Omega
    
    // assay SD in exposed individuals
    real<lower=0> Sigma_ref;    // Sigma of first study (reference class)
    vector<lower=0>[(cohortSigma ? num_study - 1 : 0)] Sigma_rel;    // Sigma of study relative to reference class

    // case reporting
    row_vector<lower=0.0>[(len_binned_phi_a > 1 ? len_binned_phi_a - 1 : 0)] binned_phi_a_rel; // age-specific reporting relative to youngest age class
    row_vector<lower=0.0,upper=1.0>[len_binned_phi_t] binned_phi_t; // time-varying reporting rate
    vector<lower=0.0,upper=1.0>[3] p_dhf_rel; // proportion of i-th infection that is DHF (severe) relative to 2nd infection

    // over dispersion of cases
    real<lower=0> psi;
    
}
transformed parameters {
    // probability of showing positive given susceptibility state
    // S=0; S=1; S>1
    matrix[3, num_thresh] p_pos[num_study, 2];  // (study x (long-term, short-term)) array of matrix (S=0, S=1, S>1) x number of thresholds
    // // proportion of i-th infection that is DHF (severe) x p_dhf epoch
    // matrix[4,num_epoch_p_dhf] p_dhf = rep_matrix(1.0, 4, num_epoch_p_dhf);
    vector[num_study] Omega_short = append_row(1.0, Omega_short_mult) * Omega_short_ref + Omega_long;
    real Omega0 = Omega_long * Omega0_rel;
    vector[num_study] Sigma = append_row(1.0, (cohortSigma ? Sigma_rel : rep_vector(1.0, (num_study > 0 ? num_study - 1 : 0)))) * Sigma_ref;

    // expand binned parameters
    row_vector[len_cohorts] lambda_birthcohort = expandBinnedByIndex( lambda, epoch_time ); // infection risk in case data
    row_vector[len_ages] phi_a = expandBinnedByIndex(append_col(1.0, binned_phi_a_rel), index_phi_a);
    row_vector[len_times_case] phi_t = expandBinnedByIndex(binned_phi_t, index_phi_t);

    // expand binned kappa
    matrix[num_epoch_kappa, len_binned_kappa] binned_kappa;
    binned_kappa[1, ] = append_col(1.0, binned_kappa_rel[1, ]);
    for(i in 2:num_epoch_kappa){
        binned_kappa[i, ] = binned_kappa[1, ] .* append_col(1.0, binned_kappa_rel[i, ]);
    }

    // probability of testing positive for cross-sectional bleeds
    // and at 2nd bleed of each study    
    for(i in 1:num_study){
        // long term
        p_pos[i,1][1, ] = 1.0 - Phi((titer_thresh[1,] - Omega0    ) / Sigma[i]);
        p_pos[i,1][2, ] = 1.0 - Phi((titer_thresh[1,] - Omega_long) / Sigma[i]);
        p_pos[i,1][3, ] = rep_row_vector(1.0, num_thresh);

        // short term
        p_pos[i,2][1, ] = 1.0 - Phi((titer_thresh[1,] - Omega0        ) / Sigma[i]);
        p_pos[i,2][2, ] = 1.0 - Phi((titer_thresh[1,] - Omega_short[i]) / Sigma[i]);
        p_pos[i,2][3, ] = rep_row_vector(1.0, num_thresh);
    }

    // proportion of i-th infection that is DHF (severe)
    vector[4] p_dhf = [p_dhf_rel[1], 1.0, p_dhf_rel[2], p_dhf_rel[3]]';

}

model {
    int pointer;
    // expected values of observed data
    matrix[num_indiv, num_thresh] expect_prev;      // seropositivity
    matrix[num_intervals, num_thresh] expect_inc;   // seroconversion
    matrix[len_Ag , len_times_case] expect_case;    // case counts
    
    // .... seroincidence
    vector[num_intervals] pPre;
    vector[num_intervals] p;
    matrix[num_intervals, 3] Spre;
    matrix[num_intervals, 3] S0;
    matrix[num_intervals, 3] S1;

    // .... seroprevalence
    vector[num_indiv] p_prev;
    matrix[num_indiv, 3] S;

    // .... case data
    // expected counts for each age
    matrix[ len_ages , len_times_case ] ageCountExpected = rep_matrix(0.0, len_ages , len_times_case);
    // expected counts for each age group
    matrix[ len_Ag , len_times_case ] countExpected;        


    pointer = 1;
    for(i in 1:num_intervals){
        // probability of escaping a particular serotype BEFORE the interval, by individual
        pPre[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], prePresence[i], break_epoch_kappa[2, ]) * lambda');
        // probability of escaping a particular serotype DURING the interval, by individual
        p[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence[i], break_epoch_kappa[2, ]) * lambda');
        if(i == sum(num_intervals_study[ :pointer])){ pointer += 1; }
    }
    if(num_intervals){
        // probability of being in each susceptibility state at START of the interval
        Spre[ ,1] = pow_vector(pPre, 4);
        Spre[ ,2] = 4 * pow_vector(pPre, 3) .* (1-pPre);
        Spre[ ,3] = 1 - Spre[ ,1] - Spre[ ,2];
        // probability of being in each susceptibility state at END of the interval, given S_t1 = 0
        S0[ ,1] = pow_vector(p, 4);
        S0[ ,2] = 4 * pow_vector(p, 3) .* (1-p);
        S0[ ,3] = 1 - S0[ ,1] - S0[ ,2];
        // probability of being in each susceptibility state at END of the interval, given S_t1 = 1
        S1[ ,1] = rep_vector(0.0, num_intervals);
        S1[ ,2] = pow_vector(p, 3);
        S1[ ,3] = 1 - S1[ ,1] - S1[ ,2];        
    }

    
    pointer = 1;
    for(i in 1:num_indiv){
        // probability of escaping a particular serotype, by individual
        p_prev[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence_prev[i], break_epoch_kappa[2, ]) * lambda');
        if(i == sum(num_indiv_study[ :pointer])){ pointer += 1; }
    }
    if(num_indiv){
        // probability of being in each susceptibility state
        S[ ,1] = pow_vector(p_prev, 4);
        S[ ,2] = 4 * pow_vector(p_prev, 3) .* (1-p_prev);
        S[ ,3] = 1 - S[ ,1] - S[ ,2];        
    }


    // for each birth cohort in case data
    for (h in 1:len_cohorts){
        // maximum time index to evaluate for this birth cohort
        int Tmax = ((h - 1 + len_ages) > len_cohorts ? len_cohorts : h - 1 + len_ages);
        // susceptibility state: 0,1,2,3,4
        matrix[Tmax - h + 1 + 1, 5] S_cohort = rep_matrix(0.0, Tmax - h + 1 + 1, 5);
        // infection
        matrix[Tmax - h + 1, 4] I;
        // probability of this cohort escaping each serotype at each age
        row_vector[Tmax - h + 1] p_esc = exp(-cumulative_sum(
            lambda_birthcohort[h:Tmax] .* 
            kappaPresenceRowvector(binned_kappa, rep_row_vector(1.0, len_binned_kappa), Presence_birthcohort[h], break_epoch_kappa[1, ])[h:Tmax]            
        ));
        
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
    

    // calculate expected values: seroincidence
    pointer = 1;
    for(i in 1:num_study){
        if(num_intervals_study[i] == 0){ continue; }
        for(it in 1:num_thresh){
            expect_inc[pointer : (pointer + num_intervals_study[i] - 1), it] =
                (((Spre[pointer : (pointer + num_intervals_study[i] - 1),1]  * (1-p_pos[i,1][1  ,it])) .* (S0[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[i,2][ ,it])) +
                ((Spre[pointer : (pointer + num_intervals_study[i] - 1),2:3] * (1-p_pos[i,1][2:3,it])) .* (S1[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[i,2][ ,it]))) ./ (Spre[pointer : (pointer + num_intervals_study[i] - 1), ] * (1-p_pos[i,1][ ,it]));
        }
        pointer += num_intervals_study[i];
    }


    // calculate expected values: seroprevalence
    pointer = 1;
    for(i in 1:num_study){
        for(it in 1:num_thresh){
            expect_prev[pointer : (pointer + num_indiv_study[i] - 1) ,it] = 
                S[pointer : (pointer + num_indiv_study[i] - 1), ] * p_pos[i,1][ ,it];
        }
        pointer += num_indiv_study[i];
    }


    // calculate expected values: case data
    if(len_cohorts){
        expect_case = mapAge * ageCountExpected;    
    }
    
  
    // likelihood: seroincidence
    for(it in 1:num_thresh){
        for(i in 1:num_intervals){
            if(include_interval[i,it]){
                obs_inc[i,it] ~ bernoulli(expect_inc[i,it]);    
            }
        }
    }

    // likelihoood: seroprevalence
    for(it in 1:num_thresh){
        obs_prev[,it] ~ bernoulli(expect_prev[,it]);
    }

    // likelihood: case data
    for(T in 1:len_times_case){
        for(ag in 1:len_Ag){
            if(obs_case[ag,T]==naVal){ continue; }
            obs_case[ag,T] ~ neg_binomial_2(expect_case[ag,T], psi) ;
        }
    }


    // priors: infection parameters
    lambda ~ exponential(2); // Very weak prior slightly favoring smaller values
    binned_kappa_rel[1,] ~ lognormal(0.0, 0.1); // E(x) close to 1
    for(i in 2:rows(binned_kappa_rel)){
        binned_kappa_rel[i,] ~ lognormal(-0.045, 0.3); // E(x) close to 1, but with more relaxed priors
    }
    for(i in 1:num_study){
        kappa_study[i] ~ lognormal(-0.045, 0.3); // E(x) close to 1, but with more relaxed priors
    }

    
    // priors: assay
    Omega_long ~ gamma(Omega_long_shape_prior, Omega_long_rate_prior);
    Omega_short_ref ~ gamma(Omega_short_shape_prior, Omega_short_rate_prior);
    Omega_short_mult ~ lognormal(0.0, 0.1); // E(x) close to 1
    Omega0_rel ~ beta(1,9); // E(x) = 0.1
    Sigma_ref ~ normal(Sigma_prior, 0.1) T[0,];
    Sigma_rel ~ lognormal(0.0, 0.1); // E(x) close to 1


    // priors: case reporting and pathogenesis
    binned_phi_a_rel ~ lognormal(0.0, 0.1); // E(x) close to 1
    binned_phi_t ~ beta(2,2);  // E(x) = 0.5
    p_dhf_rel[1] ~ beta(1, 9);  // E(x) = 0.1 with zero being possible
    p_dhf_rel[2:3] ~ beta(1, 19);  // E(x) = 0.05 with zero being possible

    // priors: hyperparameter
    psi ~ normal(0, 100);
}
generated quantities {
    // // expected values of observed data
    // matrix[num_indiv, num_thresh] expect_prev; // seropositivity
    // matrix[num_intervals, num_thresh] expect_inc; // seroconversion

    // // log-likelihood for each data point
    // matrix[num_indiv, num_thresh] log_lik_prev = rep_matrix(0.0, num_indiv, num_thresh); // seropositivity
    // matrix[num_intervals, num_thresh] log_lik_inc = rep_matrix(0.0, num_intervals, num_thresh); // seroconversion

    // if(1){
    //     int pointer;
        
    //     // .... seroincidence
    //     vector[num_intervals] pPre;
    //     vector[num_intervals] p;
    //     matrix[num_intervals, 3] Spre;
    //     matrix[num_intervals, 3] S0;
    //     matrix[num_intervals, 3] S1;

    //     // .... seroprevalence
    //     vector[num_indiv] p_prev;
    //     matrix[num_indiv, 3] S;


    //     pointer = 1;
    //     for(i in 1:num_intervals){
    //         // probability of escaping a particular serotype BEFORE the interval, by individual
    //         pPre[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], prePresence[i], break_epoch_kappa[2, ]) * lambda');
    //         // probability of escaping a particular serotype DURING the interval, by individual
    //         p[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence[i], break_epoch_kappa[2, ]) * lambda');
    //         if(i == sum(num_intervals_study[ :pointer])){ pointer += 1; }
    //     }
    //     // probability of being in each susceptibility state at START of the interval
    //     Spre[ ,1] = pow_vector(pPre, 4);
    //     Spre[ ,2] = 4 * pow_vector(pPre, 3) .* (1-pPre);
    //     Spre[ ,3] = 1 - Spre[ ,1] - Spre[ ,2];
    //     // probability of being in each susceptibility state at END of the interval, given S_t1 = 0
    //     S0[ ,1] = pow_vector(p, 4);
    //     S0[ ,2] = 4 * pow_vector(p, 3) .* (1-p);
    //     S0[ ,3] = 1 - S0[ ,1] - S0[ ,2];
    //     // probability of being in each susceptibility state at END of the interval, given S_t1 = 1
    //     S1[ ,1] = rep_vector(0.0, num_intervals);
    //     S1[ ,2] = pow_vector(p, 3);
    //     S1[ ,3] = 1 - S1[ ,1] - S1[ ,2];

        
    //     pointer = 1;
    //     for(i in 1:num_indiv){
    //         // probability of escaping a particular serotype, by individual
    //         p_prev[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence_prev[i], break_epoch_kappa[2, ]) * lambda');
    //         if(i == sum(num_indiv_study[ :pointer])){ pointer += 1; }
    //     }
    //     // probability of being in each susceptibility state
    //     S[ ,1] = pow_vector(p_prev, 4);
    //     S[ ,2] = 4 * pow_vector(p_prev, 3) .* (1-p_prev);
    //     S[ ,3] = 1 - S[ ,1] - S[ ,2];


    //     // calculate expected values: seroincidence
    //     pointer = 1;
    //     for(i in 1:num_study){
    //         if(num_intervals_study[i] == 0){ continue; }
    //         for(it in 1:num_thresh){
    //             expect_inc[pointer : (pointer + num_intervals_study[i] - 1), it] =
    //                 (((Spre[pointer : (pointer + num_intervals_study[i] - 1),1]  * (1-p_pos[i,1][1  ,it])) .* (S0[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[i,2][ ,it])) +
    //                 ((Spre[pointer : (pointer + num_intervals_study[i] - 1),2:3] * (1-p_pos[i,1][2:3,it])) .* (S1[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[i,2][ ,it]))) ./ (Spre[pointer : (pointer + num_intervals_study[i] - 1), ] * (1-p_pos[i,1][ ,it]));
    //         }
    //         pointer += num_intervals_study[i];
    //     }


    //     // calculate expected values: seroprevalence
    //     pointer = 1;
    //     for(i in 1:num_study){
    //         for(it in 1:num_thresh){
    //             expect_prev[pointer : (pointer + num_indiv_study[i] - 1) ,it] = 
    //                 S[pointer : (pointer + num_indiv_study[i] - 1), ] * p_pos[i,1][ ,it];
    //         }
    //         pointer += num_indiv_study[i];
    //     }
        

    //     // likelihood: seroincidence
    //     for(i in 1:rows(expect_inc)){
    //         for(it in 1:num_thresh){
    //             if(include_interval[i,it]){
    //                 log_lik_inc[i,it] = bernoulli_lpmf(obs_inc[i,it] | expect_inc[i,it]);
    //             }
    //         }
    //     }


    //     // likelihoood: seroprevalence
    //     for(i in 1:rows(expect_prev)){
    //         for(it in 1:num_thresh){
    //             log_lik_prev[i,it] = bernoulli_lpmf(obs_prev[i,it] | expect_prev[i,it]);
    //         }            
    //     }
        
    // }
}

