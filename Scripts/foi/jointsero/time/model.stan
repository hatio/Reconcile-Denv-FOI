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
    int<lower=1> num_thresh; // number of titer thresholds for seroprevalence
    int<lower=1> num_study;  // number of studies in this dataset
    // length of binned parameters
    int len_binned_kappa;
    // whether to allow cohort studies to have risks that differs from the general population
    int<lower=0,upper=1> cohortKappa;
    // whether to have study-specific assay parameters
    int<lower=0,upper=1> cohortSigma;

    // ... seroincidence
    int<lower=1> num_intervals;
    int<lower=0> num_intervals_study[num_study];  // number of interval in each study
    int<lower=0> include_interval[num_intervals, num_thresh];
    matrix<lower=0>[len_binned_kappa, num_epoch] Presence[num_intervals]; // Time spent facing each FOI
    matrix<lower=0>[len_binned_kappa, num_epoch] prePresence[num_intervals]; // Time spent facing each FOI
    int<lower=0,upper=1> obs_inc[num_intervals, num_thresh]; // observed seropositivity

    // ... seroprevalence
    int<lower=1> num_indiv; // number of individuals
    int<lower=1> num_indiv_study[num_study]; // number of individuals by study
    matrix<lower=0>[len_binned_kappa, num_epoch] Presence_prev[num_indiv]; // Time spent facing each FOI
    int<lower=0,upper=1> obs_prev[num_indiv, num_thresh]; // observed seropositivity


    // // index/breaks of binned parameters
    int<lower=1> break_epoch_kappa[2, num_epoch_kappa + 1]; // row 1: for cases, row 2: for serology
    
    matrix<lower=0>[1,num_thresh] titer_thresh;     // titer threshold (log-scaled)
    real<lower=0> Omega_long_shape_prior;      // priors: permanent titer upon 1st infection
    real<lower=0> Omega_long_rate_prior;
    real<lower=0> Omega_short_shape_prior;      // priors: temporary titer upon 1st infection
    real<lower=0> Omega_short_rate_prior;
    real<lower=0> Sigma_prior;      // assay SD in exposed individuals

}

transformed data {

}

parameters {
    // infection parameters
    row_vector<lower=0,upper=1>[num_epoch] lambda;
    matrix<lower=0.0>[num_epoch_kappa, len_binned_kappa - 1] binned_kappa_rel; // two epoches of K(a) with K(a,epoch2) = K_rel(a,epoch1) * K_rel(a,epoch2)
    row_vector<lower=0.0>[(cohortKappa ? len_binned_kappa : 0)] kappa_study[num_study];  // risk multiplier compared to first study (reference class)
    
    // expected monotypic titer
    real<lower=0> Omega_long; // Long-lived Omega
    real<lower=0> Omega_short_ref; // Short-lived Omega of the first study (reference class)
    vector<lower=0>[num_study - 1] Omega_short_mult;
    // expected cross-reactive titer in naive individuals
    real<lower=0,upper=1> Omega0_rel; // average CXR titer relative to long-lived DENV Omega
    
    // assay SD in exposed individuals
    real<lower=0> Sigma_ref;    // Sigma of first study (reference class)
    vector<lower=0>[(cohortSigma ? num_study - 1 : 0)] Sigma_rel;    // Sigma of study relative to reference class
    
}
transformed parameters {
    // probability of showing positive given susceptibility state
    // S=0; S=1; S>1
    matrix[3, num_thresh] p_pos[num_study, 2];  // (study x (long-term, short-term)) array of matrix (S=0, S=1, S>1) x number of thresholds
    // // proportion of i-th infection that is DHF (severe) x p_dhf epoch
    // matrix[4,num_epoch_p_dhf] p_dhf = rep_matrix(1.0, 4, num_epoch_p_dhf);
    vector[num_study] Omega_short = append_row(1.0, Omega_short_mult) * Omega_short_ref + Omega_long;
    real Omega0 = Omega_long * Omega0_rel;
    vector[num_study] Sigma = append_row(1.0, (cohortSigma ? Sigma_rel : rep_vector(1.0, num_study - 1))) * Sigma_ref;

    // expand binned parameters
    matrix[num_epoch_kappa, len_binned_kappa] binned_kappa;
    
    // kappa
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

}

model {
    int pointer;
    // expected values of observed data
    matrix[num_indiv, num_thresh] expect_prev;      // seropositivity
    matrix[num_intervals, num_thresh] expect_inc;   // seroconversion
    // matrix[len_Ag , len_times_case] expect_case;    // case counts
    
    // .... seroincidence
    vector[num_intervals] pPre;
    vector[num_intervals] p;
    matrix[num_intervals, 3] Spre;
    matrix[num_intervals, 3] S0;
    matrix[num_intervals, 3] S1;

    // .... seroprevalence
    vector[num_indiv] p_prev;
    matrix[num_indiv, 3] S;


    pointer = 1;
    for(i in 1:num_intervals){
        // probability of escaping a particular serotype BEFORE the interval, by individual
        pPre[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], prePresence[i], break_epoch_kappa[2, ]) * lambda');
        // probability of escaping a particular serotype DURING the interval, by individual
        p[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence[i], break_epoch_kappa[2, ]) * lambda');
        if(i == sum(num_intervals_study[ :pointer])){ pointer += 1; }
    }
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

    
    pointer = 1;
    for(i in 1:num_indiv){
        // probability of escaping a particular serotype, by individual
        p_prev[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence_prev[i], break_epoch_kappa[2, ]) * lambda');
        if(i == sum(num_indiv_study[ :pointer])){ pointer += 1; }
    }
    // probability of being in each susceptibility state
    S[ ,1] = pow_vector(p_prev, 4);
    S[ ,2] = 4 * pow_vector(p_prev, 3) .* (1-p_prev);
    S[ ,3] = 1 - S[ ,1] - S[ ,2];


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


    // priors: infection parameters
    lambda ~ beta( 2, 38); // E(x) = 0.05
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

}
generated quantities {
    // expected values of observed data
    matrix[num_indiv, num_thresh] expect_prev; // seropositivity
    matrix[num_intervals, num_thresh] expect_inc; // seroconversion

    // log-likelihood for each data point
    matrix[num_indiv, num_thresh] log_lik_prev = rep_matrix(0.0, num_indiv, num_thresh); // seropositivity
    matrix[num_intervals, num_thresh] log_lik_inc = rep_matrix(0.0, num_intervals, num_thresh); // seroconversion

    if(1){
        int pointer;
        
        // .... seroincidence
        vector[num_intervals] pPre;
        vector[num_intervals] p;
        matrix[num_intervals, 3] Spre;
        matrix[num_intervals, 3] S0;
        matrix[num_intervals, 3] S1;

        // .... seroprevalence
        vector[num_indiv] p_prev;
        matrix[num_indiv, 3] S;


        pointer = 1;
        for(i in 1:num_intervals){
            // probability of escaping a particular serotype BEFORE the interval, by individual
            pPre[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], prePresence[i], break_epoch_kappa[2, ]) * lambda');
            // probability of escaping a particular serotype DURING the interval, by individual
            p[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence[i], break_epoch_kappa[2, ]) * lambda');
            if(i == sum(num_intervals_study[ :pointer])){ pointer += 1; }
        }
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

        
        pointer = 1;
        for(i in 1:num_indiv){
            // probability of escaping a particular serotype, by individual
            p_prev[i] = exp(-kappaPresenceRowvector(binned_kappa, kappa_study[pointer], Presence_prev[i], break_epoch_kappa[2, ]) * lambda');
            if(i == sum(num_indiv_study[ :pointer])){ pointer += 1; }
        }
        // probability of being in each susceptibility state
        S[ ,1] = pow_vector(p_prev, 4);
        S[ ,2] = 4 * pow_vector(p_prev, 3) .* (1-p_prev);
        S[ ,3] = 1 - S[ ,1] - S[ ,2];


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
        

        // likelihood: seroincidence
        for(i in 1:rows(expect_inc)){
            for(it in 1:num_thresh){
                if(include_interval[i,it]){
                    log_lik_inc[i,it] = bernoulli_lpmf(obs_inc[i,it] | expect_inc[i,it]);
                }
            }
        }


        // likelihoood: seroprevalence
        for(i in 1:rows(expect_prev)){
            for(it in 1:num_thresh){
                log_lik_prev[i,it] = bernoulli_lpmf(obs_prev[i,it] | expect_prev[i,it]);
            }            
        }
        
    }
}

