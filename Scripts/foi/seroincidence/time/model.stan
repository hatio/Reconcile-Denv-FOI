functions {
    vector pow_vector(vector x, real y){
        vector[num_elements(x)] out;
        for(i in 1:num_elements(out)){
            out[i] = x[i] ^ y;
        }
        return out;
    }
    
}
data {  
    int<lower=1> num_intervals;
    int<lower=1> num_times;
    
    int<lower=1> num_study;  // number of studies in this dataset
    int<lower=0> num_intervals_study[num_study];  // number of interval in each study
    
    matrix<lower=0>[num_intervals, num_times] Presence;
    matrix<lower=0>[num_intervals, num_times] prePresence;
    int<lower=0,upper=1> obs[num_intervals]; // observed seropositivity

    // probability of showing positive given susceptibility state
    // (S=0; S=1; S>1) x (long-term, short-term given S_t1=0 by study)
    matrix<lower=0,upper=1>[3, 1 + num_study] p_pos;

}
transformed data {
    
}

parameters {
    vector<lower=0,upper=1>[num_times] lambda;
}
transformed parameters {
}

model {
    int pointer;
    // probability of escaping a particular serotype BEFORE the interval, by individual
    vector[num_intervals] pPre = exp(-prePresence * lambda);
    // probability of escaping a particular serotype DURING the interval, by individual
    vector[num_intervals] p = exp(-Presence * lambda);
    // probability of being in each susceptibility state at START of the interval
    matrix[num_intervals, 3] Spre;
    // probability of being in each susceptibility state at END of the interval, given S_t1 = 0
    matrix[num_intervals, 3] S0;
    // probability of being in each susceptibility state at END of the interval, given S_t1 = 1
    matrix[num_intervals, 3] S1;
    
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


    // posterior
    // expected seroincidence
    pointer = 1;
    for(i in 1:num_study){
        if(num_intervals_study[i] == 0){ continue; }
        obs[pointer : (pointer + num_intervals_study[i] - 1)] ~ bernoulli((
            ((Spre[pointer : (pointer + num_intervals_study[i] - 1) ,1]   * (1-p_pos[1  ,1]))   .* (S0[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[ ,1+i])) +
            ((Spre[pointer : (pointer + num_intervals_study[i] - 1) ,2:3] * (1-p_pos[2:3,1]))   .* (S1[pointer : (pointer + num_intervals_study[i] - 1), ] * p_pos[ ,1  ]))
        ) ./ (Spre[pointer : (pointer + num_intervals_study[i] - 1), ] * (1-p_pos[ ,1]))
        );
        pointer += num_intervals_study[i];
    }

    // prior
    lambda ~ beta( 2, 38); // E(x) = 0.05
}
generated quantities {

}
