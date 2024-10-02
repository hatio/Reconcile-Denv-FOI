functions {
    vector pow_vector(vector x, real y){
        vector[num_elements(x)] out;
        for(i in 1:num_elements(out)){
            out[i] = x[i] ^ y;
        }
        return out;
    }
    
    vector rowSums(matrix x){
        vector[rows(x)] out;
        for(i in 1:rows(x)){
            out[i] = sum(x[i, ]);
        }
        return out;
    }

}
data {  
    int<lower=1> num_times; // number of time periods to estimate FOI
    int<lower=1> num_indiv; // number of individuals
    
    matrix<lower=0>[num_indiv, num_times] Presence; // Time spent facing each FOI    
    int<lower=0,upper=1> obs[num_indiv]; // observed seropositivity
    
    // probability of showing positive given susceptibility state
    // S=0; S=1; S>1
    vector<lower=0,upper=1>[3] p_pos;

}
transformed data {
    
}

parameters {
    vector<lower=0,upper=1>[num_times] lambda;
}
transformed parameters {
}

model {
    // probability of escaping a particular serotype, by individual
    vector[num_indiv] p = exp(-Presence * lambda);
    // probability of being in each susceptibility state
    matrix[num_indiv, 3] S;
    
    S[ ,1] = pow_vector(p, 4);
    S[ ,2] = 4 * pow_vector(p, 3) .* (1-p);
    S[ ,3] = 1 - S[ ,1] - S[ ,2];
    
    // likelihoood    
    obs ~ bernoulli(S * p_pos);

    // priors
    lambda ~ exponential(2); // Very non-informative prior slightly favoring smaller values
}
generated quantities {

}
