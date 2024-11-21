#' ---
#' title: "Supplementary mathematical analysis: biases in force of infection inferred from seroincidence data"
#' author: ""
#' date: ""
#' header-includes:
#'   \pagenumbering{gobble}
#'   \setcounter{table}{0}  \renewcommand{\thetable}{S\arabic{table}} \setcounter{figure}{11} \renewcommand{\thefigure}{S\arabic{figure}}
#' ---


#+ include = FALSE
knitr::opts_chunk$set(include = TRUE, echo = FALSE, warning = FALSE)


library(tidyverse); theme_set(theme_classic())
library(ggpubr)


# Salje, 2018 estimates
Omega = 1.33
Gamma = 5.39
Delta = 0.017 * 365.25
Sigma = 0.49
# Arbitrary relative amount of CXR titer
Omega0_rel = 0.2



#' Notations used throughout this document are defined as follows.

#' \begin{align*}
#' P(\ominus_t)     &= \text{Probability of testing negative at time $t$} \\
#' P(\oplus_t)      &= \text{Probability of testing positive at time $t$} \\
#' P(\ominus^T_t)   &= \text{Probability of being true negative at time $t$} \\
#' P(\oplus^T_t)    &= \text{Probability of being true positive at time $t$} \\
#' P(\ominus^F_t)   &= \text{Probability of being false negative at time $t$} \\
#' P(\oplus^F_t)    &= \text{Probability of being false positive at time $t$} \\
#' % ---- suspectibiliyy
#' P(S_{t} = i)     &= \text{Probability of having acquired $i$ infections at time $t$} \\
#' p_{a}            &= \text{Probability of escaping a particular serotype up till age $a$} \\
#' \bar{\lambda}    &= \text{Average historical per-serotype force of infection} \\
#' \lambda_t        &= \text{Per-serotype force of infection during time interval $(t, t + \Delta{t}]$} \\
#' % ---- titer-related
#' \omega           &= \text{Long-term persistent titer rise acquired upon primary dengue infection} \\
#' \gamma           &= \text{Temporary titer rise upon primary dengue infection} \\
#' \delta           &= \text{Rate of exponential decay in temporary titers} \\
#' \sigma           &= \text{Standard deviation of assay measurements from true underlying titer of a DENV-exposed individual} \\
#' \sigma_0         &= \text{Standard deviation of assay measurements from true underlying titer of a DENV-naive individual} \\
#' \nu              &= \text{Seropositivity titer cut-off} \\
#' \Phi(x)          &= \text{Cumulative density of a standard normal distribution at value $x$} \\
#' \end{align*}
#' 



#' Longitudinal serology is typically viewed as the gold standard for measuring levels of transmission within a time period as it directly measures serological changes of individuals. 
#' We can express the
#' \textbf{true probability of first infection incidence} during a time interval $(t, t + \Delta{t}]$ as

#' \begin{equation}
#' \begin{split}
#' P(S_{t + \Delta{t}} \ge 1 | S_{t} = 0) = P(S_{t + \Delta{t}} = 1 | S_{t} = 0) + P(S_{t + \Delta{t}} > 1 | S_{t} = 0)
#' \end{split}
#' \end{equation}
#'


#' \textbf{Probability of observing seroconversion} in an individual, $P(\oplus_{(t, t + \Delta{t}]} | \ominus_{t})$, can be expanded to a weighted average between probabilities of seroconverting given being truely negative at the first time point, $t$, and being falsely negative. 
#'

#' \begin{equation} \label{eq:seroconvert}
#' \begin{split}
#' P(\oplus_{(t, t + \Delta{t}]} | \ominus_{t}) =
#'       \frac{P( \ominus^T_{t} )}{P(\ominus_{t})}{[ P( \oplus^T_{t + \Delta{t}} | \ominus^T_{t} ) + P( \oplus^F_{t + \Delta{t}} | \ominus^T_{t} ) ]} +
#'       \frac{P( \ominus^F_{t} )}{P(\ominus_{t})}{[ P( \oplus^T_{t + \Delta{t}} | \ominus^F_{t} ) + P( \oplus^F_{t + \Delta{t}} | \ominus^F_{t} ) ]}
#' \end{split}
#' \end{equation}
#'





#' \textbf{Biases in seroconversion probabilities when truly negative at t}
#' from the true probability of first infection incidence is then
#  ..................

#' \begin{equation} \label{eq:bias_tneg}
#' \begin{split}
#' & [ P( \oplus^T_{(t, t + \Delta{t}]} | \ominus^T_{t} ) + P( \oplus^F_{(t, t + \Delta{t}]} | \ominus^T_{t} ) ]
#'     - P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) 
#'     - P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) \\
#' & =   P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) ~ P(\oplus^T_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} = 1)
#'       + P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) ~ P(\oplus^T_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} > 1) \\
#' &     + P(S_{(t, t + \Delta{t}]} = 0 | S_{t} = 0) ~ P(\oplus^F_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} = 0)
#'       - P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) 
#'       - P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) \\
#' & =   P(S_{(t, t + \Delta{t}]} = 0 | S_{t} = 0) ~ P(\oplus^F_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} = 0) \\
#' &     + P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) [ P(\oplus^T_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} = 1) - 1 ] \\
#' &     + P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) [ P(\oplus^T_{(t, t + \Delta{t}]} | S_{(t, t + \Delta{t}]} > 1) - 1 ] \\
#' \end{split}
#' \end{equation}
#'

#' We can see that the bias is a weighted average between biases when no infection occurred during the interval ($S_{(t, t + \Delta{t}]} = 0, S_{t} = 0$), when one infection occurred during the interval ($S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0$), and when more than one infection occurred during the interval ($S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0$).
#' Assuming long-lived protection against infecting serotypes and no cross-protection between the serotypes, we can write the probabilities of these scenarios (i.e., the weights) for an individual of age $a$ at time $(t, t + \Delta{t}]$ as
#'

#' \begin{equation}
#' P(S_{(t, t + \Delta{t}]} = i) = {4 \choose i} ~ (e^{-\bar{\lambda}a})^{4-i} ~ (1 - e^{-\bar{\lambda}a})^{i}
#' \end{equation}
#'

#' Previous studies on the kinetics of anti-DENV antibodies have shown that the rise in antibody levels after having acquired one infection, i.e., monotypic titers, consists of a portion which diminishes within a year and  a portion which persists for longer (cite: Leah, Henrik). Following Salje et al, the kinetics of $H_{t}$, the true titer at time $t$, can be characterized as 
#' 

#' \begin{equation}
#' H_t = \omega + \gamma e^{-(t - \star{t}) \delta}
#' \end{equation}
#'

#' where $\omega$ is the long-term persistent titer acquired upon infection, $\gamma$ is the temporary titer rise, $\delta$ is the rate of exponential decay in temporary titer rise, and $\star{t}$ is time at which the infection occurred. For an assay with measurement noise which follows a normal distribution $N(0, \sigma^2)$,
#' using seropositivity cut-off $\nu$,
#' the probability of testing positive after being infected once is 
#'

#' \begin{equation}
#' P(\oplus^T_t | S_t = 1)
#' = 1 - \Phi(\frac{\nu - H_{t}}{\sigma})
#' \end{equation}
#'


#' In DENV-naive individuals, although anti-DENV antibodies are absent, cross-reactive titers against other \textit{Flaviviruses} in circulation could plausibly result in non-zero titers against DENV, $\omega_0$.
#' We can, in a similar manner, express the probability of testing positive when DENV-naive as
#'

#' \begin{equation}
#' P(\oplus^F_t | S_t = 0)
#' = 1 - \Phi(\frac{\nu - \omega_0}{\sigma})
#' \end{equation}
#'


#' Let probability of escaping a particular serotype during time interval $e^{-\lambda_t \Delta{t}} = p$,
#' we can express (\ref{eq:bias_tneg}) as

#' \begin{equation} \label{eq:bias_tneg_p}
#' \begin{split}
#' & p^4 ~ [ 1 - \Phi({\frac{\nu - \omega_0}{\sigma}}) ]
#'     + {4 \choose 1} ~ p^3 ~ (1-p) ~ [ (1 - \Phi({\frac{\nu - H_t}{\sigma}})) - 1 ]
#'     + [1 - p^4 - {4 \choose 1} ~ p^3 ~ (1-p) ] ~ [ 1 - 1 ] \\
#' & = p^4 ~ [ 1 - \Phi({\frac{\nu - \omega_0}{\sigma}}) ]
#'     - 4 ~ p^3 ~ (1-p) ~ \Phi({\frac{\nu - H_t}{\sigma}}))
#' \end{split}
#' \end{equation}
#'


# plot bias from true- at t
gTN_simple = 
    expand_grid(
        delta_t = c(0.25, 0.5, 1)
        , lambda_t = 0.03
        , nu = c(1,2)
        , sigma = Sigma
        , omega0 = Omega0_rel * Omega
    ) %>%
    mutate(
        product = delta_t * lambda_t
        , p = exp(-product)
        , Bias.lower = p^4 * (1 - pnorm(nu, omega0, sigma)) - 4 * p^3 * (1-p) * pnorm(nu, Omega + Gamma * exp(-Delta * delta_t), sigma)
        , Bias.upper = p^4 * (1 - pnorm(nu, omega0, sigma)) - 4 * p^3 * (1-p) * pnorm(nu, Omega + Gamma, sigma)
        , delta_t = factor(delta_t)
        , nu = factor(nu)
        , omega0 = factor(omega0)
    ) %>%
    ggplot(aes(x = delta_t, color = nu))+
    geom_hline(yintercept = 0, linetype = 3)+
    geom_errorbar(aes(
            ymin = Bias.lower, ymax = Bias.upper
            , group = paste(nu, omega0, delta_t)
        )
        , width = 0.3
        , linewidth = 1.5
        , position = position_dodge(0.5)
    )+
    xlab("Interval length (yrs)")+
    ylab("Bias if TN at t")+
    scale_color_brewer("Cut-off", palette = "Dark2", guide = "legend")



#' For the first case in Equation (\ref{eq:bias_tneg}) where no infection occurred, the bias increases when the seropostivity threshold is lowered.
#' In the second case where only one infection occurred, the bias becomes less negative when time since infection is small or when the seropostivity threshold is lowered.
#' Since time since infection in this case is bounded by the interval length, shorter bleeding intervals are more likely to be less negative.
#' Titers in individuals who have acquired multiple DENV infections, i.e., multitypic titers, are more robust and can be assumed to test negative at negligible levels. It follows that minimal bias would arise from the last case where more than one infection occurred during the interval.
#' 




#' \textbf{Biases in seroconversion probabilities when falsely negative at t},
#' under these described processes, 
#' can be expressed as
#  ..................

#' \begin{equation} \label{eq:bias_fneg}
#' \begin{split}
#' & P( \oplus^T_{(t, t + \Delta{t}]} | \ominus^F_{t}) + P( \oplus^F_{(t, t + \Delta{t}]} | \ominus^F_{t} ) 
#'     - P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) - P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) \\
#' & = P( \oplus^T_{(t, t + \Delta{t}]} | \ominus_{t}, S_{t} = 1) + 0
#'     - P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) - P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) \\
#' & = P( \oplus_{(t, t + \Delta{t}]}, S_{(t, t + \Delta{t}]} = 1 | \ominus_{t}, S_{t} = 1)
#'     + P( \oplus^T_{(t, t + \Delta{t}]}, S_{(t, t + \Delta{t}]} > 1 | \ominus^F_{t}, S_{t} = 1) \\
#' &   - P(S_{(t, t + \Delta{t}]} = 1 | S_{t} = 0) - P(S_{(t, t + \Delta{t}]} > 1 | S_{t} = 0) \\
#' \end{split}
#' \end{equation}
#'


#' The first case where $S_{(t, t + \Delta{t}]} = S_{t} = 1$ means that no infection occurred during the interval. $P(\oplus_{(t, t + \Delta{t}]})$ in this case has a lowerbound of $1 - \Phi(\frac{\nu - \omega}{\sigma})$. 
#' In the second case where infection(s) occurred during the interval, $S_{t} = 1$ and $S_{(t, t + \Delta{t}]} > 1$, $P(\oplus_{(t, t + \Delta{t}]}) \approx 1$.
#' Hence, $P(\oplus^T_{(t, t + \Delta{t}]} | \ominus^F_{t})$ is a weighted average between two quantities, $1 - \Phi(\frac{\nu - \omega}{\sigma})$ and $1$. Longer time between blood draws, $\Delta{t}$, and higher FOI, $\lambda_t$, during the interval tends the quantity towards $1$.
#' Again, for probability of escaping a particular serotype during the interval $p = e^{-\lambda_t \Delta{t}}$,
#' we can express (\ref{eq:bias_fneg}) as

#' \begin{equation} \label{eq:bias_fneg_p}
#' \begin{split}
#' & p^3 * (1 - \Phi(\frac{\nu - H_t}{\sigma})) + (1 - p^3) - (1 - p^4) \\
#' & = p^4 - \Phi(\frac{\nu - H_t}{\sigma}) p^3
#' \end{split}
#' \end{equation}
#'

# plot bias from false- at t
gFN_simple = 
    expand_grid(
        delta_t = c(0.25, 0.5, 1)
        , lambda_t = 0.03
        , nu = c(1,2)
    ) %>%
    mutate(
        product = delta_t * lambda_t
        , p = exp(-product)
        , Bias.lower = p^4 - pnorm(nu, Omega + Gamma * exp(-Delta * delta_t), Sigma) * p^3
        , Bias.upper = p^4 - pnorm(nu, Omega + Gamma, Sigma) * p^3
        , nu = factor(nu)
        , delta_t = factor(delta_t)
    ) %>%
    ggplot(aes(x = delta_t, color = nu))+
    geom_hline(yintercept = 0, linetype = 3)+
    geom_errorbar(aes(
            ymin = Bias.lower, ymax = Bias.upper
            , group = paste(nu, delta_t)
        )
        , width = 0.3
        , linewidth = 1.5
        , position = position_dodge(0.5)
    )+
    xlab("Interval length (yrs)")+
    ylab("Bias if FN at t")+
    scale_color_brewer("Cut-off", palette = "Dark2", guide = "legend")


#' \textbf{Relative contributions of the biases} can be expressed as
#  ..................

#' \begin{equation} \label{eq:rel}
#' \begin{split}
#' \frac{P(\ominus^F_{t}) }{ P(\ominus^T_{t})}
#' & = \frac{ P(S_{t} = 1) ~ P(\ominus^F_{t} | S_{t} = 1) }{ P(S_{t} = 0) ~ P(\ominus^T_{t} | S_{t} = 0) } \\
#' & = \frac{ {4 \choose 1} ~ (e^{-\bar{\lambda}a})^3 ~ (1 - e^{-\bar{\lambda}a}) }{ (e^{-\bar{\lambda}a})^4 } ~ \frac{ P(\ominus^F_{t} | S_{t} = 1) }{ P(\ominus^T_{t} | S_{t} = 0) } \\
#' & = 4 ~ \frac{ (1 - e^{-\bar{\lambda}a}) }{ (e^{-\bar{\lambda}a}) } ~ \frac{ P(\ominus^F_{t} | S_{t} = 1) }{ P(\ominus^T_{t} | S_{t} = 0) } \\
#' \end{split}
#' \end{equation}


#' The quantity is a product of the constant 4 and two variable components. The first variable component is zero when $a$ is zero or the average per-serotype force of infection $\bar{\lambda}$ is zero and grows as the product $\bar{\lambda}a$ increases.
#' The second variable component depends on the recency of the dengue infection with bounds
#'

#' \begin{equation}
#'     \frac{ \Phi(\frac{\nu - \omega - \gamma}{\sigma}) }{ \Phi(\frac{\nu - \omega_0}{\sigma}) }
#' \le \frac{ P(\ominus^F_{t} | S_{t} = 1) }{ P(\ominus^T_{t} | S_{t} = 0) } 
#' \le \frac{ \Phi(\frac{\nu - \omega}{\sigma}) }{ \Phi(\frac{\nu - \omega_0}{\sigma}) }
#' \le 1
#' \end{equation}
#' 


#' Assuming cross-reactive titer $\omega_0$ is at most $\omega$, the quantity is bounded at one.
#' Intuitively, lowering the positivity cut-off $\nu$ would decrease contribution of $\ominus^F_{t}$ to $\ominus_{t}$. However, the contribution can only be minimized down to $1 / (1 + { \Phi(\frac{\nu - \omega_0}{\sigma}) / \Phi(\frac{\nu - \omega - \gamma}{\sigma}) })$ and the efficiency of such reduction declines with age and historical FOI $\bar{\lambda}$.
#'

gRelative_simple = 
    expand_grid(
        nu = c(1,2)
        , avgLambda = 0.03
        , age = seq(5,30, by = 5)
    ) %>%
    mutate(val = (1 - exp(-avgLambda * age)) / exp(-avgLambda * age)) %>%    
    mutate(Upper = pnorm((nu - Omega)/Sigma) / pnorm((nu)/Sigma)) %>%
    mutate(Lower = pnorm((nu - Omega - Gamma)/Sigma) / pnorm((nu)/Sigma)) %>%
    mutate(
        age = factor(age)
        , nu = factor(nu)
    ) %>%
    ggplot(aes(x = age, fill = nu))+
    geom_col(aes(y = 4 * val * Upper, group = paste(nu, age))
        , position = "dodge"
    )+
    scale_y_continuous(
        "Relative contribution\n(FN:TN)", trans = "log2"
    )+
    xlab("Age (yrs)")+
    scale_fill_brewer("Cut-off"
        , palette = "Dark2"
        , guide = "legend"
        , labels = c(10,20)
    )

    
# .... plot the biases

#+ fig.height = 6
#+ fig.cap="\\label{fig:demo} Biases in seroconversion probabilities illustrated using antibody kinetics estimates from Salje, 2018. Biases in the case where the individual was a) true negative (TN) at pre-interval bleed versus b) false negative (FN) at pre-interval bleed. c) Upper bound of relative contribution between the FN and TN case to the seroconversion probability as a function of pre-interval age for seropositivity cut-offs 10 and 20."
ggarrange(
    ggarrange(
        gTN_simple + coord_cartesian(ylim = c(-0.2, 1))
        , gFN_simple + coord_cartesian(ylim = c(-0.2, 1))
        , nrow = 1
        , ncol = 2
        , common.legend = T
        , legend = "none"
        , align = 'h'
        , labels = "auto"
    )
    , gRelative_simple
    , ncol = 2
    , nrow = 1
    , widths = c(2,1.5)
    , labels = c("", "c")
)

#' Figure \ref{fig:demo} illustrates ranges of the biases and their relative contributions using mean parameter estimates of $\omega=`r Omega`$ and $\sigma=`r Sigma`$ from Salje, 2018 under various interval lengths and seropositivity cut-offs.
#' The amount of cross-reactive titers relative to DENV-specific titers was artibrarily set to 20\% ($\omega_0 = 0.2 \omega$). Per-serotype FOI (historical and during the interval) is set to 0.03.

#' Under these parameters, we can see that the magnitude of positive bias from false negative scenario at $t$ is far higher biases from the true negative scenario. Increasing the cut-off reduces the lower bound of such inflation (leading to greater uncertainty in the amount of bias introduced) and increases the contribution of the false negative scenario.

    
    



# ................

#' \textbf{Effects of seasonality on test positivity.}
#' The derivations, so far, assumes that FOI is constant throughout the year. 
#' However, it is conceivable that timings of the bleeds in relation to phases in the dengue seasons would lead to non-uniform distribution of the infections, and hence, may influence the expected biases.
#'

#' Let the per-serotype FOI at age $a$ of an individual $\lambda(a) = \bar{\lambda} + \delta_\lambda cos(2 \pi (a + \delta_t))$ where
#' $\delta_\lambda$ is the magnitude of fluctuation and $\delta_t$ governs the phase of the fluctuation in relation to age.
#' The per-serotype infection risk faced up till age $A$ is then
#'

#' \begin{equation}
#' \Lambda(A)
#' = \int_0^A \bar{\lambda} + \delta_\lambda cos(2 \pi (a + \delta_t)) \,da
#' = \bar{\lambda}A + \frac{\delta_\lambda}{2\pi} [ sin(2 \pi (A + \delta_t)) - sin(2 \pi \delta_t) ]
#' \end{equation}
#'

#' We can write the probability that the individual was infected by a particular serotype at age $A$ as $e^{-\Lambda(A)}\lambda(A)$.
#' Thus, the expected underlying titer of the individual at age $A_1$ given the individual was monotypically infected at that point is

#' \begin{equation}
#' \begin{split}
#' \mathbb{E}[H_{A_1} | S_{A_1} = 1]
#' & = \frac{3~e^{-\Lambda(A_1)} ~ \int_0^{A_1} (\omega + \gamma e^{-\delta(A_1 - A)}) ~ e^{-\Lambda(A)}\lambda(A) \,dA}{3~e^{-\Lambda(A_1)}~(1 - e^{-\Lambda(A_1)})} \\
#' & = \frac{\int_0^{A_1} (\omega + \gamma e^{-\delta(A_1 - A)}) ~ e^{-\Lambda(A)}\lambda(A) \,dA}{1 - e^{-\Lambda(A_1)}}
#' \end{split}
#' \end{equation}
#'

#' The growing contribution of $\bar{\lambda}A$ in $\Lambda(A)$ as age increases shrinks the impact of timing of the pre-interval blood draw in relation to seasonality.
#' Figure \ref{fig:season}a-c illustrates this shrinkage under parameter estimates from Salje, 2018.
#' On the other hand, for infections that occurred within the intervals, seasonality and the interval lengths dictates the distribution of time since infection at the second bleed (Figure \ref{fig:season}d). These timing variations affects the amount of titers waned, and consequentially, the probability of falsely testing negative at post-interval bleed. Figure \ref{fig:season}e illustrates this nuance under parameter estimates from Salje, 2018.


# FOI at time t
foi.t = function(time, phaseDiff, avgFoi = 0.03, deltaFoi = 0.02){
    # phaseDiff:
    # 1   = start at trough
    # 0.5 = start at midpoint between trough and peak
    # 0   = start at peak
    avgFoi + deltaFoi * cos((time -0.5 + phaseDiff) * 2*pi)
}

# FOI faced up till time t
foi.cum.t = function(time, phaseDiff, avgFoi = 0.03, deltaFoi = 0.02){
    # phaseDiff:
    # 1   = start at trough
    # 0.5 = start at midpoint between trough and peak
    # 0   = start at peak
    avgFoi*time + deltaFoi/(2*pi) * (
        sin( (time -0.5 + phaseDiff) * 2*pi ) -
        sin( (     -0.5 + phaseDiff) * 2*pi )
    )
}

# Probability of infection occurring at time t
probTime = function(time, phaseDiff, avgFoi = 0.03, deltaFoi = 0.02){
    # phaseDiff:
    # 0    = start at trough
    # 0.25 = start at midpoint between trough and peak (1/4 of season)
    # 0.5  = start at peak (1/2 of season)
    
    # probability of infection happening at time t
    # = probability of escaping up till t - dt *
    exp(-foi.cum.t(time, phaseDiff, avgFoi, deltaFoi)) *
    #   probability of being infected within dt
    foi.t(time, phaseDiff, avgFoi, deltaFoi)
}
    # # --- check that the math matches up -----
    # 1 - exp(-foi.cum.t(10.2, 0.5))
    # integrate(probTime, 0, 10.2, phaseDiff = 0.5)
    # # ----------------------------------------

# Probability of being monotypically infected at time t
probMono = function(time, phaseDiff, avgFoi = 0.03, deltaFoi = 0.02){
    # probability of escaping a particular serotype up till time t
    p = exp(-foi.cum.t(time, phaseDiff, avgFoi, deltaFoi))
    choose(4,1) * p^3 * (1-p)
}

# Expected titer at time t given infected once
monoTiter = function(times.eval, ...){
    sapply(times.eval, function(time.eval){
        integrate( function(time.infect){
                probTime(time.infect, ...) *
                (Omega + Gamma * exp(-Delta * (time.eval - time.infect)))
            }
            , 0
            , time.eval
        )$value /
        (1 - exp(-foi.cum.t(time.eval, ...)))
    })
    
}


Phase = c(
    '0' = 'Start'
    , '0.25' = 'Rise'
    , '0.5' = 'Peak'
)

Ages = seq(0, 10, by = 0.05)
gFoi =
    expand_grid(
        Age = Ages
        , phaseDiff = c(0, 0.25, 0.5)
    ) %>%
    mutate(
        FOI = foi.t(Age, phaseDiff)
        , Phase = factor(Phase[as.character(phaseDiff)], levels = Phase)
    ) %>%
    ggplot(aes(x = Age, y = FOI, group = Phase, color = Phase))+
    geom_line()

gTiter =
    expand_grid(
        Age = Ages
        , phaseDiff = c(0, 0.25, 0.5)
        , nu = c(1,2)
    ) %>%
    mutate(
        Titer = mapply(monoTiter, times.eval = Age, phaseDiff = phaseDiff)
        , pFN = pnorm((nu - Titer)/Sigma) *
            # scale for plotting
            (Omega + Gamma)
        , Phase = factor(Phase[as.character(phaseDiff)], levels = Phase)
    ) %>%
    ggplot(aes(x = Age, y = Titer, group = Phase, color = Phase))+
    geom_line(alpha = 0.3, linewidth = 1.2)+
    geom_line(aes(y = pFN, group = paste(Phase, nu), linetype = factor(nu)))+
    scale_y_continuous(
        'Expected titer\nat 1st bleed'
        , sec.axis = sec_axis(
            ~ . / (Omega + Gamma)
            , name = "Prob. FN"
        )
    )+
    scale_linetype_manual('Cut-off'
        , values = c(1,2)
        , labels = c(10,20)
    )+
    theme(
        axis.title.y.left = element_text(color = 'grey')
        , axis.text.y.left = element_text(color = 'grey')
        , axis.line.y.left = element_line(color = 'grey')
        , axis.ticks.y.left = element_line(color = 'grey')
        , legend.box = 'horizontal'
    )

gMono =
    expand_grid(
        Age = Ages
        , phaseDiff = c(0, 0.25, 0.5)
    ) %>%
    mutate(
        pMono = mapply(probMono, time = Age, phaseDiff = phaseDiff)
        , Phase = factor(Phase[as.character(phaseDiff)], levels = Phase)
    ) %>%
    ggplot(aes(x = Age, y = pMono, group = Phase, color = Phase))+
    geom_line()+
    ylab('Prob. mono')
    
gFoi.interval =
    expand_grid(
        delta_t = seq(0,1,by=0.01)
        , phaseDiff = c(0, 0.25, 0.5)
    ) %>%
    mutate(
        FOI = foi.t(delta_t, phaseDiff)
        , Phase = factor(Phase[as.character(phaseDiff)], levels = Phase)
    ) %>%
    ggplot(aes(x = delta_t, y = FOI, group = Phase, color = Phase))+
    geom_line()
    
gTiter.interval =
    expand_grid(
        delta_t = seq(0,1,by=0.01)
        , phaseDiff = c(0, 0.25, 0.5)
        , nu = c(1,2)
    ) %>%
    mutate(
        Titer = mapply(monoTiter, times.eval = delta_t, phaseDiff = phaseDiff)
        , pFN = pnorm((nu - Titer)/Sigma) *
            # scale for plotting
            (Omega + Gamma)
        , Phase = factor(Phase[as.character(phaseDiff)], levels = Phase)
    ) %>%
    ggplot(aes(x = delta_t, y = Titer, group = Phase, color = Phase))+
    geom_line(alpha = 0.3, linewidth = 1.2)+
    geom_line(aes(y = pFN, group = paste(Phase, nu), linetype = factor(nu)))+
    scale_y_continuous(
        'Expected titer\nat 2nd bleed'
        , sec.axis = sec_axis(
            ~ . / (Omega + Gamma)
            , name = "Prob. FN"
        )
    )+
    xlab('Interval length (yrs)')+
    theme(
        axis.title.y.left = element_text(color = 'grey')
        , axis.text.y.left = element_text(color = 'grey')
        , axis.line.y.left = element_line(color = 'grey')
        , axis.ticks.y.left = element_line(color = 'grey')
    )



#+ fig.cap="\\label{fig:season} Effects of timings of bleeds in relation to seasonality illustrated using antibody kinetics estimates from Salje, 2018. a) Per-serotype seasonal force of infection (FOI) in age. b) Expected pre-interval titer of a monotypically infected individual (tinted) and the corresponding probability of falsely testing negative given the titers at positivity threshold of 10 (solid) and 20 (dashed). c) Probability that the individual was monotypically infected. d) Relationship between per-serotype seasonal force of infection (FOI) and years since pre-interval bleed. e) Expected post-interval titer of individuals that acquired their first infection during the interval (tinted) and the corresponding probability of falsely testing negative given the titers at positivity threshold of 10 (solid) and 20 (dashed)."
ggarrange(
    gFoi
    , gFoi.interval
    , gTiter
    , gTiter.interval
    , gMono
    # , ggplot()+theme_void()
    , as_ggplot(get_legend(gTiter))
    , ncol = 2
    , nrow = 3
    , heights = c(1.5, 3, 2)
    , widths = c(5,3.5)
    , labels = c('a','d','b','e','c','')
    , hjust = 0
    , align = 'hv'
    , common.legend = T
    , legend = "none"
)

