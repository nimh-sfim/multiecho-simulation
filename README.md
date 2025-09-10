# Tools to simulate multi-echo fMRI time series

The core approach of this simulation is:  
The user inputs time series, baseline parameters including echo time,
and the proportion of $\Delta R_2^*$ and $\Delta S_0$.
The code can then generate time series for echo times given the baseline parameters.
There can be many sources of the intial time series (pure simulation or real data),
and the resulting multi-echo data will be based on a ground-truth calculation
given the inputted proportions of $\Delta R_2^*$ and $\Delta S_0$.
That is, if signal fluctuations are purely $\Delta S_0$,
then the % change from mean time series will be identical at all echo times.
If signal fluctuations are purely $\Delta R_2^*$,
then the amplitude of fluctuations will scale with echo time.
Proportions in-between will show a mix.

This also includes some functions to help generate the intial time series.


## Underlying math for the simulation

$S = (\overline{S_0} + \Delta S_0) * e^{-TE * (\overline{R_2^*} + \Delta R_2^*)}$  
$X = (\overline{S_0} + \Delta S_0)$  
$Y = e^{-TE * (\overline{R_2^*} + \Delta R_2^*)}$  
$S = X*Y$  
$var(S) = var(XY)$   
$E[X] = mean(X)=\overline{S_0}$  
$E[Y] = mean(Y) = e^{-TE * \overline{R_2^*}}$   
$var(X) = E[(\Delta S_0)^2]$  
$var(Y) = E[(e^{-TE * (\overline{R_2^*} + \Delta R_2^*)} - e^{-TE * (\Delta R_2^*)})^2]=E[(e^{-TE * (\Delta R_2^*)}*(e^{-TE * \overline{R_2^*}}-1))^2]$  
Using identitiy  
$var(XY) = (E[X])^2*var(Y) + (E[Y])^2*var(X) + var(X)*var(Y)$


$prop(\frac{\Delta S_0 }{ \Delta R_2^*})$ = $p_{S_0}$ Proportion of $\Delta S_0$ to $\Delta R_2^*$.  
1 = pure $\Delta S_0$ and 0= pure $\Delta R_2^*$

By the initial specificiation of the goal:  
$C = constant$  
$p_{S_0}*var(X) + (1-p_{S_0})*var(Y) = C$
The above cannot be correct because $p_{S_0}*var(X)*(1-p_{S_0})*var(Y)$ would inherantly be constant and $(E[X])^2*p_{S_0}*var(Y) + (E[Y])^2*(1-p_{S_0})*var(X)$ could not also be a constant.


I think we'd need to solve for
$\frac{p_{S_0}*var(X)}{(E[Y])^2} + \frac{(1-p_{S_0})*var(Y)}{(E[X])^2}=C$  
The $var(X)*var(Y)$ term would still be constant, but this would then be solvable. (need to think more if this is breaking the specification I'm tying to solve for)

Identity if the above is a constant  
$p*A + (1-p)*B = (1-p)*A + pB$  
$p*A +B -pB = A - pA +pB$  
$2p*A-A = 2p*B-B$  
$A*(2p-1) = B*(2p-1)$  
Cancelling out p and claiming unequal things are equal so I'm doing something wrong.
 
That is, the variance explained by X and Y
$(1-p_{S_0})*(E[X])^2*var(Y) + p_{S_0}*(E[Y])^2*var(X) = C$


older stuff

$prop(\frac{\Delta S_0 }{ \Delta R_2^*})$ = $p_{S_0}$ Proportion of $\Delta S_0$ to $\Delta R_2^*$.  
1 = pure $\Delta S_0$ and 0= pure $\Delta R_2^*$

$TE = echo\space  time$

$S = (S_0 + p_{S_0} * \Delta S_0) * e^{-TE * (R_2^* + (1 - p_{S_0}) * \Delta R_2^*)}$

If $S$, a baseline $S_0$, a baseline $R_2^*$, a baseline $TE$ and $p_{S_0}$ are provided, there should be stable solutions for $\Delta R_2^*$ and $\Delta S_0$.
Thus, new $S$ time series can be calculated using those $\Delta R_2^*$ and $\Delta S_0$ values and a different $TE$

$ln(S_0 + p_{S_0} * \Delta S_0) - ln(S)  = TE * R_2^* + TE * ((1 - p_{S_0}) * \Delta R_2^*)$

$(1 - p_{S_0}) * \Delta R_2^* = \frac{ln(S_0 + p_{S_0} * \Delta S_0) - ln(S)}{TE} - R_2^*$

Since only want to solve for $0 \leq p_{S_0} \leq1$ first constraining the above equation by setting $p_{S_0}=1$

$0 = \frac{ln(S_0 + \Delta S_0) - ln(S)}{TE} - R_2^*$

$TE * R_2^*  = ln(\frac{S_0 + \Delta S_0}{S})$

$\frac{S_0 + \Delta S_0}{S} = e^{TE * R_2^*}$

$\Delta S_0 = S*e^{TE * R_2^*} - S_0$

With a known $\Delta S_0$, the inputted $p_{S_0}$ can be used to solve for $R_2^*$ for $p_{S_0}\lt1$

$(1 - p_{S_0}) * \Delta R_2^* = \frac{ln(S_0 + p_{S_0} * \Delta S_0) - ln(S)}{TE} - R_2^*$

$\Delta R_2^* = \frac{\frac{ln(S_0 + p_{S_0} * \Delta S_0) - ln(S)}{TE} - R_2^*}{1 - p_{S_0}}$

Replacing ${\Delta S_0}$:  
$ln(S_0 + p_{S_0} * \Delta S_0) - ln(S) =$
$ln(S_0 + p_{S_0} * S*e^{TE * R_2^*} - S_0) - ln(S) =$
$ln(p_{S_0} * S*e^{TE * R_2^*}) - ln(S) =$
$ln(\frac{p_{S_0} * S*e^{TE * R_2^*}}{S}) =$
$ln(p_{S_0} * e^{TE * R_2^*})= $
$ln(p_{S_0}) +  TE * R_2^*$

$\Delta R_2^* = \frac{\frac{ln(p_{S_0}) +  TE * R_2^*}{TE} - R_2^*}{1 - p_{S_0}}= $
$\frac{\frac{ln(p_{S_0})}{TE}}{1 - p_{S_0}}= $
$\frac{(1 - p_{S_0})* ln(p_{S_0})}{TE}$


Working here & trying to resolve some issues:

    For proportion_s0_r2s==1:
        0 = (np.log(s0_baseline + delta_s0_scale) - te*r2s_baseline - np.log(S)) / TE
        np.log(s0_baseline + delta_s0_scale) = te*r2s_baseline + np.log(S)
        delta_s0_scale = np.exp(te*r2s_baseline + np.log(S)) - s0_baseline
        The above eq is used to calculate delta_s0_scale in the function
    With a known delta_s0_scale, the above eq with proportion_s0_r2s can be used to solve for delta_r2s_scale
        (1-proportion_s0_r2s)*delta_r2s_scale = (np.log(s0_baseline + proportion_s0_r2s*delta_s0_scale) - te*r2s_baseline- np.log(S)) / TE
        delta_r2s_scale = (np.log(s0_baseline+proportion_s0_r2s*delta_s0_scale) - te*r2s_baseline - np.log(S))/((1-proportion_s0_r2s)*te)
        The above eq is used to calculate delta_r2s_scale in the function, but it fails for 1-proportion_s0_r2s==0, but,
        in that case (1-proportion_s0_r2s)*delta_r2s_scale should be 0 so it's just set to 0
    The returned delta_s0 and delta_r2s values are the scaled values multiplied by their proportions