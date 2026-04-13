As we were working through the simulation, we attempted may different formations of the core simulation equations.
Most of these tried to solve the system of equations without any approximations.
The following is effectively my scrap notes from those attempts.
They are (temporarily) worth backin up, but aren't actually used for anything.


Starting with the assumption that this proportional relationship can be model with a linear scaling of the two $\Delta$ terms:

$S(t,TE) = (\overline{S_0(TE)} + p*\Delta S_0(t,TE)) * e^{-TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))}$

That is, when p=1, temporal fluctuations of $S(t,TE)$ are fully modeled by $\Delta S_0(t,TE)$,
when p=0, temporal fluctuations of $S(t,TE)$ are fully modeled by $\Delta R_2^*(t,TE)$,
and intermediate values of p represent a mixes of changes from each.

Constraint: for a given $S(t,TE)$, $\overline{S_0(TE)}$, $\overline{R_2^*(TE)}$, & $TE$:

$p \in [0, 1]$

$(abs(\Delta S_0(t,TE))\ when\ p=1)\gt(abs(\Delta S_0(t,TE))\ when\ p\lt1)$

$(abs(\Delta R_2^*(t,TE))\ when\ p=0)\gt(abs(\Delta S_0(t,TE))\ when\ p\gt0)$

There are an infinite number of solutions without this constraint.
This constraint is saying that $\Delta S_0(t,TE))$ should be biggest when fully modeled by just that source.

$S(t,TE) = (\overline{S_0(TE)} + p*\Delta S_0(t,TE)) * e^{-TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))}$

$(\overline{S_0(TE)} + p*\Delta S_0(t,TE)) = S(t,TE) * e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))}$

$\Delta S_0(t,TE) = \frac{S(t,TE) * e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} - \overline{S_0(TE)}}{p}$

When $p=1$

$\Delta S_0(t,TE) = S(t,TE) * e^{TE * \overline{R_2^*(TE)}} - \overline{S_0(TE)}$

For $p\lt1$

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{\overline{S_0(TE)} + p*\Delta S_0(t,TE)}{S(t,TE)}$

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{\overline{S_0(TE)} + p*(S(t,TE) * e^{TE * \overline{R_2^*(TE)}} - \overline{S_0(TE)})}{S(t,TE)} $

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{\overline{S_0(TE)} + p*S(t,TE) * e^{TE * \overline{R_2^*(TE)}} - p*\overline{S_0(TE)}}{S(t,TE)} $

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{(1-p)*\overline{S_0(TE)} + p*S(t,TE) * e^{TE * \overline{R_2^*(TE)}}}{S(t,TE)}$

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)} + p * e^{TE * \overline{R_2^*(TE)}}$

$e^{TE * \overline{R_2^*(TE)}}*e^{(1-p)*\Delta R_2^*(t,TE)} = \frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)} + p * e^{TE * \overline{R_2^*(TE)}}$

$e^{TE * \overline{R_2^*(TE)}}*(e^{(1-p)*\Delta R_2^*(t,TE)} - p) = \frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)}$

$e^{(1-p)*\Delta R_2^*(t,TE)} = \frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)}*e^{-TE * \overline{R_2^*(TE)}}+p$

$(1-p)*\Delta R_2^*(t,TE) = ln(\frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)}*e^{-TE * \overline{R_2^*(TE)}}+p)$

$\Delta R_2^*(t,TE) = \frac{ln(\frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)}*e^{-TE * \overline{R_2^*(TE)}}+p)}{1-p}$

Sanity check of solving for $p=0$ from original equation and then from full solution for $\Delta R_2^*(t,TE)$

$e^{TE * (\overline{R_2^*(TE)} + (1-p)*\Delta R_2^*(t,TE))} = \frac{\overline{S_0(TE)} + p*\Delta S_0(t,TE)}{S(t,TE)}$ Original eq

$e^{TE * (\overline{R_2^*(TE)} + \Delta R_2^*(t,TE))} = \frac{\overline{S_0(TE)}}{S(t,TE)}$

$TE * (\overline{R_2^*(TE)} + \Delta R_2^*(t,TE)) = ln(\frac{\overline{S_0(TE)}}{S(t,TE)})$

$\bm{\Delta R_2^*(t,TE) = \frac{ln(\frac{\overline{S_0(TE)}}{S(t,TE)})}{TE} - \overline{R_2^*(TE)}}$

$\Delta R_2^*(t,TE) = \frac{ln(\frac{(1-p)*\overline{S_0(TE)}}{S(t,TE)}*e^{-TE * \overline{R_2^*(TE)}}+p)}{1-p}$ for $p=0$

$\Delta R_2^*(t,TE) = ln(\frac{\overline{S_0(TE)}}{S(t,TE)}*e^{-TE * \overline{R_2^*(TE)}})$

$\bm{\Delta R_2^*(t,TE) = ln(\frac{\overline{S_0(TE)}}{S(t,TE)}) -TE * \overline{R_2^*(TE)}}$

***
***
***
***
***
***
***
***
***
***
***
***

$S = (\overline{S_0} + \Delta S_0) * e^{-TE * (\overline{R_2^*} + \Delta R_2^*)}$

Not included to keep equations simple, but $S$
can have a different value for each voxel (sample) & timepoint and is TE sensitive.
$\overline{S_0}\ \&\ \overline{R_2^*}$ are the means across time so they can have a distinct value for each voxel.
$\Delta S_0, \&\ \Delta R_2^*$ have different values for each voxel and timepoint and are NOT TE sensitive.

Assuming $\Delta R_2^*=0$: 

$S = (\overline{S_0} + \Delta S_{0_{full}}) * e^{-TE * \overline{R_2^*}}$

$\Delta S_{0_{full}} = \frac{S}{e^{-TE * \overline{R_2^*}}} - \overline{S_0}$

$\Delta S_{0_{full}}$ is the maximum possible values of $\Delta S_0$ for a given $S, \overline{S_0}, \overline{R_2^*}, \&\ TE$.
That is, it is if the changes across time in $S$ are fully modeled by $\Delta S_0$

The $\Delta S_{0_{full}}$ can then be scaled from $proportion=[0, 1]$.

**NOTE**: $proportion$ is **not** a division.
$\Delta S_0$ and $\Delta R_2^*$ have different ranges and the scaling to go from fully one to fully the other is not equal.

$\Delta S_{0_{prop}} = proportion * \Delta S_{0_{full}}$

$S = (\overline{S_0} + \Delta S_{0_{prop}}) * e^{-TE * (\overline{R_2^*} + \Delta R_{2_{prop}}^*)}$

$\Delta R_{2_{prop}}^* = \frac{-ln(\frac{S}{\overline{S_0} + \Delta S_{0_{prop}}})}{TE} - \overline{R_2^*}$


Obseration is: $ \overline{R_2^*} + \Delta R_{2_{prop}}^* = A*(\overline{S_0} + \Delta S_{0_{prop}}) + B$, where $A \lt 0$

***
***
***


$\frac{(\overline{S_0} + \Delta S_{0_{prop}})}{S} =  e^{TE * (\overline{R_2^*} + \Delta R_{2_{prop}}^*)}$

$\ln{(\overline{S_0} + \Delta S_{0_{prop}})} - ln{S} =  TE * (\overline{R_2^*} + \Delta R_{2_{prop}}^*)$

$\ln{(\overline{S_0} + proportion * \Delta S_{0_{full}})} - ln{S} =  TE * (\overline{R_2^*} + \Delta R_{2_{prop}}^*)$

$\ln{(\overline{S_0} + proportion * (\frac{S}{e^{-TE * \overline{R_2^*}}} - \overline{S_0}))} - ln{S} =  TE * (\overline{R_2^*} + \Delta R_{2_{prop}}^*)$



***
***
***

$\Delta R_{2_{prop}}^* = \frac{-ln(S) + ln(\overline{S_0} + proportion * \Delta S_{0_{full}})}{TE} - \overline{R_2^*}$

$\Delta R_{2_{prop}}^* = \frac{ ln(\overline{S_0} + proportion * \Delta S_{0_{full}})}{TE} - \overline{R_2^*} - \frac{ln(S)}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ ln(\overline{S_0} + proportion * (\frac{S}{e^{-TE * \overline{R_2^*}}} - \overline{S_0}))}{TE} - \overline{R_2^*} - \frac{ln(S)}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\overline{S_0}*(1-proportion) + proportion * \frac{S}{e^{-TE * \overline{R_2^*}}})}{TE} - \overline{R_2^*} - \frac{ln(S)}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\overline{S_0}*(1-proportion) + proportion * \frac{(\overline{S_0} + \Delta S_0) * e^{-TE * (\overline{R_2^*} + \Delta R_2^*)}}{e^{-TE * \overline{R_2^*}}})}{TE} - \overline{R_2^*} - \frac{ln(S)}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\overline{S_0}*(1-proportion) + proportion * (\overline{S_0} + \Delta S_0) * e^{-TE * \Delta R_2^*})}{TE} - \overline{R_2^*} - \frac{ln(S)}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\overline{S_0}*(1-proportion) + proportion * (\overline{S_0} + \Delta S_0) * e^{-TE * \Delta R_2^*})}{TE} - \overline{R_2^*} - \frac{ln((\overline{S_0} + \Delta S_{0_{full}}) * e^{-TE * \overline{R_2^*}})}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\overline{S_0}*(1-proportion) + proportion * (\overline{S_0} + \Delta S_0) * e^{-TE * \Delta R_2^*})}{TE} - \overline{R_2^*} - \frac{ln(\overline{S_0} + \Delta S_{0_{full}})}{TE} + \frac{-TE*\overline{R_2^*}}{TE}$




looks either like a mess or I'm heading to $\Delta R_{2}^* = \Delta R_{2}^*$

$\Delta R_{2_{prop}}^* = \frac{ln(proportion * S)}{TE} - \frac{ln(S)}{TE} - \overline{R_2^*} - \frac{ln(e^{-TE * \overline{R_2^*}})}{TE}$

$\Delta R_{2_{prop}}^* = \frac{ln(\frac{proportion * S}{S})}{TE} - \overline{R_2^*} + \overline{R_2^*}$
    
$\Delta R_{2_{prop}}^* = \frac{ln(proportion)}{TE}$


The initial $S$ can then be generated with a user-defined $TE, \overline{S_0}\ \&\ \overline{R_2^*}$
and, for a user-defined $\frac{\Delta S_0}{\Delta R_2^*} proportion$
a new $S$ can be calculated for a different $TE$

***
***
***
***
***
***
***
***
Everything below here was messing around & may eventually be deleted.

The $\Delta S_0$ time series
Can linearly scale the proportion of $\Delta S_0$ from 1 to 0 and calculate $\Delta R_2^*$




$X = (\overline{S_0} + \Delta S_0)$  
$Y = e^{-TE * (\overline{R_2^*} + \Delta R_2^*)}$  
$S = X*Y$  
Assuming independent variables (a bad assumption, but an approximation)
$var(XY) = E[X^2]*E[Y^2] - E[X]^2*E[Y]^2$

$E[X] = mean(X)=\overline{S_0}$ (a predefined constant)

$E[Y] = mean(Y) = e^{-TE * \overline{R_2^*}}$ (a predefined constant)

What happens if $p*X^2 = (1-p)*Y^2$ ?

******
NEED TO TAKE A STEP BACK AND FIGURE OUT WHAT I'M TRYING TO CALCULATE.
$\Delta R_2^*$ and $Delta S_0$ ARE COMPLETELY INTERCONNECTED BECAUSE A LARGER VALUE OF ONE WILL DIRECTLY AFFECT HOW MUCH A CHANGE IN THE OTHER WILL ALTER THE OVERALL SIGNAL CHANGE.
IF THE GOAL IS A PROPORTIONAL PROXY DO I ACTUALLY WANT SOMETHING SIMPLER WHERE THE VARIANCE OF THE FINAL T2* TIME SERIES AND THE VARIANCE OF THE FINAL S0 TIME SERIES ARE PROPORTIONAL TO EACH OTHER?
(i.e. CREATE TIME SERIES WITH MANY ECHOES AND FIT TO GET PRECISE VALUES FOR EACH)
IF SO THEN GO BACK TO AN EARLIER STAGE AND GET THE PURE R2* and PURE S0 CHANGES TO CALCULATE THE MAX VARIANCE FOR EACH.
THEN FIND WAYS TO SHIFT BETWEEN TO SEE IF THE SCALED VARIANCE OF ONE CAN CHANGE IN PROPORTION TO THE OTHER



$var(S) = var(XY) = E[X^2*Y^2] - E[X]^2*E[Y]^2$

$E[X] = mean(X)=\overline{S_0}$ (a predefined constant)

$E[Y] = mean(Y) = e^{-TE * \overline{R_2^*}}$ (a predefined constant)

Therefore $E[X]^2*E[Y]^2$ is also a predefined constant and $var(S)$ is varied only by $E[X^2*Y^2]$

That means if $p*X^2 = (1-p)*Y^2$ then

$E[X^2*Y^2] = E[\frac{(1-p)*Y^2*Y^2}{p}] = E[\frac{X^2*X^2*p}{1-p}]$

Thus $E[\frac{(1-p)}{p}*Y^4] = E[\frac{p}{1-p}*X^4]$

$E[\frac{(1-p)}{p}*Y^4] = E[\frac{p}{1-p}*X^4]$


*************
$S = (\overline{S_0} + \Delta S_0) * e^{-TE * \overline{R_2^*}}*e^{-TE * \Delta R_2^*}$

$S = \overline{S_0}*e^{-TE * \overline{R_2^*}}*e^{-TE * \Delta R_2^*} + \Delta S_0*e^{-TE * \overline{R_2^*}}*e^{-TE * \Delta R_2^*}$ 
********
****


$var(S) = var(XY) = (E[X])^2*var(Y) + (E[Y])^2*var(X) + var(X)*var(Y)$




   
$var(X) = E[(\Delta S_0)^2]$  
$var(Y) = E[(e^{-TE * (\overline{R_2^*} + \Delta R_2^*)} - e^{-TE * (\Delta R_2^*)})^2]=E[(e^{-TE * (\Delta R_2^*)}*(e^{-TE * \overline{R_2^*}}-1))^2]$  




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