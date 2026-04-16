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

$S(v,t,TE) = (\overline{S_0(v,TE)} + \Delta S_0(v,t,TE)) * e^{-TE * (\overline{R_2^*(v,TE)} + \Delta R_2^*(v,t,TE))}$

For everything below, we will get different values for each voxel (v),
but removing the v to make the equation slightly cleaner

$S(t,TE) = (\overline{S_0(TE)} + \Delta S_0(t,TE)) * e^{-TE * (\overline{R_2^*(TE)} + \Delta R_2^*(t,TE))}$

The overall goal is for user specified values for $S(t,TE)$, $\overline{S_0(TE)}$, $\overline{R_2^*(TE)}$,
for a reference $TE$, create time series that are
a roughly linear proportion of contributions from $\Delta S_0(t,TE)$ and $\Delta R_2^*(t,TE)$.
That is, for the reference TE, regardless of the proportion, $S(t,TE)$ is exactly the prespecified value,
but at other echo times, a new $S(t,TE)$ can be calculated using the same
$\overline{S_0(TE)}$, $\overline{R_2^*(TE)}$, $\Delta S_0(t,TE)$, and $\Delta R_2^*(t,TE)$

The reason this is useful is that a user could then specify exactly what their data look
like at a specific echo time, even by providing real data,
and then model or simulate what the time series would be at other echo times
depending on the relative contributions of $\Delta S_0(t,TE)$ and $\Delta R_2^*(t,TE)$.
If done correctly, we should be able to estimate the $\Delta S_0(t,TE)$ and $\Delta R_2^*(t,TE)$ from simulated data (using kappa and rho or other metrics) while also knowing the
true values.

$S(t,TE) = (\overline{S_0(TE)} + \Delta S_0(t,TE)) * e^{-TE * (\overline{R_2^*(TE)} + \Delta R_2^*(t,TE))}$

$S(t,TE) = \overline{S_0(TE)}*(1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * e^{-TE * \overline{R_2^*(TE)}}*e^{-TE * \Delta R_2^*(t,TE)}$

$\overline{S(TE)} = \overline{S_0(TE)} * e^{-TE * \overline{R_2^*(TE)}}$

$S(t,TE) = \overline{S(TE)}*(1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * e^{-TE * \Delta R_2^*(t,TE)}$

Percent change from the mean

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} = (1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * e^{-TE * \Delta R_2^*(t,TE)} - 1$

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} = (1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * e^{-TE * \Delta R_2^*(t,TE)} -1$

First order Taylor Explansion:
$e^{-TE * \Delta R_2^*(t,TE)} \approx 1 - TE * \Delta R_2^*(t,TE)$

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}}\approx (1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * (1 - TE * \Delta R_2^*(t,TE)) - 1$

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} \approx \frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}-TE * \Delta R_2^*(t,TE) - \frac{TE * \Delta R_2^*(t,TE) * \Delta S_0(t,TE)}{\overline{S_0(TE)}}$

Common assumption. This is small enough that it can be dropped: $\frac{TE * \Delta R_2^*(t,TE) * \Delta S_0(t,TE)}{\overline{S_0(TE)}}$

This is not a great assumption.
If mean $\Delta R_2^*(t,TE)$ is 30 and the $\Delta {T_2^*}=-2$ then
$TE * \Delta R_2^*(t,TE) = TE * (1/30-1/28) = - TE * 0.002$.
For a third echo echo time of $TE=50ms$ that is -0.12.
That is a scaling factor of 1.12 is added to the the reduced $\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}$ term.

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} \approx 1.12*\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}-TE * \Delta R_2^*(t,TE)$

If $\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}=0.1$ (a 10% signal change) then, with the above assumption:
$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} \approx 0.1 + 0.12 = 0.22$
Without the above assumption,
$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} \approx 1.12*0.1 + 0.12 = 0.232$
That's about a 5% misapproximation of $\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}}$

Going forward this with approximation/assumption for now and will test overall error in the simulation results later.

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} \approx \frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}-TE * \Delta R_2^*(t,TE)$

$\frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} = S_{spc}(t,TE) = signal\ percent\ change\ from\ mean$

For scaling the proportion of each:
Let p be the proportion of signal driven by $\Delta S_0(t,TE)$ changes

$p*S_{spc}(t,TE) \approx \frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}$

$(1-p)*S_{spc}(t,TE) \approx -TE * \Delta R_2^*(t,TE)$

That is the $\Delta S_0(t,TE)$ and the $\Delta R_2^*(t,TE)$ parts of the equation should each model a specified fraction of the total $S_{spc}(t,TE)$

$\Delta S_0(t,TE) \approx p*S_{spc}(t,TE)*\overline{S_0(TE)}$

$\Delta R_2^*(t,TE) \approx -\frac{(1-p)*S_{spc}(t,TE)}{TE}$

Note: The approximation above means that the approximation error will be smallest for p=0 or 1 and largest somewhere where there is a mix

For any given $S(t,TE)$, $\overline{S_0(TE)}$, $\overline{R_2^*(TE)}$, $TE$ & $p$ we can thus approximate $\Delta S_0(t,TE)$ and $\Delta R_2^*(t,TE)$.
There are several levels of tests to see how far this approximation diverges from truth.

Full decay curve: 
$S_{spc} = \frac{S(t,TE)-\overline{S(TE)}}{\overline{S(TE)}} = (1+\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}) * e^{-TE * \Delta R_2^*(t,TE)} -1$

Without any approximations, the above should be identical to the prespecified $S(t,TE)$ for the pre-specified $TE$.
For other $TEs$ the above is our closest value to truth can can be compared to the other levels of approximation.

First order Taylor approximation:
$S_{spc} \approx \frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}-TE * \Delta R_2^*(t,TE) - \frac{TE * \Delta R_2^*(t,TE) * \Delta S_0(t,TE)}{\overline{S_0(TE)}}$

Taylor approximation with dropping $\frac{TE * \Delta R_2^*(t,TE) * \Delta S_0(t,TE)}{\overline{S_0(TE)}}$

$S_{spc} \approx \frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}}-TE * \Delta R_2^*(t,TE)$

By definition, the above should show a linear relationship with p at all echo times,
but the previous two will be non-linear, & a Q is how non-linear.

Alternative: If we want proportional variance, then we'd need to solve the above for variance.

$p*var(S_{spc}(t,TE)) \approx var(\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}})$

$p*S_{spc}(t,TE)^2 \approx (\frac{\Delta S_0(t,TE)}{\overline{S_0(TE)}})^2$

$\Delta S_0(t,TE)^2\approx p*S_{spc}(t,TE)^2 \overline{S_0(TE)}^2$

$\bold{\Delta S_0(t,TE) \approx \sqrt{p}*S_{spc}(t,TE) \overline{S_0(TE)}}$

$(1-p)*var(S_{spc}(t,TE)) \approx var(-TE * \Delta R_2^*(t,TE))$

$(1-p)*S_{spc}(t,TE)^2 \approx TE^2 * \Delta R_2^*(t,TE)^2$

$\Delta R_2^*(t,TE) \approx \sqrt{\frac{(1-p)*S_{spc}(t,TE)^2}{TE^2}}$

$\bold{\Delta R_2^*(t,TE) \approx \frac{\sqrt{1-p}*S_{spc}(t,TE)}{TE}}$