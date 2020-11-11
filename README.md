SBG.jl
=

A `julia` package for the fast simulation of the extrema of Lévy processes. Main methods:
* `rand_G` - Draws a sample from the marginal value 
    <a href="https://www.codecogs.com/eqnedit.php?latex=X_t^{(\kappa)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_t^{(\kappa)}" title="X_t^{(\kappa)}" /></a> 
  which is the Gaussian approximation of the target Lévy process 
    <a href="https://www.codecogs.com/eqnedit.php?latex=X" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X" title="X" /></a> 
  at time 
    <a href="https://www.codecogs.com/eqnedit.php?latex=t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t" title="t" /></a> 
  with cutoff level 
    <a href="https://www.codecogs.com/eqnedit.php?latex=\kappa" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\kappa" title="\kappa" /></a>.
* `rand_Gχ` - Draws an exact sample of the Gaussian approximation vector 
    <a href="https://www.codecogs.com/eqnedit.php?latex=(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" title="(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" /></a> 
  consisting of the state 
    <a href="https://www.codecogs.com/eqnedit.php?latex=X_t^{(\kappa)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_t^{(\kappa)}" title="X_t^{(\kappa)}" /></a>, 
  the supremum 
    <a href="https://www.codecogs.com/eqnedit.php?latex=\overline{X}_t^{(\kappa)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overline{X}_t^{(\kappa)}" title="\overline{X}_t^{(\kappa)}" /></a> 
  and the time the supremum is attained 
    <a href="https://www.codecogs.com/eqnedit.php?latex=\overline{\tau}_t(X^{(\kappa)})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overline{\tau}_t(X^{(\kappa)})" title="\overline{\tau}_t(X^{(\kappa)})" /></a>.
* `rand_SBG` - Draws a fast sample of the Gaussian approximation vector <a href="https://www.codecogs.com/eqnedit.php?latex=(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" title="(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" /></a>.

Table of Contents
==

1. [Input](#input) 
2. [Remarks and References](#references)
3. [Examples](#examples)
4. [Author and Contributor List](#authors)

<a name="input"/>

Input
==

* `b_κ` (or `b_κ1`, `b_κ2`) - real number - drift associated to the cutoff function 
    <a href="https://www.codecogs.com/eqnedit.php?latex=x\mapsto&space;1_{(-\kappa,\kappa)}(x)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x\mapsto&space;1_{(-\kappa,\kappa)}(x)" title="x\mapsto 1_{(-\kappa,\kappa)}(x)" /></a>.
* `σ` - nonnegative number - Brownian component.
* `sν_p` - function - magnitude batch sampler from the positive tail measure 
  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap[x,\infty))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap[x,\infty))" title="\nu(\cdot\cap[x,\infty))" /></a>, 
i.e. `sν_p(x)` gives a 
  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathrm{Poisson}(\nu([x,\infty)))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathrm{Poisson}(\nu([x,\infty)))" title="\mathrm{Poisson}(\nu([x,\infty)))" /></a>
  number of positive variables with law 
  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap[x,\infty))/\nu([x,\infty))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap[x,\infty))/\nu([x,\infty))" title="\nu(\cdot\cap[x,\infty))/\nu([x,\infty))" /></a>.
* `sν_n` - function - magnitude batch sampler from the negative tail measure <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap(-\infty,-x]))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap(-\infty,-x])" title="\nu(\cdot\cap(-\infty,-x])" /></a>, i.e. `sν_n(x)` gives a 
  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathrm{Poisson}(\nu((-\infty,x])))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathrm{Poisson}(\nu((-\infty,x]))" title="\mathrm{Poisson}(\nu((-\infty,x]))" /></a> number of positive variables with law 
  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap(-\infty,--x])/\nu((-\infty,-x])" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap(-\infty,--x])/\nu((-\infty,-x])" title="\nu(\cdot\cap(-\infty,--x])/\nu((-\infty,-x])" /></a>.
* `sum_sν_p` - function - sampler from the sums of the jumps from the positive tail measure <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap[x,\infty))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap[x,\infty))" title="\nu(\cdot\cap[x,\infty))" /></a>, i.e. `sum_sν_p(x)` has the same law as `sum(sν_p(x))`.
* `sum_sν_n` - function - sampler from the sums of the jumps from the negative tail measure <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\nu(\cdot\cap(-\infty,-x])" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\nu(\cdot\cap(-\infty,-x])" title="\nu(\cdot\cap(-\infty,-x])" /></a>, i.e. `sum_sν_n(x)` has the same law as `sum(sν_n(x))`.
* `σ_κ` (or `σ_κ1`, `σ_κ2`) - nonnegative number - standard deviation of small jumps 
    <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\sqrt{\int_{(-\kappa,\kappa)}x^2\nu(\mathrm{d}x)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\sqrt{\int_{(-\kappa,\kappa)}x^2\nu(\mathrm{d}x)}" title="\sqrt{\int_{(-\kappa,\kappa)}x^2\nu(\mathrm{d}x)}" /></a>.
* `t` - nonnegative number - time horizon.
* `κ` (or `κ1`, `κ2`) - nonnegative number - the cutoff level for small jumps.

The meaning behind the arguments `sum_sν_p` and `sum_sν_n` is that there may be faster ways of coding such functions instead of just calling `sum(sν_p(x))` and `sum(sν_n(x))`. For instance, multithreading may be easilly used.

Methods
==

```julia
rand_G
```
This method has two versions:
* `rand_G(b_κ,σ,sum_sν_p,sum_sν_n,σ_κ,t,κ)` - whose output is 
  <a href="https://www.codecogs.com/eqnedit.php?latex=X_t^{(\kappa)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_t^{(\kappa)}" title="X_t^{(\kappa)}" /></a>.
* `rand_G(b_κ1,b_κ2,σ,sum_sν_p,sum_sν_n,σ_κ1,σ_κ2,t,κ1,κ2)` - whose output is the pair 
  <a href="https://www.codecogs.com/eqnedit.php?latex=(X_t^{(\kappa_1)},X_t^{(\kappa_2)})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(X_t^{(\kappa_1)},X_t^{(\kappa_2)})" title="(X_t^{(\kappa_1)},X_t^{(\kappa_2)})" /></a> 
under a synchronous coupling.

```julia
rand_Gχ
```
This method has two versions:
* `rand_Gχ(b_κ,σ,sν_p,sν_n,σ_κ,t,κ)` - whose output is 
  <a href="https://www.codecogs.com/eqnedit.php?latex=(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" title="(X_t^{(\kappa)},\overline{X}_t^{(\kappa)},\overline{\tau}_t(X^{(\kappa)}))" /></a>.
* `rand_Gχ(b_κ1,b_κ2,σ,sν_p,sν_n,σ_κ1,σ_κ2,t,κ1,κ2)` - whose output is the pair 
  <a href="https://www.codecogs.com/eqnedit.php?latex=(X_t^{(\kappa_1)},\overline{X}_t^{(\kappa_1)},\overline{\tau}_t(X^{(\kappa_1)}),X_t^{(\kappa_2)},\overline{X}_t^{(\kappa_2)},\overline{\tau}_t(X^{(\kappa_2)}))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(X_t^{(\kappa_1)},\overline{X}_t^{(\kappa_1)},\overline{\tau}_t(X^{(\kappa_1)}),X_t^{(\kappa_2)},\overline{X}_t^{(\kappa_2)},\overline{\tau}_t(X^{(\kappa_2)}))" title="(X_t^{(\kappa_1)},\overline{X}_t^{(\kappa_1)},\overline{\tau}_t(X^{(\kappa_1)}),X_t^{(\kappa_2)},\overline{X}_t^{(\kappa_2)},\overline{\tau}_t(X^{(\kappa_2)}))" /></a> 
under a relatively poor coupling.

```julia
rand_SBG
```
This method has two versions:
* `rand_SBG(b_κ,σ,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_κ,n,t,κ)` - a faster version of `rand_Gχ(b_κ,σ,sν_p,sν_n,σ_κ,t,κ)` where it is suggested to take 
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;n=5&plus;\lfloor\nu(\mathbb{R}\setminus(-\kappa,\kappa))\rfloor" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;n=5&plus;\lfloor\nu(\mathbb{R}\setminus(-\kappa,\kappa))\rfloor" title="n=5+\lfloor\nu(\mathbb{R}\setminus(-\kappa,\kappa))\rfloor" /></a>.
* `rand_SBG(b_κ1,b_κ2,σ,sν_p,sν_n,sum_sν_p,sum_sν_n,σ_κ1,σ_κ2,n,t,κ1,κ2)` - a faster version of `rand_Gχ(b_κ1,b_κ2,σ,sν_p,sν_n,σ_κ1,σ_κ2,t,κ1,κ2)` with an improved coupling, where it is suggested to take 
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;n=5&plus;\lfloor\nu(\mathbb{R}\setminus(-\kappa_2,\kappa_2))\rfloor" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;n=5&plus;\lfloor\nu(\mathbb{R}\setminus(-\kappa_2,\kappa_2))\rfloor" title="n=5+\lfloor\nu(\mathbb{R}\setminus(-\kappa_2,\kappa_2))\rfloor" /></a>.

<a name="references"/>

## Remarks and References

The details behind the algorithms can be found in the article: 
Jorge González Cázares and Aleksandar Mijatović, *Simulation of the drawdown and its duration in Lévy models via stick-breaking Gaussian approximation*, [arXiv:1806.01870v2](https://arxiv.org/abs/1806.01870v2) (2020).

<a name="examples"/>

## Examples
The file 'Examples.jl' in the current repository includes the input functions for tempered stable and Watanabe processes.  


<a name="authors"/>

## Author and Contributor List
Jorge I. González Cázares

Aleksandar Mijatović 
