# MCTP simulation on fake methylation data
_Simulation of the influence of asymptotic approximations used in MCTP on fake methylation data._  

## Introduction
The rank-based multiple contrast test procedure (MCTP) proposed in [Kimihiro Noguchi et al. (2019)](https://link.springer.com/article/10.3758/s13428-019-01247-9) has been shown to be an advisable alternative to the 2 steps Kruskal-Wallis (KW) & Wilcoxon procedure to address the nontransitive paradoxes when doing multiple groups statistical comparison. The MCTP, being non-parametric, makes it suitable to non-normally distributed data. Thus, the MCTP should be applicable to methylation data, more precisely beta values, to perform multiple samples (or groups) comparisons.  
In this repository, I propose a simulation of the MCTP applied to randomly generated beta values. The simulation has been carried out onto two distinct custom distributions: a classical bimodal beta distribution, and a trimodal distribution more representative of some bulk methylation data distribution. All asymptotic approximations supported by the MCTP were used to assess their impact on the type-I error rate of the procedure following the two custom distributions.  

## Random generation of fake beta values for the simulation
In order to simulate the effect of the different asymptotic approximations on the type-I error of the MCTP, I first needed to generate truthful distributions, representative of the usual distribution observed for methylation beta values.  

A notorious distribution, representative of the usual methylation beta values distribution, is the bimodal beta distribution. The beta distribution is a continuous probability distribution defined on the interval (0-1), and governed by two parameters, α and β, controlling the shape of the distribution. Maxima of a beta distribution are usually expected around x = 0 and x = 1.  
In order to create such distribution I made use of the `rbeta()` function in R as follows:  
```R
rbimodal_meth <- function(n){
  bimodal_distrib <- rbeta(n = n, shape1 = 0.2, shape2 = 0.2, ncp = 0)
  return(bimodal_distrib)
}
```
Sometimes, it happens in bulk methylation data that, instead of following a bimodal distribution, methylation beta values follow a trimodal distribution. This distribution is also defined on the interval (0-1), but it possesses 3 local maxima, usually found around x = 0, x = 0.5 and x = 1, respectively.  
There is no straight forward way to generate such trimodal distribution in R. In order to create it, I made use of the `rtruncnorm()` function from the [truncnorm R package (Olaf Mersmann et al. 2023)](https://CRAN.R-project.org/package=truncnorm) as follows:  
```R
rtrimodal_meth <- function(n){
  low_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 0, sd = 0.05)
  high_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 1, sd = 0.05)
  mid_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.05)
  trimodal_distrib <- sample(
    x = c(low_distrib, mid_distrib, high_distrib), size = n)
  return(trimodal_distrib)
}
```

## Simulation of asymptotic approximations effects on the MCTP type-I error rate

In the [nparcomp R package (Frank Konietschke et al. 2015)](https://www.jstatsoft.org/article/view/v064i09) I used the asymptotic approximations supported by the `mctp()` function: **"fisher", "mult.t", "log.odds and "normal"** to test what effect each one of these had over the MCTP type-I error rate.  
`mctp()` has been used to perform a **"many-to-many"** comparison between the 3 groups, using the **"Tukey" contrast**.  
For each asymptotic approximation supported, and each type of fake methylation data distributions (`rbimodal_meth()` and `rtrimodal_meth()`) I have ran `mctp()` **1,000 times**, creating the 3 groups using different sample sizes (**N=10, N=100 and N=1000**).   

I collected all p.values returned by the `mctp()` function for each group-to-group comparisons, and calculated the percentage of those falling below 0.05 (**p.value <= 0.05**) at random.  
Results of the simulation are shown in the table below:  

|Distributions      |Asympt. approx. |Sample size |Sim. err. rate |
|:------------------|:---------------|:-----------|:--------------|
|`rbimodal_meth()`  |fisher          |10          |1.53%          |
|`rbimodal_meth()`  |fisher          |100         |1.53%          |
|`rbimodal_meth()`  |fisher          |1000        |2.03%          |
|`rbimodal_meth()`  |mult.t          |10          |2.83%          |
|`rbimodal_meth()`  |mult.t          |100         |1.63%          |
|`rbimodal_meth()`  |mult.t          |1000        |1.47%          |
|`rbimodal_meth()`  |log.odds        |10          |1.83%          |
|`rbimodal_meth()`  |log.odds        |100         |2%             |
|`rbimodal_meth()`  |log.odds        |1000        |2.13%          |
|`rbimodal_meth()`  |normal          |10          |3.97%          |
|`rbimodal_meth()`  |normal          |100         |1.8%           |
|`rbimodal_meth()`  |normal          |1000        |1.27%          |
|`rtrimodal_meth()` |fisher          |10          |0.43%          |
|`rtrimodal_meth()` |fisher          |100         |0.37%          |
|`rtrimodal_meth()` |fisher          |1000        |0.57%          |
|`rtrimodal_meth()` |mult.t          |10          |0.73%          |
|`rtrimodal_meth()` |mult.t          |100         |0.23%          |
|`rtrimodal_meth()` |mult.t          |1000        |0.53%          |
|`rtrimodal_meth()` |log.odds        |10          |0.53%          |
|`rtrimodal_meth()` |log.odds        |100         |0.3%           |
|`rtrimodal_meth()` |log.odds        |1000        |0.47%          |
|`rtrimodal_meth()` |normal          |10          |1.57%          |
|`rtrimodal_meth()` |normal          |100         |0.57%          |
|`rtrimodal_meth()` |normal          |1000        |0.8%           |

Results from the simulation suggest that the MCTP using a normal approximation was associated with the lowest type-I error rate on bimodal methylation distributions. However, on trimodal methylation distributions, the MCTP using a log.odds approximation was associated with the lowest type-I error.  

## Discussion
Generating fake credible methylation beta values is not trivial, and the approaches I chosed here are most probably perfectible. These do not take into account common batch effects that occur during the beta values aquisition which often skew distributions. This random generation approach also doesn't necessarily represent very well potential artifacts during methylation measurements. it also doesn't represent very well the potential inter-sample variabilities existing within a given cohort. Moreover, specific conditions like cancer or auto-immune diseases are known to affect the overall methylation data distribution in a sample or an entire cohort.  
In methylation data analysis the bimodal distribution is more often encountered than the trimodal distribution, which only may occur occasionnaly, in specific dataset, in specific conditions, when using low-resolution methods (e.g. methylation array data) on biopsies containing heterogeous cell types.  
For these reasons, I would prefer the results obtained on bimodal random distributions over those obtained on trimodal distributions.  
The parameters (`shape1`, `shape2`, and `ncp`) used to generate the random bimodal distributions could have also been teaked, to consolidate the simulation results and assess the influence these parameters had over the MCTP type-I error rate.  

## Code access and reproducibility
All scripts used in the simulation are available in this repository. As the simulation is based on randomly generated data, it thus won't be exactly reproductible, and the results won't exactly be the same as reported here.  

## Acknowledgement
I would like to thank Dr. Kerstin Rubarth for her expertise and her support over the use of the MCTP, and her suggestions on carrying out and sharing the results of this simulation.  

## References
1. [Noguchi, K., Abel, R.S., Marmolejo-Ramos, F. et al. Nonparametric multiple comparisons. Behav Res 52, 489–502 (2020). https://doi.org/10.3758/s13428-019-01247-9.](https://link.springer.com/article/10.3758/s13428-019-01247-9)  
2. [Olaf Mersmann, Heike Trautmann, Detlef Steuer, Björn Bornkamp. truncnorm: Truncated Normal Distribution (2023).](https://CRAN.R-project.org/package=truncnorm)  
3. [Konietschke F, Placzek M, Schaarschmidt F, Hothorn LA (2015). “nparcomp: An R Software Package for Nonparametric Multiple Comparisons and Simultaneous Confidence Intervals.” Journal of Statistical Software, 64(9), 1–17. http://www.jstatsoft.org/v64/i09/.](https://www.jstatsoft.org/article/view/v064i09)  

