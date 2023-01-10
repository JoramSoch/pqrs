# pqrs

### The pqrs model for recognition-confidence ratings from episodic memory experiments

**Note: This is the repository containing the general toolbox code. For the data set and analysis scripts belonging to the paper describing the model, please see the repository "[FADE_pqrs](https://github.com/JoramSoch/FADE_pqrs)".**

These functions provide methods for maximum likelihood estimation or Bayesian estimation of pqrs models based on subsequent memory reports (also called "recognition-confidence ratings") acquired during episodic memory retrieval. A manuscript describing this novel computational model, entitled "A novel approach for modelling subsequent memory reports by separating decidedness, recognition and confidence", is [available from psyArXiv](https://psyarxiv.com/u5t82/) and described in [another repository](https://github.com/JoramSoch/FADE_pqrs).


## Description of the model

Given trial-wise subsequent memory reports on a Likert scale (ranging from 1 to 5) and a binary variable indicating true item status (old vs. new), the pqrs model separates behavioral responses into three distinct cognitive processes, namely *decidedness*, i.e. the tendency to give a
neutral or non-neutral response; *recognition*, i.e. the ability to label previously seen items as "old"; and *confidence*, i.e. the act of labeling an item as "definitely" rather than "probably" old or new.


## Implementation of the model

Let `y` be an `n x 1` vector of recognition-confidence ratings ranging from 1 (= this item is definitely new) over 3 (= don't know) to 5 (= this is item is definitely old) and let `x` be an `n x 1` vector coding item status as 1 (= old item) or 2 (= new item).

*Maximum likelihood estimation* of the pqrs model proceeds by calling `ME_pqrs_MLE.m`

```matlab
[p_MLE, p_lab, MLL, k] = ME_pqrs_MLE(y, x, m)
```

and *Bayesian estimation* for the pqrs model proceeds by calling `ME_pqrs_Bayes.m`

```matlab
[ab_post, p_lab, LME, k] = ME_pqrs_Bayes(y, x, m, ab_prior)
```

where `m` is a string specifying the pqrs model to be estimated (see [ME_pqrs_MLE.m](https://github.com/JoramSoch/pqrs/blob/main/ME_pqrs_MLE.m#L19-L33) or [ME_pqrs_Bayes.m](https://github.com/JoramSoch/pqrs/blob/main/ME_pqrs_Bayes.m#L21-L35)), `p_lab` is a `k x 1` cell array indicating names of model parameters resulting from `m`; `p_MLE` are [maximum likelihood estimates](https://statproofbook.github.io/D/mle) and `ab_post` are hyperparameters of the [posterior distributions](https://statproofbook.github.io/D/post); `MLL` is the [maximum log-likelihood](https://statproofbook.github.io/D/mll), `LME` is the [log model evidence](https://statproofbook.github.io/D/lme) and `k` is the number of free model parameters.

Given that model parameters have been pre-specified (or estimated), one can also generate simulated data from the model using `ME_pqrs_Sim.m`

```matlab
[y, x, m] = ME_pqrs_Sim(p, p_lab, n, o)
```

where `p` are values of ground truth parameters `p_lab`, `n` is the number of trials and `o` is the fraction of old items.
