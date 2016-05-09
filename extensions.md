# Extending glmmADMB

I started to try to extend glmmADMB to allow non-trivial models for
the zero-inflation and dispersion components of the models (e.g.
allow them to differ between groups via ~group).

I came up with what I think is a reasonably nice extension of the
interface: instead of the 'formula' argument being a single formula
(e.g. `counts~lat+long`), it can be a list with components referring
to the different parts of the model, e.g. `list(eta=counts~lat+long,
zi~species,alpha~species)` (I would like a better name than `eta`
for the formula specifying the linear model for the mean of the
linear predictor ... `fixed` isn't quite right, because all of the
model components we're talking about here are fixed effects (and
in fact the `random` argument could be extended to allow lists in
the same way).

However, I get into a little more trouble trying to extend the
back end while keeping it both backward-compatible and clean.
Right now the TPL file has parameters `pz` (zero-inflation probability)
and `log_alpha` (log of the overdispersion parameter), which are
appropriately bounded (`pz` bounds = [1e-6,0.999]; `log_alpha` bounds = [-5,6]
unless fitting NB1, in which case [0.001,6]), and whose phase
is set negative when fitting models that don't involve these parameters.

My initial thought was just to turn these parameters into vectors; in
the standard backward-compatible case (single parameters), these would
just be "intercept-only" models (we could take shortcuts and not bother
to multiply the (constant) design matrix by the (single) parameter in
the trivial case, but conceptually it would be the same).  *However*, there's
a little bit of a hitch with the way that `pz` and `log_alpha` are currently
implemented/bounded.

`pz` is estimated on the raw [0,1] scale, and bounded at the edges of the range.  This
makes a little less sense when `pz` is based on a non-trivial model. Either (1) we
should switch to a transformed scale (the most obvious choice is to fit logit(pz)),
or (2) we need to bound the result of `pz_X*pz` (i.e. the design matrix times the
parameter vector) between [0,1], presumably using `fpos` to penalize results
outside of the allowable range.  Choice #1 is simpler, but is not backward
compatible and may lead to convergence problems if the best fit is at `pz=0` (although
bounding logit(`pz`) at -14 should be about the same as bounding `pz` at 1e-6).

A similar although slightly easier situation holds for `log_alpha`.  The dispersion
parameter is already being estimated on a log scale, which makes life a little easier.
The only real issue is what to do about the NB1 case, where the lower bound on the
parameter is 1, not 0: should we fit on a transformed scale `alpha=1+exp(alpha)`?

Suggestions?
