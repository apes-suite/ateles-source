title: Postprocess Filtering

High-order discretizations suffer from Gibbs oscillations in the proximity of
discontinuities.
The idea with postprocess filtering is to leave those oscillations in the
simulation and only remove them in a postprocessing step. This allows the
computation to maintain strong gradients as long as we do not run into
unphysical values due to the oscillations.
Theoretically, if the position of the discontinuity is known, it is also
possible to recover a high-order convergence despite the discontinuity.
However, positions are generally not known and not too trivial to find.

In `atl_harvesting` we offer two solutions to removing oscillations.
The first is to simply apply a strong spectral filtering (possibly with
covolume filtering) on the loaded data before refinement.
These are controlled by the usual `scheme.stabilization` table.
However, you need to explicitly activate them with the additional
`use_post_filter = true` setting.
This is the simple most option, but it also quite clearly smoothes out
discontinuities.

An alternative is invited by the adaptive refinement, as it should be
able to narrow down discontinuities.
We implemented here the `filter_element` option that allows us to
filter high modes in each refinement step based on the polynomial data.
The filtering is controlled by the `filter_element` table within the
`ply_sampling` part:

```
  filter_element = {
    strategy = 'oddfract',
    min_order = 2,
    max_order = 10,
    fract_exponent = 3
  }
```

Currently we only have the 'oddfract' strategy, which uses the energy contained
in the odd modes compared to the overall spectral energy.
Only the odd modes contribute to the average derivative and would be stronger
if there is a discontinuity there.
The filter order is based on this fraction, the larger the fraction of odd
modes, the larger the filter order of the used spectral filter.
We then apply a spectral filtering in each refinement.
The filter order for this filtering on each level is then ranging between
`min_order` (small fraction of odd modes) and `max_order` (high fraction of
odd modes).
The formula for the spectral filter order is given by:
`min + (max-min+1)*odd_fract^fract_exponent`.
It should be noted, that this filter is applied in each direction and thus,
there is a multiplicative effect for multiple dimensions.
Also keep in mind that there is a multiplicative effect from refinement to
refinement, as the filter is applied in each step.
Note, that the refinement will reduce the last mode by `2^(degree+1)` and
the filtering will anyway reduce it to 0. Thus, removing the last mode will
not drastically affect the refined solution.
To just remove the last mode during refinement from further consideration,
use the decrement strategy for dof reduction.
Note that this will be quite costly an can result in an greatly increased
memory consumption. But cutting too many modes during refinement will result
in less smooth adaptation to the shocks.
