title: Neumann boundary conditions

In our implementation of boundary conditions we prescribe the "outer"
state and use that to compute a flux on the boundary.
For Neumann boundary conditions it is not so obvious, which value to
use for this outer state.
The default approach is to simply use the solution in the domain
at the boundary also for the "outer" state.
However, this does not really imply any gradient as actually required
for the Neumann boundaries.
To get a better approximation for the boundary condition, we need to
take the gradient normal to the boundary into account.
This can be activated for individual boundary conditions via the
`enforce_zero_grad` option.

For example like this:
~~~~~~~~.lua
boundary_condition = {
    {
         label = 'west',
         kind = 'outflow',
         enforce_zero_grad = true,
         pressure = p
    }
}
~~~~~~~~

When a boundary is set to enforce a zero gradient, we compute the
gradient at the side and instead of using the original value at the
boundary we use the one that we get when changing the last mode to obtain
a zero gradient.

For increased stability it might be desirable to cut off even more
modes for the extrapolation.
This can be controlled by the `neumann_mode_fraction` \(\phi\) setting.
The fraction defines how many modes \(N\) of the original polynomial
 \(p\) in boundary normal direction should be used for the extrapolation:
\[ N = \phi \cdot p \]
A setting of `neumann_mode_fraction = 0` results in just the integral
mean value being used for the boundary value.
The default of `neumann_mode_fraction = 1` results in the value from
changing just the last Legendre mode of the original polynomial to meet
the zero gradient constraint at the boundary.

To find the suitable extrapolation value for the boundary, taking the
gradient into account, we exploit some properties of the Legendre
polynomials.
All Legendre polynomials are 1 at the right end of the interval:
\[ L_i(1) = 1 \]
They are symmetric (even modes) or antisymmetric (odd modes):
\[ L_i(-x) = (-1)^i L_i(x) \]
Thus, we get at the left end of the interval:
\[ L_i(1) = (-1)^i \]

With the Legendre series of the form
\[ P_N(x) = \sum_{i=0}^N c_i L_i(x) \]
we therefore, obtain at the element boundaries:
\[ P_N(-1) = \sum_{i=0}^N (-1)^i c_i \]
and
\[ P_N(1) = \sum_{i=0}^N c_i \]

The derivative of a Legendre polynomial can be expressed recursively
by:
\[ L_i'(x) = (2i-1) L_{i-1}(x) + L_{i-2}'(x) \]
and the first two are:
\[ L_0'(x) = 0 \]
\[ L_1'(x) = L_0(x) = 1 \]
As can be seen from this recursion, the derivative only consists of
even Legendre modes for odd Legendre modes to derive and vice-versa.
The derivative at the right end of the interval is given by:
\[ L_i'(1) = \frac{i \cdot (i+1)}{2} \]
And due to the symmetry property for the derivatives, on the left end
the value is accordingly given by:
\[ L_i'(-1) = (-1)^{i+1} \frac{i \cdot (i+1)}{2} \]

For a polynomial series \( P_N(x) \) we now can also compute the derivative
at the two end points:
\[ P_N'(1) = \sum_{i=0}^{N} c_i L_i'(1) \]
\[ P_N'(1) = \sum_{i=0}^{N} c_i \frac{i \cdot (i+1)}{2} \]
and
\[ P_N'(-1) = \sum_{i=0}^{N} c_i L_i'(-1) \]
\[ P_N'(-1) = \sum_{i=0}^{N} c_i (-1)^{i+1} \frac{i \cdot (i+1)}{2} \]

For Neumann boundaries we assume a zero gradient at the boundary, and
we prescribe the value that we get when we use a polynomial, for which
this condition \( P_N'(1) = 0 \) or \( P_N'(-1) = 0 \) is met.
We can find this polynomial with minimal modifications by changing the
last mode of the polynomial.
Thus, we compute a new polynomial \( \rho_N(x) \), which has all the
lower modes from the original polynomial but the last mode taken such,
that the gradient at the boundary is \(0\).
For this we introduce a correcting term:
\[ \chi(1) = \rho_{N-1}'(1) = P_{N-1}'(1)\]
\[ \chi(-1) = \rho_{N-1}'(-1) = P_{N-1}'(-1)\]

For the newly computed polynomial \( \rho_N(x) \), we want to ensure
\[ \rho_N'(1) = \chi(1) + c_R \frac{N \cdot (N+1)}{2} = 0 \]
or
\[ \rho_N'(-1) = \chi(-1) + c_L (-1)^{N+1} \frac{N \cdot (N+1)}{2} = 0 \]
Here, \(c_L\) is the last mode that is computed according to this
constraint if the zero gradient should be at the left side and \(c_R\) if
should be on the right side. From this we get:
\[ c_R = - \chi(1) \frac{2}{N \cdot (N+1)} \]
and
\[ c_L = (-1)^N \chi(-1) \frac{2}{N \cdot (N+1)} \]

Finally the value of this new, corrected polynomial at the boundary
can be computed by plugging this into the series evaluation:
\[ \rho_N(1) = \rho_{N-1}(1) + c_R \]
\[ \rho_N(1) = \rho_{N-1}(1) - \chi(1) \frac{2}{N \cdot (N+1)} \]
and on the left side:
\[ \rho_N(-1) = \rho_{N-1}(-1) + (-1)^N c_L \]
\[ \rho_N(-1) = \rho_{N-1}(-1) + (-1)^N (-1)^N \chi(-1) \frac{2}{N \cdot (N+1)}\]
\[ \rho_N(-1) = \rho_{N-1}(-1) + \chi(-1) \frac{2}{N \cdot (N+1)}\]

To use the same code for left and right, we introduce a factor to take the
side (left or right) into account. This sidefactor \(s_f\) is \(1\) for
right and \(-1\) for left.
With this we get:

\[ \chi = \sum_{i=1}^{N-1} s_f^{i+1} \cdot c_i \frac{(i+1) \cdot i}{2} \]
and
\[ \rho_{N-1} = c_0 + \sum_{i=1}^{N-1} s_f^i \cdot c_i \]
resulting in the end in:
\[ \rho_{N} = \rho_{N-1} - \frac{2 s_f}{(N+1) \cdot N} \chi \]

Note that the factor of two within the sum cancels out with the
final multiplication and can be left out completely.
