# Compound Memory Initialisation

## Mathematical Recurrence

The chemical memory or internal burden ($B_t$) follows a discrete recurrence relation tracking the retained memory of past exposures and newly absorbed ambient concentration ($C_t$):

$$ B_t = \rho B_{t-1} + (1-\rho) K C_t $$

where:
- $B_t$ = retained internal compound burden
- $\rho$ = monthly retention / carryover factor ($0 \le \rho < 1$)
- $K$ = bioaccumulation or internal magnification factor
- $C_t$ = ambient concentration at month $t$

## Constant-Background Spin-up

To avoid initialising models from an arbitrary zero-burden state (which assumes no prior exposure history), we provide an analytical solution for the internal burden after $n$ months of spin-up under a constant background concentration $C_{bg}$:

$$ B_n = \rho^n B_0 + K C_{bg} (1-\rho^n) $$

Where $B_0$ is the initial burden (usually 0) before the spin-up period.

## Inverse Target-Burden Background

If a target internal burden ($B_{target}$) needs to be achieved after $n$ spin-up months, we can invert the above equation to determine the necessary constant background concentration $C_{bg}$:

$$ C_{bg} = \frac{B_{target} - \rho^n B_0}{K(1-\rho^n)} $$

This is particularly useful when starting from a fractional assumption (e.g. 10% of NOEC).

## Periodic-Cycle Warm-up

For a repeating ambient concentration cycle $C_1, C_2, \dots, C_n$, we want to find the initial burden $B_0$ such that, after running one full cycle, the system returns exactly to $B_0$. The closed-form analytical solution is:

$$ B_0 = \frac{(1-\rho) K \sum_{i=1}^{n} \rho^{n-i} C_i}{1-\rho^n} $$

**Convention:** `analytical_periodic_initial_burden` returns the burden *immediately before* applying the first month of the cycle ($C_1$).

## Conceptual Distinctions

It is critical to distinguish between different state variables in the modeling framework:
- **Chemical memory ($B_t$)** represents the physical accumulation or retention of a specific toxicant (internal burden).
- **$B_t$ is NOT physiological condition memory ($Z_t$)**. Condition memory tracks the biological state (e.g., organism health or resource buffers), which may recover more slowly than the chemical clears out of the organism.
- **$B_t$ is NOT DEBtox scaled damage ($D_t$)**. Scaled damage in DEBtox is a continuous dynamic variable with its own clearance rate parameter ($k_r$), and relates directly to hazard rates.

These distinct variables are why compound memory requires independent warm-up utilities.
