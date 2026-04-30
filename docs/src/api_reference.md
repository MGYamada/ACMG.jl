# API Reference

This reference emphasizes stable public APIs.  Experimental functions are
listed in [API Stability](@ref) and documented in their docstrings.

## Cyclotomic and conductor utilities

```@docs
CyclotomicContext
field
zeta
conductor
cond_S
cond_T
cond_F
galois_action
galois_orbit
frobenius
reduce_mod_p
```

## Modular data and Verlinde checks

```@docs
ModularData
modular_data
semion_modular_data
fibonacci_modular_data
ising_modular_data
toric_code_modular_data
validate_modular_data
validate_exact_modular_data
check_modular_relations
check_unitarity
check_verlinde_integrality
verlinde_coefficient
extract_and_lift
```

## Higher central charges

```@docs
HigherCentralChargeResult
gauss_sum_plus
gauss_sum_minus
normalized_gauss_sum
higher_central_charge
higher_central_charges
central_charge
```

## F/R data and equations

```@docs
FusionRule
FRData
FREquationSystem
FiniteFieldEquationSystem
fr_equation_system
pentagon_equations
hexagon_equations
validate_fr_system
F_symbol
R_symbol
semion_fr_data
fibonacci_fr_data
ising_fr_data
```

## Gauge and braid representations

```@docs
GaugeTransform
GaugeParameters
GaugeChoice
GaugeFixingResult
gauge_parameters
gauge_fixing_plan
is_gauge_fixed
gauge_fix
FusionPath
FusionTreeBasis
BraidRepresentation
braid_representation
braid_generator
braid_generators
check_braid_relations
```

## Pipeline

```@docs
ClassifiedMTC
classify_mtcs_at_conductor
classify_mtcs_auto
fr_status
select_admissible_primes
recommend_primes
recommend_skip_FR
save_classification
load_classification
write_report
```
