# Gauge Fixing

The gauge layer records gauge parameters, gauge choices, and conservative
normalization metadata for F/R data.

Representative functions:

- `GaugeTransform`, `GaugeParameters`, `GaugeChoice`, `GaugeFixingResult`
- `gauge_parameters`, `gauge_transform`
- `gauge_fixing_plan`, `is_gauge_fixed`
- `gauge_variables`, `gauge_fix`

Notes: high-level records and accessors are intended for normal users.
Low-level normal-form helpers and finite-field gauge internals are
experimental/internal implementation details.
