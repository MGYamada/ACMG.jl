# Block-U Search

The search layer enumerates conductor strata and performs finite-field block-U
searches used by the conductor-first pipeline.

Representative functions:

- `enumerate_strata`, `count_strata`, `describe_stratum`
- `find_mtcs_at_prime`
- `classify_mtcs_at_conductor`, `classify_mtcs_auto`
- `estimate_search_complexity`, `recommend_primes`

Notes: the public search entry points are the conductor-level pipeline
functions.  Low-level enumeration helpers in `Search/BlockU.jl` are classified
as experimental/internal in v0.8.5 even if exported for compatibility.
