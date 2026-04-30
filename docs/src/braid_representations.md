# Braid Representations

The braid-representation layer constructs braid group matrices from F/R data.

Representative functions:

- `FusionPath`, `FusionTreeBasis`
- `braid_representation`
- `braid_generator`, `braid_generators`
- `check_braid_relations`
- `reduce_mod_p` for braid representations

Notes: exact braid representation construction and braid-relation checking are
the stable surface.  Finite-field reductions and image diagnostics are
experimental and should be version-pinned when used in research scripts.
