The following empirical anchors are sourced from R/03_define_priors.R:

| Anchor | Value |
|---|---:|
| `semantics_pooled` | 0.47 |
| `s_type_active_interaction` | -0.31 |
| `s_type_pseudo_passive_interaction` | -0.36 |
| `semantics_min` | 0.27 |
| `semantics_max` | 0.8 |

Prior regimes are sourced from R/03_define_priors.R, where N denotes a normal prior
and t a Student-t prior:

| Regime | default | semantics | active | pseudo | Intercept | sd | cor |
|---|---|---|---|---|---|---|---|
| primary | N(0, 1.5) | N(0, 0.5) | N(0, 0.5) | N(0, 0.6) | t(3, 0, 2.5) | t(3, 0, 1) | lkj(2) |
| weak | N(0, 2) | N(0, 1) | N(0, 1) | N(0, 1) | t(3, 0, 2.5) | t(3, 0, 2) | lkj(1) |
| literature_centred | N(0, 1.5) | N(0.47, 0.35) | N(-0.31, 0.4) | N(-0.36, 0.5) | t(3, 0, 2.5) | t(3, 0, 1) | lkj(2) |
| heavy_tailed | t(3, 0, 1.5) | t(3, 0, 0.5) | t(3, 0, 0.5) | t(3, 0, 0.6) | t(3, 0, 2.5) | t(3, 0, 1) | lkj(2) |
