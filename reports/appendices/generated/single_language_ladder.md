The following ladder definitions are sourced from R/04_model_formulas.R 
via build_model_ladder(has_pseudo_passive = TRUE).

### `L5_correlated_maximal`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled | 
    Participant) + (1 + S_Type | Verb)
```

### `L4_uncorrelated_maximal`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled || 
    Participant) + (1 + S_Type || Verb)
```

### `L3_no_participant_interaction_slope`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type + Semantics_scaled || 
    Participant) + (1 + S_Type || Verb)
```

### `L2_sentence_type_slopes_only`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type || Participant) + 
    (1 + S_Type || Verb)
```

### `L1_random_intercepts_plus_participant_semantics`

```r
Response ~ S_Type * Semantics_scaled + (1 + Semantics_scaled || 
    Participant) + (1 | Verb)
```

### `L0_random_intercepts_only`

```r
Response ~ S_Type * Semantics_scaled + (1 | Participant) + (1 | 
    Verb)
```

