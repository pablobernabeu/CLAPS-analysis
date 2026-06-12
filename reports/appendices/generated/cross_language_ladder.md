The following ladder definitions are sourced from R/04_model_formulas.R 
via build_multilanguage_ladder().

### `L5_cross_maximal`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled | 
    Participant) + (1 + S_Type * Semantics_scaled | Verb) + (1 + 
    S_Type * Semantics_scaled | Language)
```

### `L4_cross_uncorrelated`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type * Semantics_scaled || 
    Participant) + (1 + S_Type * Semantics_scaled || Verb) + 
    (1 + S_Type * Semantics_scaled || Language)
```

### `L3_cross_no_participant_interaction`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type + Semantics_scaled || 
    Participant) + (1 + S_Type * Semantics_scaled || Verb) + 
    (1 + S_Type * Semantics_scaled || Language)
```

### `L2_cross_stype_participant_verb`

```r
Response ~ S_Type * Semantics_scaled + (1 + S_Type || Participant) + 
    (1 + S_Type || Verb) + (1 + S_Type * Semantics_scaled || 
    Language)
```

### `L1_cross_intercepts_only_ppt_verb`

```r
Response ~ S_Type * Semantics_scaled + (1 | Participant) + (1 | 
    Verb) + (1 + S_Type || Language)
```

### `L0_cross_intercepts_only`

```r
Response ~ S_Type * Semantics_scaled + (1 | Participant) + (1 | 
    Verb) + (1 | Language)
```

