# Oracle Consistency Analysis Results

## Summary Metrics

```
    config mean_optimal sd_optimal min_optimal max_optimal range_optimal
  baseline           49      29.23          20          90            70
    min_50           42      26.58          20          90            70
   min_100           42      20.98          20          90            70
  equal_75           41      25.14          20         100            80
  equal_50           29      11.97          20          50            30
 equal_100           42      20.98          20          90            70
 mad_optimal pct_at_10 pct_at_100 pct_at_boundaries n_comparisons
       29.65         0          0                 0            10
       14.83         0          0                 0            10
        7.41         0          0                 0            10
        7.41         0         10                10            10
        7.41         0          0                 0            10
        7.41         0          0                 0            10
```

## Recommendation

**Best configuration: equal_50**

- Optimal nbasis: 29.0 ± 12.0 (mean ± SD)
- Range: 20 - 50
- Boundary selections: 0.0%

## Plots

![Oracle Curves](oracle_curves_faceted.png)

![Optimal nbasis Distribution](optimal_nbasis_boxplot.png)

