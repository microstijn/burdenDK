## Mixture Aggregation Update

Exact analytical DEB-axis mixture aggregation was added to the runtime model, maintaining strict conservatism:
- **Exact same-axis summation only**: `s_{a,t} = sum_j s_{j,a,t}`
- **No Independent Action (IA)**
- **No Concentration Addition (CA) over raw external concentrations**
- **No bounded nonlinear aggregation**
- **No synergism or antagonism**
- **No low-effect approximation flags**
- **No physiological condition memory Z_t**
- **No raster integration**

The methods implemented (`additive_axis_burden` and `axis_toxic_unit_sum`) precisely compute the manual aggregation while providing deterministic contribution diagnostics without introducing nonlinear thresholds.
