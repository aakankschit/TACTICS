# Search Performance Metrics for TACTICS

## Tracking, Visualization, and Improvement of Chemical Space Exploration

---

## Abstract

This document defines a comprehensive set of metrics for tracking and evaluating the performance of TACTICS selection strategies during chemical space exploration. We categorize metrics into five dimensions: (1) exploration vs exploitation balance, (2) search efficiency, (3) component-level dynamics, (4) posterior quality, and (5) diversity measures. Each metric is formally defined with mathematical specifications, implementation guidance, and visualization recommendations. These metrics provide actionable insights for developers to improve search strategies in future releases and enable users to understand algorithm behavior during screening campaigns.

---

## Table of Contents

1. [Motivation and Design Goals](#1-motivation-and-design-goals)
2. [Exploration vs Exploitation Metrics](#2-exploration-vs-exploitation-metrics)
3. [Search Efficiency Metrics](#3-search-efficiency-metrics)
4. [Component-Level Dynamics](#4-component-level-dynamics)
5. [Posterior Quality Metrics](#5-posterior-quality-metrics)
6. [Diversity and Coverage Metrics](#6-diversity-and-coverage-metrics)
7. [Comparative Metrics (TACTICS vs Legacy)](#7-comparative-metrics-tactics-vs-legacy)
8. [Visualization Recommendations](#8-visualization-recommendations)
9. [Implementation Guide](#9-implementation-guide)
10. [Diagnostic Dashboards](#10-diagnostic-dashboards)

---

## 1. Motivation and Design Goals

### 1.1 Why Metrics Matter

Effective chemical space exploration requires understanding:
1. **Is the algorithm exploring enough?** (avoiding local optima)
2. **Is the algorithm exploiting enough?** (not wasting budget on poor reagents)
3. **Is criticality detection working?** (CATS adaptation effectiveness)
4. **How does TACTICS compare to legacy approaches?** (improvement quantification)
5. **Where are the bottlenecks?** (reagent exhaustion, duplicate pressure)

### 1.2 Design Principles

1. **Interpretability**: Metrics should have clear physical meaning
2. **Actionability**: Metrics should suggest improvements
3. **Efficiency**: Metrics should be computable online during search
4. **Comparability**: Metrics should enable apples-to-apples comparison across runs

### 1.3 Metric Categories

| Category | Purpose | Key Questions | Ground Truth Required? |
|----------|---------|---------------|----------------------|
| **Exploration/Exploitation** | Balance assessment | Are we stuck? Too random? | **No** - computable online |
| **Efficiency** | Budget utilization | How fast do we find top compounds? | **Yes** - needs true optima |
| **Component Dynamics** | Per-component analysis | Which components are critical? | **No** - computable online |
| **Posterior Quality** | Uncertainty calibration | Are posteriors reliable? | **Partial** - some metrics need truth |
| **Diversity** | Chemical space coverage | Are we seeing diverse chemistry? | **No** - computable online |

### 1.4 Ground Truth Requirements

An important distinction for practical use:

**Metrics computable during any search** (no ground truth needed):
- Selection entropy, exploitation ratio, effective temperature
- Criticality time series, CATS multipliers
- Coverage, reagent entropy, scaffold diversity
- Best-so-far trajectory, unique compound rate

**Metrics requiring ground truth** (benchmark/retrospective analysis only):
- Top-K hit rate (needs to know true top-K)
- Cumulative regret (needs to know optimal score)
- Posterior mean correlation (needs true reagent means)
- Budget efficiency (needs to know when all top-K found)

**Metrics with partial requirements**:
- Uncertainty calibration: Can check internal consistency without truth, but full calibration needs true values

---

## 2. Exploration vs Exploitation Metrics

### 2.1 Effective Temperature Distribution

**Definition 2.1 (Effective Temperature)**: Track the distribution of effective temperatures across iterations.

$$T_c^{\text{eff}}(t) = T_c^{\text{base}}(t) \cdot m_c^{\text{eff}}(t)$$

**Metric**: Temperature entropy over time:
$$H_T(t) = -\sum_{c=1}^{C} \frac{T_c^{\text{eff}}(t)}{\sum_{c'} T_{c'}^{\text{eff}}(t)} \ln \frac{T_c^{\text{eff}}(t)}{\sum_{c'} T_{c'}^{\text{eff}}(t)}$$

**Interpretation**:
- Low $H_T$: Temperatures are concentrated (one component dominates)
- High $H_T$: Temperatures are uniform (balanced exploration)

**Visualization**: Time series of $T_c^{\text{eff}}$ for each component, with shaded regions for phases.

### 2.2 Selection Entropy

**Definition 2.2 (Selection Entropy)**: Measure the entropy of reagent selections over a window.

For a sliding window of $W$ selections from component $c$:
$$p_{c,i}^{(W)} = \frac{\#\{r_{c,i} \text{ selected in window}\}}{W}$$

$$H_{\text{sel},c}^{(W)} = -\sum_{i=1}^{n_c} p_{c,i}^{(W)} \ln(p_{c,i}^{(W)} + \epsilon)$$

**Normalized Selection Entropy**:
$$\hat{H}_{\text{sel},c}^{(W)} = \frac{H_{\text{sel},c}^{(W)}}{\ln(n_c)}$$

**Interpretation**:
- $\hat{H}_{\text{sel}} \approx 1$: Uniform selection (maximum exploration)
- $\hat{H}_{\text{sel}} \approx 0$: Concentrated selection (maximum exploitation)

**Recommended Window**: $W = 100$ iterations

### 2.3 Exploitation Ratio

**Definition 2.3 (Exploitation Ratio)**: Fraction of selections from top-$k$% posterior reagents.

Let $\mathcal{T}_c(k)$ be the set of reagents in the top $k$% by posterior mean.

$$\rho_{\text{exploit}}^{(k)}(t) = \frac{\#\{r_{c,i}(t) \in \mathcal{T}_c(k)\}}{C}$$

**Interpretation**:
- Early search: $\rho_{\text{exploit}}^{(10)} \approx 0.1$ (random)
- Late search: $\rho_{\text{exploit}}^{(10)} \to 1$ (exploitation dominates)

**Recommendation**: Track $k \in \{10, 25, 50\}$ to see exploitation at different thresholds.

### 2.4 Effective Exploration Bonus

**Definition 2.4**: For Bayes-UCB, track the average exploration bonus:
$$B_{\text{explore}}(t) = \frac{1}{|S_t|} \sum_{i \in S_t} \frac{\sigma_i \cdot t_{\nu_i}(p^{\text{eff}})}{\sqrt{n_i}}$$

where $S_t$ is the set of selected reagents at iteration $t$.

**For RWS**, approximate via temperature effect:
$$B_{\text{explore}}^{\text{RWS}}(t) = \frac{1}{C} \sum_{c=1}^{C} T_c^{\text{eff}}(t)$$

---

## 3. Search Efficiency Metrics

### 3.1 Top-K Hit Rate

**Definition 3.1 (Top-K Hit Rate)**: Fraction of true top-$K$ compounds found by iteration $t$.

Let $\mathcal{G}_K$ be the global top-$K$ compounds (if known from exhaustive evaluation).

$$\text{HitRate}_K(t) = \frac{|\mathcal{E}_t \cap \mathcal{G}_K|}{K}$$

where $\mathcal{E}_t$ is the set of evaluated compounds up to iteration $t$.

**Use Case**: Benchmark comparison when ground truth is available.

### 3.2 Cumulative Regret

**Definition 3.2 (Cumulative Regret)**:
$$R_t = \sum_{\tau=1}^{t} (f^* - f(\mathbf{r}^{(\tau)}))$$

For minimization mode:
$$R_t = \sum_{\tau=1}^{t} (f(\mathbf{r}^{(\tau)}) - f_{\min})$$

**Visualization**: Plot $R_t / t$ (average regret) over time. Should decrease.

### 3.3 Best-So-Far Trajectory

**Definition 3.3 (Best-So-Far)**:
$$f^{\text{best}}(t) = \begin{cases}
\max_{\tau \leq t} f(\mathbf{r}^{(\tau)}) & \text{maximize} \\
\min_{\tau \leq t} f(\mathbf{r}^{(\tau)}) & \text{minimize}
\end{cases}$$

**Metric**: Time to reach $\alpha$% of optimal:
$$T_\alpha = \min\{t : f^{\text{best}}(t) \geq \alpha \cdot f^*\}$$

### 3.4 Budget Efficiency

**Definition 3.4 (Budget Efficiency)**: Compounds evaluated to find top-$K$:
$$\text{Efficiency}_K = \frac{K}{T_{\text{HitRate}_K = 1}}$$

**Interpretation**: Higher is better. Efficiency = 1 means perfect oracle.

### 3.5 Unique Compound Rate

**Definition 3.5**: Under DisallowTracker, this should be 100%. Track as sanity check:
$$\text{UniqueRate}(t) = \frac{|\mathcal{E}_t|}{t}$$

For legacy (without DisallowTracker):
$$\text{DuplicateRate}(t) = 1 - \text{UniqueRate}(t)$$

---

## 4. Component-Level Dynamics

### 4.1 Criticality Time Series

**Definition 4.1**: Track component criticality over time using z-score softmax with IPR:
$$\kappa_c(t) = 1 - \frac{1/\text{IPR}_c(t)}{n_c^{\text{active}}}$$

where $\text{IPR}_c = \sum_i p_{c,i}^2$ and $p_{c,i} = \text{softmax}(z_{c,i})$ are SNR-dampened, N-sharpened z-score probabilities. See [thompson_sampling_equations.md, Section 5.1](thompson_sampling_equations.md#51-component-criticality-via-z-score-softmax-with-ipr) for the full derivation.

> **Note: Dual Role of Criticality**
>
> Criticality ($\kappa_c$) serves two distinct purposes in TACTICS:
>
> 1. **Algorithm component**: Drives CATS temperature/percentile modulation (see [thompson_sampling_equations.md, Section 5.1](thompson_sampling_equations.md))
> 2. **Diagnostic metric**: Reveals component convergence dynamics for monitoring
>
> When tracking criticality as a metric, you are observing the same quantity that CATS uses internally to make decisions. This provides direct insight into *why* CATS is adjusting exploration for each component.
>
> The `track_diagnostics=True` config option enables automatic collection of criticality and all intermediate values (SNR, IPR, effective_n, multipliers) per cycle. Use `sampler.get_diagnostics()` to retrieve the collected data as a Polars DataFrame.

**Visualization**: Stacked area chart or multi-line plot of $\kappa_c(t)$ for all components.

**Interpretation**:
- Rising $\kappa_c$: Component is converging (posterior sharpening) → CATS will reduce exploration
- Stable low $\kappa_c$: Component remains flexible → CATS will maintain/increase exploration
- Oscillating $\kappa_c$: Unstable criticality (may indicate insufficient data or CATS phase transitions)

### 4.2 Component Selection Frequency

**Definition 4.2**: Track how often each reagent is selected:
$$f_{c,i}(t) = \#\{r_{c,i} \text{ selected in } [0, t]\}$$

**Derived Metrics**:
- Gini coefficient of $f_{c,i}(t)$: Measures selection inequality
- Top-1 concentration: $\max_i f_{c,i}(t) / t$

### 4.3 CATS Multiplier Distribution

**Definition 4.3**: Distribution of CATS multipliers across components, using the relative neutral-point mapping (see [thompson_sampling_equations.md, Section 5.2.2](thompson_sampling_equations.md#522-relative-neutral-point-multiplier)):

$$m_c(t) = \begin{cases}
m_{\max} + \frac{\kappa_c}{\bar{\kappa}} \cdot (1 - m_{\max}) & \text{if } \kappa_c \leq \bar{\kappa} \\[4pt]
1 + \frac{\kappa_c - \bar{\kappa}}{1 - \bar{\kappa}} \cdot (m_{\min} - 1) & \text{if } \kappa_c > \bar{\kappa}
\end{cases}$$

where $\bar{\kappa}$ is the mean criticality across components.

**Visualization**: Violin plots of $m_c$ distribution at different phases (warmup, early, mid, late).

### 4.4 Thermal Cycling Effectiveness

**Definition 4.4**: Measure the "surprise" when component is heated vs cooled:
$$\Delta_c(t) = \mathbb{E}[\text{score} | c \text{ heated at } t] - \mathbb{E}[\text{score} | c \text{ cooled at } t]$$

**Interpretation**: Large $\Delta_c$ indicates thermal cycling is finding better combinations when exploring component $c$.

### 4.5 Reagent Exhaustion Tracking

**Definition 4.5**: Track exhaustion progress:
$$\text{ExhaustionProgress}_c(t) = \frac{\sum_i |D_{c,i}(t)|}{n_c \cdot E_c}$$

where $D_{c,i}(t)$ is the disallow set for reagent $i$ at component $c$, and $E_c = \prod_{c' \neq c} n_{c'}$.

**Visualization**: Progress bars or time series showing approach to full coverage.

---

## 5. Posterior Quality Metrics

### 5.1 Posterior Mean Correlation

**Definition 5.1**: Compare posterior means to true means (if available):
$$\rho_{\mu}(t) = \text{Corr}(\mu_{c,i}(t), \theta_{c,i}^*)$$

**Use Case**: Validation with synthetic data or retrospective analysis.

### 5.2 Uncertainty Calibration

**Definition 5.2**: Check if posterior uncertainty matches empirical variance.

For reagent $i$ with multiple observations:
$$\text{Z-score}_{c,i} = \frac{\bar{x}_{c,i} - \mu_{c,i}}{\sigma_{c,i} / \sqrt{N_{c,i}}}$$

If well-calibrated, Z-scores should follow $\mathcal{N}(0, 1)$.

**Metric**: Kolmogorov-Smirnov statistic against standard normal.

### 5.3 Posterior Variance Trajectory

**Definition 5.3**: Average posterior variance over time:
$$\bar{\sigma}^2(t) = \frac{1}{\sum_c n_c} \sum_{c,i} \sigma_{c,i}^2(t)$$

**Expected Behavior**: $\bar{\sigma}^2(t) \propto 1/t$ (variance decreases as observations accumulate).

### 5.4 Prior-Posterior KL Divergence

**Definition 5.4**: Measure information gain:
$$D_{KL}(t) = \sum_{c,i} D_{KL}(\text{Prior} \| \text{Posterior}_{c,i}(t))$$

For Gaussians:
$$D_{KL} = \frac{1}{2}\left[\frac{\sigma_0^2}{\sigma_i^2} + \frac{(\mu_i - \mu_0)^2}{\sigma_i^2} - 1 + \ln\frac{\sigma_i^2}{\sigma_0^2}\right]$$

---

## 6. Diversity and Coverage Metrics

### 6.1 Chemical Space Coverage

**Definition 6.1**: Fraction of reagents sampled at least once:
$$\text{Coverage}_c(t) = \frac{|\{i : N_{c,i}(t) > 0\}|}{n_c}$$

**Global Coverage**:
$$\text{Coverage}(t) = \frac{1}{C} \sum_{c=1}^{C} \text{Coverage}_c(t)$$

### 6.2 Pairwise Diversity (Chemical Fingerprints)

**Definition 6.2**: Average Tanimoto dissimilarity of selected compounds:
$$\text{Diversity}(t) = 1 - \frac{2}{|\mathcal{E}_t|(|\mathcal{E}_t|-1)} \sum_{i < j} \text{Tanimoto}(\mathbf{fp}_i, \mathbf{fp}_j)$$

**Implementation**: Use Morgan fingerprints (radius 2, 2048 bits).

### 6.3 Scaffold Diversity

**Definition 6.3**: Number of unique Murcko scaffolds in evaluated set:
$$\text{ScaffoldDiversity}(t) = |\{\text{MurckoScaffold}(\mathbf{r}) : \mathbf{r} \in \mathcal{E}_t\}|$$

**Normalized**: Divide by $|\mathcal{E}_t|$ for scaffold coverage ratio.

### 6.4 Reagent Entropy

**Definition 6.4**: Shannon entropy of reagent selection distribution:
$$H_{\text{reagent}}(t) = -\sum_{c,i} \frac{N_{c,i}(t)}{\sum_{c',j} N_{c',j}(t)} \ln \frac{N_{c,i}(t)}{\sum_{c',j} N_{c',j}(t)}$$

**Maximum**: $H_{\max} = \ln(\sum_c n_c)$ when all reagents sampled equally.

---

## 7. Comparative Metrics (TACTICS vs Legacy)

### 7.1 CATS Improvement Ratio

**Definition 7.1**: Compare CATS to baseline RWS:
$$\text{ImprovementRatio} = \frac{\text{HitRate}_K^{\text{CATS}}(t)}{\text{HitRate}_K^{\text{RWS}}(t)}$$

**Visualization**: Plot over iterations to show when CATS advantage emerges.

### 7.2 Adaptation Quantification

**Definition 7.2**: Measure how much CATS adapts compared to fixed temperature:
$$\text{AdaptationRange}(t) = \max_c T_c^{\text{eff}}(t) - \min_c T_c^{\text{eff}}(t)$$

For legacy RWS: $\text{AdaptationRange} = \alpha - \beta$ (constant).
For CATS: $\text{AdaptationRange}$ varies based on criticality.

### 7.3 Duplicate Avoidance Efficiency

**Definition 7.3**: Compare duplicate rates:
$$\text{DuplicateReduction} = \frac{\text{DuplicateRate}^{\text{Legacy}}(t) - \text{DuplicateRate}^{\text{TACTICS}}(t)}{\text{DuplicateRate}^{\text{Legacy}}(t)}$$

With DisallowTracker: $\text{DuplicateRate}^{\text{TACTICS}} = 0$, so reduction = 100%.

### 7.4 Convergence Speed

**Definition 7.4**: Time to reach convergence threshold:
$$T_{\text{converge}}^\theta = \min\left\{t : \frac{d f^{\text{best}}}{dt} < \theta\right\}$$

**Speedup**: $\text{Speedup} = T_{\text{converge}}^{\text{Legacy}} / T_{\text{converge}}^{\text{TACTICS}}$

### 7.5 Regret Comparison

**Definition 7.5**: Relative regret reduction:
$$\text{RegretReduction}(t) = \frac{R_t^{\text{Legacy}} - R_t^{\text{TACTICS}}}{R_t^{\text{Legacy}}}$$

---

## 8. Visualization Recommendations

### 8.1 Real-Time Dashboard Components

| Metric | Visualization Type | Update Frequency |
|--------|-------------------|------------------|
| Best-So-Far | Line chart | Every iteration |
| Selection Entropy | Line chart (per component) | Every 10 iterations |
| Criticality | Stacked area chart | Every 10 iterations |
| Effective Temperature | Multi-line chart | Every iteration |
| Coverage | Progress bars | Every 100 iterations |
| Posterior Means | Histogram | Every 100 iterations |

### 8.2 Post-Run Analysis Plots

1. **Phase Diagram**: 2D plot with X = exploitation ratio, Y = selection entropy. Color by iteration.

2. **Criticality Heatmap**: Component × Iteration heatmap of $\kappa_c(t)$.

3. **Cumulative Regret Curves**: Compare multiple strategies on same axes.

4. **Fingerprint PCA**: 2D projection of selected compounds, colored by iteration.

5. **Reagent Selection Heatmap**: Reagent × Iteration heatmap of selection frequency.

### 8.3 Diagnostic Alerts

| Alert | Trigger Condition | Suggested Action |
|-------|-------------------|------------------|
| Stalled Exploration | $\hat{H}_{\text{sel}} < 0.1$ for 500 iterations | Increase $\alpha$ |
| Over-Exploration | $\rho_{\text{exploit}}^{(25)} < 0.3$ in late phase | Decrease $\alpha/\beta$ ratio |
| Unbalanced Criticality | $\kappa_c$ variance > 0.2 | Check warmup balance |
| Rapid Exhaustion | ExhaustionProgress > 0.5 | Consider larger library |

---

## 9. Implementation Guide

### 9.1 Data Collection Schema

```python
@dataclass
class IterationMetrics:
    """Metrics collected at each iteration."""
    iteration: int
    timestamp: float

    # Selection
    selected_reagents: List[int]  # indices per component
    score: float
    is_best: bool

    # Component state
    criticalities: List[float]  # per component
    effective_temperatures: List[float]  # per component (or percentiles for Bayes-UCB)
    cats_multipliers: List[float]  # per component

    # Posterior state (sampled)
    mean_posterior_mean: float
    mean_posterior_std: float

    # Exploration/exploitation
    selection_from_top_10pct: List[bool]  # per component

@dataclass
class WindowMetrics:
    """Metrics computed over sliding windows."""
    window_start: int
    window_end: int

    # Selection entropy per component
    selection_entropies: List[float]

    # Coverage
    new_reagents_covered: int
    cumulative_coverage: float

    # Performance
    best_in_window: float
    mean_in_window: float
    std_in_window: float
```

### 9.2 Efficient Computation

**Online Algorithms**:
- Selection entropy: Use reservoir sampling for window
- Coverage: Use bitset for reagent tracking
- Best-so-far: Simple max/min update

```python
class OnlineMetricsTracker:
    def __init__(self, n_components: int, reagent_counts: List[int]):
        self.n_components = n_components
        self.reagent_counts = reagent_counts

        # Coverage tracking
        self.reagent_seen = [np.zeros(n, dtype=bool) for n in reagent_counts]

        # Best-so-far
        self.best_score = -np.inf  # or +np.inf for minimize

        # Selection history (circular buffer for entropy)
        self.window_size = 100
        self.selection_buffer = deque(maxlen=self.window_size)

    def update(self, selected_reagents: List[int], score: float):
        # Update coverage
        for c, idx in enumerate(selected_reagents):
            self.reagent_seen[c][idx] = True

        # Update best
        if score > self.best_score:  # or < for minimize
            self.best_score = score

        # Update selection buffer
        self.selection_buffer.append(selected_reagents)

    def compute_selection_entropy(self, component: int) -> float:
        if len(self.selection_buffer) < self.window_size:
            return np.nan

        counts = np.zeros(self.reagent_counts[component])
        for selection in self.selection_buffer:
            counts[selection[component]] += 1

        probs = counts / counts.sum()
        probs = probs[probs > 0]  # filter zeros
        entropy = -np.sum(probs * np.log(probs))
        return entropy / np.log(self.reagent_counts[component])  # normalize

    def compute_coverage(self, component: int) -> float:
        return self.reagent_seen[component].sum() / self.reagent_counts[component]
```

### 9.3 Storage and Export

**Recommended Format**: Parquet for large runs, JSON for metadata.

```python
def export_run_metrics(
    iteration_metrics: List[IterationMetrics],
    window_metrics: List[WindowMetrics],
    run_config: dict,
    output_path: Path
):
    # Convert to DataFrames
    df_iter = pd.DataFrame([asdict(m) for m in iteration_metrics])
    df_window = pd.DataFrame([asdict(m) for m in window_metrics])

    # Save
    df_iter.to_parquet(output_path / "iteration_metrics.parquet")
    df_window.to_parquet(output_path / "window_metrics.parquet")

    with open(output_path / "config.json", "w") as f:
        json.dump(run_config, f, indent=2)
```

---

## 10. Diagnostic Dashboards

### 10.1 Live Monitoring Dashboard (Marimo/Panel)

**Layout**:
```
┌─────────────────────────────────────────────────────────────┐
│  [Header: Run ID, Elapsed Time, Iterations, Best Score]     │
├─────────────────────┬───────────────────────────────────────┤
│  Best-So-Far Plot   │  Selection Entropy (per component)    │
│  (Line chart)       │  (Multi-line chart)                   │
├─────────────────────┼───────────────────────────────────────┤
│  Criticality        │  Effective Temperature                │
│  (Stacked area)     │  (Multi-line chart)                   │
├─────────────────────┼───────────────────────────────────────┤
│  Coverage Progress  │  Exploitation Ratio                   │
│  (Progress bars)    │  (Gauge charts)                       │
├─────────────────────┴───────────────────────────────────────┤
│  [Alerts Panel: Warnings and Recommendations]                │
└─────────────────────────────────────────────────────────────┘
```

### 10.2 Post-Run Analysis Dashboard

**Tabs**:
1. **Summary**: Key statistics, comparison to baselines
2. **Exploration**: Selection entropy, coverage, diversity
3. **Exploitation**: Hit rates, regret, best trajectory
4. **Components**: Per-component criticality, selection frequency
5. **Posteriors**: Calibration, mean trajectories, uncertainty
6. **Comparison**: Side-by-side with legacy or other runs

### 10.3 Comparative Analysis Tools

**Multi-Run Comparison**:
```python
def compare_runs(
    run_paths: List[Path],
    labels: List[str],
    metrics: List[str] = ["hit_rate_100", "regret", "coverage"]
) -> pd.DataFrame:
    """
    Load and compare multiple runs.

    Returns DataFrame with columns: run_label, iteration, metric_name, metric_value
    """
    all_data = []
    for path, label in zip(run_paths, labels):
        df = pd.read_parquet(path / "iteration_metrics.parquet")
        for metric in metrics:
            all_data.append({
                "run_label": label,
                "iteration": df["iteration"],
                "metric": metric,
                "value": compute_metric(df, metric)
            })
    return pd.DataFrame(all_data)
```

---

## Appendix A: Quick Reference for Developers

### A.1 Metrics to Add to New Strategies

When implementing a new selection strategy, track:
1. Selection probabilities (for debugging)
2. Exploration bonus magnitude
3. Any strategy-specific internal state

### A.2 Suggested Improvements Based on Metrics

| Observed Pattern | Likely Issue | Suggested Fix |
|-----------------|--------------|---------------|
| Low coverage, high exploitation | Too exploitative | Increase temperature/percentile |
| High coverage, low hit rate | Too explorative | Decrease temperature/percentile |
| Uneven criticality | Warmup imbalance | Use balanced warmup |
| Rapid exhaustion | Small library | Increase library size |
| Oscillating criticality | Unstable posteriors | Increase min_observations |

### A.3 Metric Computation Frequency

| Metric Type | Recommended Frequency | Reason |
|-------------|----------------------|--------|
| Per-iteration | Every 1 iteration | Real-time tracking |
| Windowed | Every 10-100 iterations | Noise reduction |
| Expensive (fingerprints) | Every 100-1000 iterations | Computational cost |
| Post-run only | End of run | Requires full data |

---

*Document Version: 2.0*
*Last Updated: March 2026*
*Authors: TACTICS Development Team*
