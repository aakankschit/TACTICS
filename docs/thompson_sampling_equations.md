# Mathematical Foundations of Thompson Sampling in TACTICS

## A Scientific Treatment of Component-Aware Thompson Sampling (CATS) with Derivations

---

## Abstract

This document provides rigorous mathematical derivations for the Thompson Sampling algorithms implemented in TACTICS (Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection). We derive the theoretical foundations of Component-Aware Thompson Sampling (CATS) and compare it systematically with legacy Roulette Wheel Selection (RWS) and standard Thompson Sampling (TS) frameworks. The key innovations of CATS include IPR-based component criticality detection with SNR dampening and N-adaptive sharpening, adaptive temperature modulation with relative neutral-point multipliers, criticality-weighted component rotation, and progressive exploration-to-exploitation transition. We demonstrate how CATS provides principled, automatic tuning of exploration/exploitation trade-offs that legacy approaches require manual intervention to achieve.

---

## Table of Contents

1. [Notation and Preliminaries](#1-notation-and-preliminaries)
2. [Problem Formulation](#2-problem-formulation)
3. [Bayesian Posterior Framework](#3-bayesian-posterior-framework)
4. [Legacy Approaches: Standard TS and RWS](#4-legacy-approaches-standard-ts-and-rws)
5. [Component-Aware Thompson Sampling (CATS)](#5-component-aware-thompson-sampling-cats)
6. [Warmup Strategies: Theoretical Analysis](#6-warmup-strategies-theoretical-analysis)
7. [Duplicate Prevention: DisallowTracker](#7-duplicate-prevention-disallowtracker)
8. [Comparative Analysis: CATS vs Legacy](#8-comparative-analysis-cats-vs-legacy)
9. [Convergence Properties](#9-convergence-properties)
10. [References](#10-references)

---

## 1. Notation and Preliminaries

### 1.1 Symbol Table

| Symbol | Description | Domain |
|--------|-------------|--------|
| $C$ | Number of reaction components | $\mathbb{Z}^+$ |
| $\mathcal{R}_c$ | Set of reagents for component $c$ | Finite set |
| $n_c$ | Number of reagents in component $c$: $|\mathcal{R}_c|$ | $\mathbb{Z}^+$ |
| $r_{c,i}$ | Reagent $i$ in component $c$ | $r_{c,i} \in \mathcal{R}_c$ |
| $\mu_{c,i}$ | Posterior mean for reagent $r_{c,i}$ | $\mathbb{R}$ |
| $\sigma_{c,i}$ | Posterior standard deviation for $r_{c,i}$ | $\mathbb{R}^+$ |
| $\sigma^2_{c,i}$ | Posterior variance: $\sigma^2_{c,i} = \sigma_{c,i}^2$ | $\mathbb{R}^+$ |
| $N_{c,i}$ | Number of observations for $r_{c,i}$ | $\mathbb{Z}^{\geq 0}$ |
| $\sigma_0^2$ | Known (prior) variance | $\mathbb{R}^+$ |
| $\alpha$ | Temperature for heated component | $\mathbb{R}^+$ |
| $\beta$ | Temperature for cooled components | $\mathbb{R}^+$ |
| $t$ | Current iteration/cycle | $\mathbb{Z}^{\geq 0}$ |
| $T$ | Total number of iterations | $\mathbb{Z}^+$ |
| $\gamma$ | Progress fraction: $\gamma = t/T$ | $[0, 1]$ |
| $\kappa_c$ | Component criticality for component $c$ | $[0, 1]$ |
| $H_c$ | Shannon entropy for component $c$ | $\mathbb{R}^{\geq 0}$ |
| $m_c$ | CATS temperature multiplier for component $c$ | $\mathbb{R}^+$ |
| $w(\gamma)$ | Progressive criticality weight | $[0, 1]$ |

### 1.2 Optimization Modes

TACTICS supports two optimization modes:

**Maximization Mode**: Seek compounds with highest scores (e.g., binding affinity, potency)
$$\mathbf{r}^* = \arg\max_{\mathbf{r} \in \mathcal{R}_1 \times \cdots \times \mathcal{R}_C} f(\mathbf{r})$$

**Minimization Mode**: Seek compounds with lowest scores (e.g., docking scores, free energy)
$$\mathbf{r}^* = \arg\min_{\mathbf{r} \in \mathcal{R}_1 \times \cdots \times \mathcal{R}_C} f(\mathbf{r})$$

---

## 2. Problem Formulation

### 2.1 Combinatorial Library Screening as Multi-Armed Bandits

Consider a multi-step chemical reaction with $C$ components. Each component $c \in \{1, \ldots, C\}$ has a pool of $n_c$ reagents $\mathcal{R}_c = \{r_{c,1}, \ldots, r_{c,n_c}\}$. The total library size is:

$$|\mathcal{L}| = \prod_{c=1}^{C} n_c$$

**Example**: A two-component reaction with 130 acids and 3,844 amines yields:
$$|\mathcal{L}| = 130 \times 3844 = 499,720 \text{ compounds}$$

### 2.2 Regret Minimization Objective

The goal is to minimize **cumulative regret** over $T$ iterations:

$$R_T = \sum_{t=1}^{T} \left( f(\mathbf{r}^*) - f(\mathbf{r}^{(t)}) \right)$$

where $\mathbf{r}^*$ is the optimal compound and $\mathbf{r}^{(t)}$ is the compound selected at iteration $t$.

### 2.3 Additive Decomposition Assumption

A key assumption in Thompson Sampling for combinatorial libraries is that compound quality can be approximated as a sum of component contributions:

$$\mathbb{E}[f(\mathbf{r})] \approx \sum_{c=1}^{C} \theta_{c,i_c}$$

where $\theta_{c,i_c}$ is the inherent quality of reagent $r_{c,i_c}$. This enables independent posterior updates per reagent.

---

## 3. Bayesian Posterior Framework

### 3.1 Normal-Normal Conjugate Model

We model each reagent's score distribution as Normal with unknown mean and known variance:

**Prior**:
$$\theta_{c,i} \sim \mathcal{N}(\mu_0, \sigma_0^2)$$

**Likelihood**: Given observation $x$:
$$x | \theta_{c,i} \sim \mathcal{N}(\theta_{c,i}, \sigma_0^2)$$

**Posterior**: After observing $x$, the posterior is:
$$\theta_{c,i} | x \sim \mathcal{N}(\mu', {\sigma'}^2)$$

### 3.2 Derivation of Posterior Update Rules

**Theorem 3.1 (Bayesian Posterior Update)**: Given prior $\mathcal{N}(\mu, \sigma^2)$ and observation $x$ with known variance $\sigma_0^2$, the posterior is $\mathcal{N}(\mu', {\sigma'}^2)$ where:

$$\mu' = \frac{\sigma^2 x + \sigma_0^2 \mu}{\sigma^2 + \sigma_0^2}$$

$${\sigma'}^2 = \frac{\sigma^2 \sigma_0^2}{\sigma^2 + \sigma_0^2}$$

**Proof**: Using Bayes' theorem:

$$p(\theta | x) \propto p(x | \theta) p(\theta)$$

The log-posterior is:
$$\log p(\theta | x) = -\frac{(x - \theta)^2}{2\sigma_0^2} - \frac{(\theta - \mu)^2}{2\sigma^2} + \text{const}$$

Expanding and completing the square:
$$= -\frac{\theta^2}{2}\left(\frac{1}{\sigma_0^2} + \frac{1}{\sigma^2}\right) + \theta\left(\frac{x}{\sigma_0^2} + \frac{\mu}{\sigma^2}\right) + \text{const}$$

Recognizing this as a Normal distribution:
$${\sigma'}^2 = \left(\frac{1}{\sigma_0^2} + \frac{1}{\sigma^2}\right)^{-1} = \frac{\sigma^2 \sigma_0^2}{\sigma^2 + \sigma_0^2}$$

$$\mu' = {\sigma'}^2 \left(\frac{x}{\sigma_0^2} + \frac{\mu}{\sigma^2}\right) = \frac{\sigma^2 x + \sigma_0^2 \mu}{\sigma^2 + \sigma_0^2} \quad \square$$

### 3.3 Precision Interpretation

Define precision as $\tau = 1/\sigma^2$. The update rules become:

$$\tau' = \tau + \tau_0$$
$$\mu' = \frac{\tau \mu + \tau_0 x}{\tau'}$$

This shows that precision accumulates additively with observations.

### 3.4 Asymptotic Behavior

After $N$ observations $\{x_1, \ldots, x_N\}$ from a reagent with true mean $\theta^*$:

$$\sigma_N^2 = \frac{\sigma_0^2}{N + 1} \xrightarrow{N \to \infty} 0$$

$$\mu_N \xrightarrow{N \to \infty} \bar{x} \xrightarrow{\text{LLN}} \theta^*$$

The posterior concentrates around the true mean as observations accumulate.

---

## 4. Legacy Approaches: Standard TS and RWS

### 4.0 Understanding Thermal Cycling: Motivation and Mechanism

Before diving into specific algorithms, it's essential to understand **thermal cycling** - a key mechanism that TACTICS inherits from RWS and extends with CATS.

#### 4.0.1 The Exploration Starvation Problem

In multi-component combinatorial optimization, a dangerous failure mode occurs:

1. **Early luck**: Component A gets a high-scoring reagent early by chance
2. **Exploitation bias**: The algorithm focuses on Component A's "winner"
3. **Neglect**: Components B and C receive less exploration
4. **Suboptimality**: Better reagents in B and C are never discovered

**Example**: In a 2-component reaction (acids × amines), if acid #42 scores well in the first few trials, a greedy algorithm will pair acid #42 with many amines, but never explore whether acid #73 might be even better.

#### 4.0.2 Thermal Cycling Solution

Thermal cycling addresses this by **rotating exploration focus** across components:

```
Iteration 1-100:    Component 1 (acids) is "heated" → more exploration
Iteration 101-200:  Component 2 (amines) is "heated" → more exploration
Iteration 201-300:  Component 1 (acids) is "heated" again
... and so on
```

**Key Insight**: By periodically giving each component enhanced exploration opportunities, thermal cycling ensures no component is permanently starved of exploration budget.

#### 4.0.3 Temperature Assignment

At each iteration, exactly **one component is "heated"** while **all others are "cooled"**:

$$T_c^{\text{base}} = \begin{cases}
\alpha & \text{if } c = c_{\text{hot}} \quad \text{(heated: higher temperature → more random)} \\
\beta & \text{if } c \neq c_{\text{hot}} \quad \text{(cooled: lower temperature → more greedy)}
\end{cases}$$

where typically $\alpha > \beta$ (e.g., $\alpha = 0.1$, $\beta = 0.05$).

**Important**: In the Boltzmann distribution used here, **higher temperature means more randomness** (flatter probability distribution), while **lower temperature means more deterministic** selection (peaked at the best option). This is consistent with statistical mechanics conventions.

### 4.1 Standard Thompson Sampling (Legacy TS)

**Algorithm (Legacy TS)**:
```
For each iteration t:
    For each component c:
        1. For each reagent i: sample θ_{c,i} ~ N(μ_{c,i}, σ_{c,i}²)
        2. Select: i* = argmax_i θ_{c,i}
    Evaluate compound and update posteriors
```

**Limitations**:
1. **Greedy exploitation**: No mechanism to increase exploration when stuck
2. **No thermal cycling**: All components treated identically
3. **No component awareness**: Cannot detect "critical" vs "flexible" components

### 4.2 Enhanced Thompson Sampling with Roulette Wheel Selection (Legacy RWS)

The RWS approach (Zhao et al., 2024) introduces Boltzmann-weighted selection and thermal cycling.

#### 4.2.1 Posterior Sampling

Sample from the posterior:
$$s_{c,i} = \mu_{c,i} + \sigma_{c,i} \cdot z_{c,i}, \quad z_{c,i} \sim \mathcal{N}(0, 1)$$

#### 4.2.2 Score Normalization

Normalize scores to zero mean and unit variance:
$$\bar{s}_c = \frac{1}{n_c} \sum_{i=1}^{n_c} s_{c,i}$$

$$\hat{\sigma}_c = \sqrt{\frac{1}{n_c} \sum_{i=1}^{n_c} (s_{c,i} - \bar{s}_c)^2}$$

$$\tilde{s}_{c,i} = \frac{s_{c,i} - \bar{s}_c}{\hat{\sigma}_c}$$

#### 4.2.3 Thermal Cycling

Apply temperature based on component status (see Section 4.0.3 for motivation):
$$T_c = \begin{cases}
\alpha & \text{if } c = c_{\text{hot}} \quad \text{(this component gets exploration focus)} \\
\beta & \text{otherwise} \quad \text{(these components are exploited)}
\end{cases}$$

Rotate heated component after each cycle (ensures fair exploration across all components):
$$c_{\text{hot}} \leftarrow (c_{\text{hot}} + 1) \mod C$$

**Note**: The rotation period is typically one iteration, meaning each component takes turns being heated in a round-robin fashion.

#### 4.2.4 Boltzmann Distribution

Convert scores to selection probabilities:
$$P(r_{c,i}) = \frac{\exp(\tilde{s}_{c,i} / T_c)}{\sum_{j=1}^{n_c} \exp(\tilde{s}_{c,j} / T_c)}$$

**Derivation of Boltzmann Distribution Properties**:

The partition function is:
$$Z_c = \sum_{j=1}^{n_c} \exp(\tilde{s}_{c,j} / T_c)$$

As $T_c \to 0$ (low temperature):
$$P(r_{c,i}) \to \begin{cases} 1 & \text{if } i = \arg\max_j \tilde{s}_{c,j} \\ 0 & \text{otherwise} \end{cases}$$

As $T_c \to \infty$ (high temperature):
$$P(r_{c,i}) \to \frac{1}{n_c}$$ (uniform distribution)

#### 4.2.5 Boltzmann-Weighted Posterior Updates (Legacy Only)

> **Important**: This section describes a feature of **Legacy RWS only**. Modern CATS uses **uniform Bayesian updates** (Section 3.2) instead. This section is retained for historical context and to explain why CATS made a different choice.

Legacy RWS uses Boltzmann-weighted moving averages for posterior updates:

$$w_i = \exp\left(\frac{x_i}{\sigma_{\text{known}}}\right)$$

$$\mu' = \mu + \frac{w_i}{\sum_j w_j} (x_i - \mu)$$

**How it works**: This weighted average gives higher influence to better-scoring observations, creating "rich get richer" dynamics that accelerate convergence to high-quality reagents.

**Theorem 4.1 (Boltzmann Weighting Effect)**: Under Boltzmann weighting with temperature $T$, the effective sample size $N_{\text{eff}}$ satisfies:
$$N_{\text{eff}} \leq N$$
with equality when all observations are identical.

**Proof**: Define effective sample size as:
$$N_{\text{eff}} = \frac{(\sum_i w_i)^2}{\sum_i w_i^2}$$

By Cauchy-Schwarz: $(\sum w_i)^2 \leq N \sum w_i^2$, giving $N_{\text{eff}} \leq N$. $\square$

**Why CATS Uses Uniform Updates Instead**:

| Aspect | Boltzmann-Weighted (Legacy) | Uniform Bayesian (CATS) |
|--------|----------------------------|-------------------------|
| Convergence speed | Faster (amplifies good scores) | Slower but steadier |
| Posterior calibration | Poor (overconfident) | Good (well-calibrated) |
| Criticality calculation | Unreliable (biased posteriors) | Reliable (true uncertainty) |
| Recovery from bad early samples | Difficult | Natural |

CATS relies on accurate posterior uncertainty for criticality calculation. Boltzmann-weighted updates distort the posterior variance, making criticality estimates unreliable. Therefore, **CATS returns to uniform Bayesian updates** as described in Section 3.2.

#### 4.2.6 Reactive Temperature Adjustment

When sampling efficiency drops below threshold $\eta$:
$$\alpha \leftarrow \alpha + \Delta\alpha$$

This is **reactive**: adjustment occurs only after performance degradation is detected.

**Limitations of Legacy RWS**:
1. **Fixed temperature ratio**: $\alpha/\beta$ is static
2. **Reactive adaptation**: Temperature adjusts only when stuck
3. **No component awareness**: Same temperature policy for all components

---

## 5. Component-Aware Thompson Sampling (CATS)

CATS addresses the limitations of legacy approaches through five innovations:

1. **Component Criticality Detection** via z-score softmax with SNR dampening and IPR metric
2. **Adaptive Temperature Modulation** with relative neutral-point multipliers
3. **Progressive Weighting** with observation-gated (p25) ramp-up
4. **Criticality-Weighted Component Rotation** (flexible components heated more often)
5. **Diagnostics API** for post-hoc SAR assessment and convergence analysis

### 5.0 CATS and Thermal Cycling: Two Orthogonal Mechanisms

**This is the most important conceptual section for understanding CATS.**

A common source of confusion is how CATS interacts with thermal cycling. The key insight is that they operate on **orthogonal dimensions**:

| Mechanism | Question Answered | Dimension |
|-----------|-------------------|-----------|
| **Thermal Cycling** | "Which component gets exploration focus THIS iteration?" | Temporal (rotation) |
| **CATS Multiplier** | "Does this component actually NEED more/less exploration?" | Component-specific (adaptation) |

#### 5.0.1 How They Combine

The effective temperature for any component is computed by **multiplying** the thermal cycling base with the CATS adjustment:

$$T_c^{\text{eff}} = \underbrace{T_c^{\text{base}}}_{\text{From thermal cycling}} \times \underbrace{m_c^{\text{eff}}}_{\text{From CATS}}$$

This means:
- **ALL components** receive a base temperature from thermal cycling (one gets $\alpha$, others get $\beta$)
- **ALL components** have their temperature adjusted by CATS based on their individual criticality
- The adjustments are **independent**: a heated component can have its temperature reduced by CATS, and a cooled component can have its temperature increased

#### 5.0.2 Four Possible Scenarios

Consider a two-component system with $\alpha = 0.1$ (hot), $\beta = 0.05$ (cold), and CATS multipliers ranging from $m_{\min} = 0.5$ to $m_{\max} = 2.0$:

| Component State | Thermal Cycling | CATS Criticality | CATS Action | Net Effect |
|-----------------|-----------------|------------------|-------------|------------|
| Heated + Flexible | $T^{\text{base}} = 0.1$ | $\kappa \approx 0$ → $m \approx 2.0$ | Increase temp | $T^{\text{eff}} = 0.2$ (very explorative) |
| Heated + Critical | $T^{\text{base}} = 0.1$ | $\kappa \approx 1$ → $m \approx 0.5$ | Decrease temp | $T^{\text{eff}} = 0.05$ (more exploitative) |
| Cooled + Flexible | $T^{\text{base}} = 0.05$ | $\kappa \approx 0$ → $m \approx 2.0$ | Increase temp | $T^{\text{eff}} = 0.1$ (compensatory exploration) |
| Cooled + Critical | $T^{\text{base}} = 0.05$ | $\kappa \approx 1$ → $m \approx 0.5$ | Decrease temp | $T^{\text{eff}} = 0.025$ (strong exploitation) |

**Key Observations**:
1. A **heated but critical** component has its exploration bonus *reduced* because CATS detects it doesn't need exploration
2. A **cooled but flexible** component has its exploitation bias *compensated* because CATS detects it needs more exploration
3. The thermal cycling ensures **fair rotation** while CATS ensures **appropriate intensity**

#### 5.0.3 Why Both Mechanisms Are Needed

**Thermal cycling alone** (Legacy RWS) treats all components identically when heated/cooled. This wastes exploration budget on components that have already converged.

**CATS alone** (without thermal cycling) might never give exploration opportunities to components that start with misleadingly good early results.

**Together**, they provide:
- **Fairness**: Every component gets periodic exploration opportunities (thermal cycling)
- **Efficiency**: Exploration intensity matches actual need (CATS)

### 5.1 Component Criticality via Z-Score Softmax with IPR

#### 5.1.1 Motivation

**Definition**: A component is **critical** if only a few reagents have high posterior means (peaked distribution). A component is **flexible** if many reagents have similar posterior means (flat distribution).

**Key Insight**: Critical components should be exploited (their best reagents are clear), while flexible components should be explored (many potentially good options remain).

#### 5.1.2 Z-Score Softmax with SNR Dampening

The criticality calculation uses z-score normalization (scale-invariant) with signal-to-noise ratio (SNR) dampening to avoid acting on noise-dominated differences.

**Step 1: Filter to active reagents and compute z-scores**

Only reagents with $N_{c,i} > 0$ are considered. For maximization mode, let $\mu_{c,i}$ be the posterior means; for minimization mode, negate them.

$$z_{c,i} = \frac{\mu_{c,i} - \bar{\mu}_c}{\sigma_{\mu,c}}$$

where $\bar{\mu}_c = \frac{1}{n_c^{\text{active}}} \sum_i \mu_{c,i}$ and $\sigma_{\mu,c} = \text{std}(\mu_{c,i})$.

**Step 2: Signal-to-noise dampening**

Each reagent's posterior standard error is $\text{SE}_{c,i} = \sigma_{c,i} / \sqrt{N_{c,i}}$. The noise level is:
$$\sigma_{\text{noise}} = \sqrt{\frac{1}{n_c^{\text{active}}} \sum_i \text{SE}_{c,i}^2}$$

The signal-to-noise ratio and dampening factor:
$$\text{SNR}_c = \frac{\sigma_{\mu,c}}{\sigma_{\text{noise}}}, \qquad \lambda_c = 1 - \exp(-\max(\text{SNR}_c - 1, 0))$$

The dampened z-scores: $z_{c,i} \leftarrow z_{c,i} \cdot \lambda_c$

**Interpretation**: When SNR $\leq 1$ (noise dominates), $\lambda_c \approx 0$ and z-scores are suppressed, yielding criticality $\approx 0$. When SNR $> 1$, z-scores are preserved proportionally.

**Step 3: N-adaptive sharpening (IPR mode only)**

For large components, the softmax distributes probability mass broadly, making Shannon entropy insensitive to genuine concentration. A sharpening factor counteracts this:

$$\text{sharpening} = \max(1, \sqrt{\ln(n_c^{\text{active}})}), \qquad z_{c,i} \leftarrow z_{c,i} \cdot \text{sharpening}$$

Applied only when `criticality_metric="ipr"` and `n_adaptive_sharpening=True` and $n_c^{\text{active}} > 2$.

**Step 4: Softmax probabilities**

$$p_{c,i} = \frac{\exp(z_{c,i} - \max_j z_{c,j})}{\sum_{j} \exp(z_{c,j} - \max_j z_{c,j})}$$

(The subtraction of $\max_j z_{c,j}$ is for numerical stability.)

#### 5.1.3 Criticality Metrics: IPR vs Shannon Entropy

**Inverse Participation Ratio (IPR)** — default metric:

$$\text{IPR}_c = \sum_{i=1}^{n_c^{\text{active}}} p_{c,i}^2, \qquad N_{\text{eff}} = \frac{1}{\text{IPR}_c}$$

$$\kappa_c^{\text{IPR}} = 1 - \frac{N_{\text{eff}}}{n_c^{\text{active}}}$$

IPR measures probability concentration directly. $N_{\text{eff}}$ is the effective number of reagents with significant probability mass. When one reagent dominates, $N_{\text{eff}} \to 1$ and $\kappa \to 1 - 1/n_c \approx 1$. When all reagents are equal, $N_{\text{eff}} = n_c$ and $\kappa = 0$.

**Shannon Entropy** — legacy metric (`criticality_metric="shannon"`):

$$\kappa_c^{\text{Shannon}} = 1 - \frac{H_c}{\ln(n_c^{\text{active}})}, \qquad H_c = -\sum_i p_{c,i} \ln(p_{c,i})$$

**Why IPR is preferred**: Shannon entropy's logarithmic compression causes it to saturate near 1.0 at large $n_c$ even when genuine concentration exists. IPR scales naturally with the number of dominant reagents, providing more dynamic range for CATS to differentiate components.

#### 5.1.4 Interpretation of Criticality

| Criticality $\kappa_c$ | Concentration | Posterior Shape | Interpretation | Action |
|------------------------|---------------|-----------------|----------------|--------|
| $\approx 0$ | Low | Uniform | Flexible (many options) | **Explore** |
| $\approx 0.5$ | Moderate | Some structure | Balanced | Moderate |
| $\approx 1$ | High | Peaked | Critical (few options) | **Exploit** |

#### 5.1.5 Observation Guard: p25 Percentile Gating

**Problem**: Early in the search, posteriors may be unreliable due to few observations. With large components, a handful of rarely-sampled reagents can pin the minimum observation count at the warmup level forever.

**Solution**: Use the 25th percentile of observation counts (p25) among active reagents rather than the minimum:

$$\kappa_c = \begin{cases}
\kappa_c^{\text{metric}} & \text{if } P_{25}(\{N_{c,i} : N_{c,i} > 0\}) \geq N_{\min} \\
0.5 & \text{otherwise}
\end{cases}$$

where $N_{\min} = 5$ by default. The p25 threshold reflects the bulk of the observation distribution, allowing CATS to ramp up as the majority of reagents accumulate data, rather than being held hostage by the long tail of rarely-sampled reagents.

### 5.2 CATS Temperature Multiplier

#### 5.2.1 Multiplier Range Derivation

The CATS multiplier range is derived from the $\alpha/\beta$ ratio:

$$m_{\max} = \frac{\alpha}{\beta}$$ (for flexible components: increase exploration)

$$m_{\min} = \frac{\beta}{\alpha}$$ (for critical components: increase exploitation)

**Derivation**:
- When $m_c = m_{\max}$: $T_c^{\text{eff}} = \beta \cdot m_{\max} = \beta \cdot (\alpha/\beta) = \alpha$
- When $m_c = m_{\min}$: $T_c^{\text{eff}} = \alpha \cdot m_{\min} = \alpha \cdot (\beta/\alpha) = \beta$

This ensures the effective temperature spans $[\beta, \alpha]$ regardless of thermal cycling state.

#### 5.2.2 Relative Neutral-Point Multiplier

**Problem**: With IPR + N-adaptive sharpening, online criticalities tend to cluster in a narrow range (e.g., [0.6, 1.0]) rather than spanning [0, 1]. A fixed linear mapping with neutral point at $\kappa = 0.5$ would place all components in the exploit zone, nullifying the exploration boost entirely.

**Solution**: Use the **mean criticality across components** as the neutral point. Components below the mean get an exploration boost (multiplier > 1); components above the mean get exploitation reduction (multiplier < 1). This ensures CATS always differentiates between components regardless of the absolute criticality scale.

**Definition 5.3 (Relative Neutral-Point CATS Multiplier)**:

Let $\bar{\kappa} = \frac{1}{C} \sum_{c=1}^{C} \kappa_c$ be the mean criticality across all components (updated each cycle by `rotate_component_weighted`).

$$m_c = \begin{cases}
m_{\max} + \frac{\kappa_c}{\bar{\kappa}} \cdot (1 - m_{\max}) & \text{if } \kappa_c \leq \bar{\kappa} \quad \text{(below mean → explore)} \\[6pt]
1 + \frac{\kappa_c - \bar{\kappa}}{1 - \bar{\kappa}} \cdot (m_{\min} - 1) & \text{if } \kappa_c > \bar{\kappa} \quad \text{(above mean → exploit)}
\end{cases}$$

**Derivation of Mapping**:
- When $\kappa_c = 0$ (most flexible): $m_c = m_{\max} + 0 = m_{\max}$
- When $\kappa_c = \bar{\kappa}$ (at mean): $m_c = m_{\max} + 1 \cdot (1 - m_{\max}) = 1.0$ (neutral)
- When $\kappa_c = 1$ (most critical): $m_c = 1 + 1 \cdot (m_{\min} - 1) = m_{\min}$

This piecewise linear mapping pivots around the mean criticality, providing adaptive differentiation even when all criticalities are high.

### 5.3 Progressive Criticality Weight

#### 5.3.1 Three-Phase Progression

**Problem**: Early criticality estimates are noisy due to limited data.

**Solution**: Gradually introduce criticality-based adjustments:

**Definition 5.4 (Progressive Weight)**:
$$w(\gamma) = \begin{cases}
0 & \text{if } \gamma < \gamma_1 \quad \text{(Phase 1: Pure exploration)} \\
\frac{\gamma - \gamma_1}{\gamma_2 - \gamma_1} & \text{if } \gamma_1 \leq \gamma < \gamma_2 \quad \text{(Phase 2: Transition)} \\
1 & \text{if } \gamma \geq \gamma_2 \quad \text{(Phase 3: Full CATS)}
\end{cases}$$

where:
- $\gamma = t/T$ is the progress fraction
- $\gamma_1 = 0.20$ (20% of iterations)
- $\gamma_2 = 0.60$ (60% of iterations)

#### 5.3.2 Effective Multiplier with Progressive Blending

$$m_c^{\text{eff}} = 1 + w(\gamma) \cdot (m_c - 1)$$

**Derivation**:
- When $w(\gamma) = 0$: $m_c^{\text{eff}} = 1$ (no CATS adjustment)
- When $w(\gamma) = 1$: $m_c^{\text{eff}} = 1 + (m_c - 1) = m_c$ (full CATS)

This blends between neutral ($m=1$) and CATS-adjusted multipliers.

### 5.4 Criticality-Weighted Component Rotation

#### 5.4.1 Motivation

Legacy RWS uses **round-robin** component rotation, cycling the heated component index sequentially: $c_{\text{hot}} \leftarrow (c_{\text{hot}} + 1) \mod C$. This gives every component equal heating time regardless of need.

**Problem**: Components with low criticality (flexible, many viable options) benefit more from exploration than components with high criticality (already converged). Equal heating wastes budget on converged components.

**Solution**: Select the next heated component **probabilistically**, weighted by flexibility (inverse of criticality).

#### 5.4.2 Weighted Rotation Algorithm

**Definition 5.5 (Criticality-Weighted Rotation)**:

Given per-component criticalities $\{\kappa_1, \ldots, \kappa_C\}$:

1. Compute flexibility scores with a minimum floor:
$$f_c = \max(1 - \kappa_c, 0.1)$$

2. Normalize to heating probabilities:
$$P(\text{heat } c) = \frac{f_c}{\sum_{c'} f_{c'}}$$

3. Sample the next heated component:
$$c_{\text{hot}} \sim \text{Categorical}(P(\text{heat } 1), \ldots, P(\text{heat } C))$$

4. Cache the mean criticality for the relative neutral-point multiplier:
$$\bar{\kappa} \leftarrow \frac{1}{C} \sum_c \kappa_c$$

**Properties**:
- Flexible components (low $\kappa$) are heated more often → more exploration where needed
- Critical components (high $\kappa$) are heated less often → less wasted exploration
- The minimum floor of 0.1 ensures no component is completely starved of exploration
- Mean criticality is updated each cycle, keeping the CATS multiplier neutral point current

### 5.5 Complete CATS Equations

#### 5.4.1 Base Temperature from Thermal Cycling

$$T_c^{\text{base}} = \begin{cases}
\alpha & \text{if } c = c_{\text{hot}} \\
\beta & \text{otherwise}
\end{cases}$$

#### 5.4.2 Effective Temperature

$$T_c^{\text{eff}} = T_c^{\text{base}} \cdot m_c^{\text{eff}}$$

#### 5.4.3 Selection Probability

$$P(r_{c,i}) = \frac{\exp(\tilde{s}_{c,i} / T_c^{\text{eff}})}{\sum_{j=1}^{n_c} \exp(\tilde{s}_{c,j} / T_c^{\text{eff}})}$$

### 5.6 Worked Example

**Setup**: Two-component reaction (acids, amines), $\alpha = 0.1$, $\beta = 0.05$.

**Derived parameters**:
- $m_{\max} = \alpha / \beta = 0.1 / 0.05 = 2.0$
- $m_{\min} = \beta / \alpha = 0.05 / 0.1 = 0.5$

**At iteration $t = 4000$ of $T = 10000$** ($\gamma = t/T = 0.4$):

Since $\gamma_1 = 0.20$ and $\gamma_2 = 0.60$, we are in Phase 2 (transition):
$$w(\gamma) = \frac{\gamma - \gamma_1}{\gamma_2 - \gamma_1} = \frac{0.4 - 0.2}{0.6 - 0.2} = \frac{0.2}{0.4} = 0.5$$

**For acids component** with $\kappa_{\text{acid}} = 0.8$ (critical - few good options identified):

Step 1: Compute CATS multiplier from criticality:
$$m_{\text{acid}} = m_{\min} + (m_{\max} - m_{\min}) \cdot (1 - \kappa_{\text{acid}})$$
$$m_{\text{acid}} = 0.5 + (2.0 - 0.5) \cdot (1 - 0.8) = 0.5 + 1.5 \cdot 0.2 = 0.5 + 0.3 = 0.8$$

Step 2: Apply progressive blending (partial CATS effect at $w = 0.5$):
$$m_{\text{acid}}^{\text{eff}} = 1 + w \cdot (m_{\text{acid}} - 1) = 1 + 0.5 \cdot (0.8 - 1) = 1 + 0.5 \cdot (-0.2) = 0.9$$

Step 3: Compute effective temperature (assuming acids are currently heated):
$$T_{\text{acid}}^{\text{eff}} = T_{\text{acid}}^{\text{base}} \times m_{\text{acid}}^{\text{eff}} = \alpha \times 0.9 = 0.1 \times 0.9 = 0.09$$

**Interpretation**: Even though acids are heated ($T^{\text{base}} = 0.1$), CATS *reduces* the effective temperature to 0.09 because the component is critical and doesn't need as much exploration.

---

**For amines component** with $\kappa_{\text{amine}} = 0.2$ (flexible - many potentially good options):

Step 1: Compute CATS multiplier:
$$m_{\text{amine}} = 0.5 + (2.0 - 0.5) \cdot (1 - 0.2) = 0.5 + 1.5 \cdot 0.8 = 0.5 + 1.2 = 1.7$$

Step 2: Apply progressive blending:
$$m_{\text{amine}}^{\text{eff}} = 1 + 0.5 \cdot (1.7 - 1) = 1 + 0.5 \cdot 0.7 = 1.35$$

Step 3: Compute effective temperature (assuming amines are currently cooled):
$$T_{\text{amine}}^{\text{eff}} = T_{\text{amine}}^{\text{base}} \times m_{\text{amine}}^{\text{eff}} = \beta \times 1.35 = 0.05 \times 1.35 = 0.0675$$

**Interpretation**: Even though amines are cooled ($T^{\text{base}} = 0.05$), CATS *increases* the effective temperature to 0.0675 because the component is flexible and needs more exploration.

---

**Summary of This Iteration**:
| Component | $T^{\text{base}}$ | $\kappa$ | $m^{\text{eff}}$ | $T^{\text{eff}}$ | Net Effect |
|-----------|-------------------|----------|------------------|------------------|------------|
| Acids (heated) | 0.10 | 0.8 (critical) | 0.9 | 0.09 | Exploration reduced |
| Amines (cooled) | 0.05 | 0.2 (flexible) | 1.35 | 0.0675 | Exploration increased |

This demonstrates how CATS **rebalances** the exploration budget based on component needs, even overriding the thermal cycling direction when appropriate.

### 5.7 CATS Algorithm Summary

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        CATS ALGORITHM SUMMARY                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  FOR EACH ITERATION t:                                                       │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ STEP 1: CRITICALITY-WEIGHTED ROTATION (assigns base temperature)     │    │
│  │                                                                      │    │
│  │   • Compute flexibility: f_c = max(1 - κ_c, 0.1)                     │    │
│  │   • Sample heated component: c_hot ~ Categorical(f / sum(f))          │    │
│  │   • ONE component is "heated": T_base = α (more exploration)         │    │
│  │   • ALL OTHER components are "cooled": T_base = β (more exploitation)│    │
│  │   • Flexible components get heated MORE OFTEN                         │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                              ↓                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ STEP 2: CRITICALITY DETECTION (for EACH component independently)     │    │
│  │                                                                      │    │
│  │   z-scores → SNR dampening → softmax → IPR → criticality              │    │
│  │   κ_c = 1 - (1/IPR_c) / N   where IPR_c = Σ p_{c,i}²               │    │
│  │                                                                      │    │
│  │   • κ ≈ 0: "Flexible" - many good options, NEEDS exploration         │    │
│  │   • κ ≈ 1: "Critical" - few good options, should EXPLOIT             │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                              ↓                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ STEP 3: RELATIVE NEUTRAL-POINT MULTIPLIER                            │    │
│  │                                                                      │    │
│  │   Neutral point = mean criticality across components (κ̄)             │    │
│  │   κ < κ̄ → m interpolates [m_max, 1.0] → INCREASE exploration        │    │
│  │   κ > κ̄ → m interpolates [1.0, m_min] → DECREASE exploration        │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                              ↓                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ STEP 4: PROGRESSIVE BLENDING (phases in CATS effect gradually)       │    │
│  │                                                                      │    │
│  │   γ = t/T (progress through search)                                  │    │
│  │                                                                      │    │
│  │   Phase 1 (γ < 0.20):  w = 0    → No CATS effect (pure exploration)  │    │
│  │   Phase 2 (0.20-0.60): w = linear → Gradual introduction             │    │
│  │   Phase 3 (γ ≥ 0.60):  w = 1    → Full CATS effect                   │    │
│  │                                                                      │    │
│  │   m_eff = 1 + w × (m_c - 1)                                          │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                              ↓                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ STEP 5: FINAL EFFECTIVE TEMPERATURE                                  │    │
│  │                                                                      │    │
│  │   T_eff = T_base × m_eff                                             │    │
│  │                                                                      │    │
│  │   This temperature is used in the Boltzmann selection distribution   │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
│  ═══════════════════════════════════════════════════════════════════════    │
│                                                                              │
│  KEY INSIGHT: Rotation and CATS are ORTHOGONAL mechanisms:                   │
│                                                                              │
│    • Weighted rotation: WHICH component is heated (flexible ones more often) │
│    • CATS multiplier: HOW MUCH to adjust temperature (relative to mean κ)   │
│                                                                              │
│  A heated component can have its temperature REDUCED by CATS if critical.    │
│  A cooled component can have its temperature INCREASED by CATS if flexible.  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 6. Warmup Strategies: Theoretical Analysis

### 6.1 Role of Warmup

The warmup phase initializes posterior distributions before the main search. Choice of warmup strategy affects:
1. **Posterior quality**: How well initial estimates reflect true reagent quality
2. **Posterior balance**: Whether all reagents have similar confidence levels
3. **Chemical space coverage**: Whether partners are sampled diversely

### 6.2 Standard Warmup

**Algorithm**: Random partner selection with replacement

**Total evaluations**:
$$E_{\text{standard}} = \sum_{c=1}^{C} n_c \times K$$

**Observations per reagent**:
$$N_{c,i} = K \quad \forall c, i$$

**Expected duplicate partners** (birthday problem):
$$\mathbb{E}[\text{duplicates}] = K - n_{c'} \left(1 - \left(1 - \frac{1}{n_{c'}}\right)^K\right)$$

**Approximate for large $n_{c'}$**: Using $(1 - 1/n)^K \approx e^{-K/n}$:
$$\mathbb{E}[\text{duplicates}] \approx K - n_{c'}(1 - e^{-K/n_{c'}})$$

### 6.3 Enhanced Warmup (Stochastic Parallel Pairing)

**Repetition factor**:
$$\rho_c = \left\lceil \frac{n_{\max}}{n_c} \right\rceil$$

**Observations per reagent (IMBALANCED)**:
$$N_{c,i} = \begin{cases}
K \times \rho_c & \text{if } n_c < n_{\max} \\
K & \text{if } n_c = n_{\max}
\end{cases}$$

**Imbalance ratio**:
$$\frac{N_{\text{small}}}{N_{\text{large}}} = \rho_{\text{small}} = \left\lceil \frac{n_{\max}}{n_{\min}} \right\rceil$$

**Problem**: This creates asymmetric posterior confidence:
$$\sigma_{\text{small}}^2 \propto \frac{1}{N_{\text{small}}} \ll \sigma_{\text{large}}^2 \propto \frac{1}{N_{\text{large}}}$$

### 6.4 Balanced Warmup (Stratified Sampling)

**Algorithm**: Exactly $K$ observations per reagent with stratified partner selection

**Stratification**: For partner component $c'$ with $n_{c'}$ reagents, divide into $K$ strata:
$$\text{Stratum}_k = \left[ \left\lfloor \frac{k \cdot n_{c'}}{K} \right\rfloor, \left\lfloor \frac{(k+1) \cdot n_{c'}}{K} \right\rfloor \right)$$

**Partner selection**:
$$j \sim \text{Uniform}(\text{Stratum}_k)$$

**Properties**:
- Balanced: $N_{c,i} = K \quad \forall c, i$
- Stratified: Partners cover full range of partner space
- No duplicates within strata

### 6.5 Per-Reagent Variance with James-Stein Shrinkage

**Definition 6.1 (James-Stein Shrinkage Estimator)**: Given per-reagent sample variance $s_{c,i}^2$ and global variance $s_{\text{global}}^2$:

$$\hat{\sigma}_{c,i}^2 = \frac{N_{c,i} \cdot s_{c,i}^2 + \lambda \cdot s_{\text{global}}^2}{N_{c,i} + \lambda}$$

where $\lambda$ is the shrinkage strength (default: 3.0).

**Derivation**: This is a weighted average with weight:
$$w = \frac{N_{c,i}}{N_{c,i} + \lambda}$$

- When $N_{c,i}$ is large: $w \to 1$, use per-reagent estimate
- When $N_{c,i}$ is small: $w \to 0$, use global estimate

**Theorem 6.1 (Shrinkage Improvement)**: The James-Stein estimator dominates the maximum likelihood estimator in terms of expected squared error when the number of reagents $n_c \geq 3$.

### 6.6 Warmup Strategy Comparison

| Property | Standard | Enhanced | Balanced |
|----------|----------|----------|----------|
| Balance | $\checkmark$ Equal $N_{c,i}$ | $\times$ Imbalanced | $\checkmark$ Equal $N_{c,i}$ |
| Diversity | $\times$ Random | $\checkmark$ Shuffled | $\checkmark$ Stratified |
| Coverage | Variable | Excellent for small | Uniform across strata |
| Per-reagent variance | $\times$ | $\times$ | $\checkmark$ With shrinkage |
| Recommended | Legacy | Small component focus | **General use** |

---

## 7. Duplicate Prevention: DisallowTracker

### 7.1 The Duplicate Sampling Problem

**Problem**: As posteriors sharpen, selection converges to the same "best" reagents, generating duplicate compounds that waste evaluation budget.

**Legacy Approach (Reactive)**:
```
Generate 100 combinations → Filter duplicates → Evaluate unique only
```

**Issues**:
1. Wasted computation (generate then discard)
2. Convergence collapse (>90% duplicates possible)
3. Early termination when stuck

### 7.2 DisallowTracker: Proactive Prevention

**Data Structure**:
$$\text{DisallowMask}: \text{PartialCombination} \to 2^{\text{ReagentIndex}}$$

**Example** (3-component reaction): After sampling $[4, 5, 3]$:
- Key $(*, 5, 3) \mapsto \{4\}$
- Key $(4, *, 3) \mapsto \{5\}$
- Key $(4, 5, *) \mapsto \{3\}$

### 7.3 Integration with Selection

The disallow mask is passed to selection strategies:

**Roulette Wheel**:
```python
probs[disallow_mask] = 0
probs = probs / sum(probs)  # Renormalize
```

**Greedy/UCB**:
```python
scores[disallow_mask] = -inf  # (or +inf for minimize)
return argmax(scores)
```

### 7.4 Exhaustion Handling

**Exhaustion count for reagent at position $c$**:
$$E_c = \prod_{c' \neq c} n_{c'}$$

When a reagent has been paired with all possible partners ($|D_c(r)| = E_c$), it is marked exhausted and automatically excluded.

---

## 8. Comparative Analysis: CATS vs Legacy

### 8.1 Key Differences Summary

| Aspect | Standard TS | Legacy RWS | CATS |
|--------|-------------|------------|------|
| **Selection** | Argmax on samples | Roulette wheel | Roulette wheel |
| **Temperature** | N/A | Fixed $\alpha$, $\beta$ | Adaptive via $m_c$ |
| **Component awareness** | None | None | IPR-based criticality |
| **Exploration control** | None | Fixed thermal cycling | Criticality-based |
| **Adaptation** | None | Reactive (when stuck) | Proactive (continuous) |
| **Posterior update** | Uniform Bayesian | Boltzmann-weighted | Uniform Bayesian |
| **Duplicate handling** | None | Reactive filtering | Proactive DisallowTracker |

### 8.2 Mathematical Comparison

**Legacy RWS Temperature**:
$$T_c^{\text{RWS}} = \begin{cases} \alpha & c = c_{\text{hot}} \\ \beta & \text{otherwise} \end{cases}$$

**CATS Effective Temperature**:
$$T_c^{\text{CATS}} = T_c^{\text{base}} \cdot \left[1 + w(\gamma)(m_c - 1)\right]$$

**Key insight**: CATS modulates temperature based on **component-specific information** (criticality), while RWS uses a **global, fixed policy**.

### 8.3 Advantages of CATS

1. **Automatic Tuning**: No need to manually adjust temperatures during search
2. **Component-Specific Adaptation**: Different components can have different exploration levels
3. **Progressive Introduction**: Avoids early-iteration instability
4. **Principled Foundation**: Information-theoretic basis (IPR criticality with SNR dampening)

### 8.4 When CATS Outperforms Legacy

| Scenario | CATS Advantage |
|----------|----------------|
| Heterogeneous components | Detects critical vs flexible components |
| Long searches | Progressive adaptation improves late-stage |
| Multi-component reactions | Per-component optimization |
| Unknown score landscape | Entropy-based adaptation is data-driven |

---

## 9. Convergence Properties

### 9.1 Posterior Concentration

**Theorem 9.1**: Under the CATS algorithm with properly initialized posteriors, the posterior variance for any reagent $r_{c,i}$ satisfies:
$$\sigma_{c,i}^2(N) = \frac{\sigma_0^2}{N + 1} = O(1/N)$$

### 9.2 Exploration Guarantee

**Theorem 9.2 (Boltzmann Exploration)**: For any $\epsilon > 0$ and temperature $T > 0$, every reagent has selection probability at least:
$$P(r_{c,i}) \geq \frac{\exp(-\|\tilde{s}\|_\infty / T)}{n_c}$$

where $\|\tilde{s}\|_\infty = \max_j |\tilde{s}_{c,j}|$.

### 9.3 Criticality Convergence

**Theorem 9.3**: As the posterior means converge to true values, criticality $\kappa_c$ converges to the true entropy-based measure of the score distribution for component $c$.

---

## 10. Tunable Parameters Reference

This section provides a complete reference for all tunable parameters in the TACTICS Thompson Sampling implementation, organized by configuration class.

### 10.1 Main Thompson Sampling Configuration

**Class**: `ThompsonSamplingConfig`

| Parameter | Type | Default | Range/Constraints | Description |
|-----------|------|---------|-------------------|-------------|
| `num_ts_iterations` | int | *required* | > 0 | Total number of Thompson Sampling iterations in the search phase |
| `num_warmup_trials` | int | 3 | > 0 | Number of warmup observations per reagent |
| `batch_size` | int | 1 | > 0 | Compounds to sample per iteration (for parallel evaluation) |
| `max_resamples` | int | None | > 0 or None | Maximum resampling attempts for duplicates (None = unlimited) |
| `use_boltzmann_weighting` | bool | False | True/False | Use legacy Boltzmann-weighted Bayesian updates (not recommended with CATS) |
| `track_diagnostics` | bool | False | True/False | Collect per-cycle diagnostics (criticality, SNR, multipliers) for post-hoc analysis |
| `hide_progress` | bool | False | True/False | Hide progress bars during execution |

**Example**:
```python
from TACTICS.thompson_sampling import ThompsonSamplingConfig

config = ThompsonSamplingConfig(
    synthesis_pipeline=pipeline,
    num_ts_iterations=5000,
    num_warmup_trials=5,
    batch_size=10,
    strategy_config=...,
    evaluator_config=...,
)
```

### 10.2 Roulette Wheel Selection (RWS/CATS) Parameters

**Class**: `RouletteWheelConfig`

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `mode` | str | "maximize" | "maximize", "minimize", "maximize_boltzmann", "minimize_boltzmann" | Optimization direction |
| `alpha` | float | 0.1 | > 0 | Base temperature for **heated** component (higher = more exploration) |
| `beta` | float | 0.05 | > 0 | Base temperature for **cooled** components (lower = more exploitation) |
| `exploration_phase_end` | float | 0.20 | (0, 1] | Fraction of iterations before CATS starts ($\gamma_1$) |
| `transition_phase_end` | float | 0.60 | (0, 1] | Fraction of iterations when CATS is fully active ($\gamma_2$) |
| `min_observations` | int | 5 | > 0 | Minimum observations per reagent before trusting criticality ($N_{\min}$) |
| `criticality_metric` | str | "ipr" | "ipr", "shannon" | Criticality metric: IPR (default, recommended) or Shannon entropy (legacy) |
| `n_adaptive_sharpening` | bool | True | True/False | Apply $\sqrt{\ln N}$ sharpening to z-scores before softmax (IPR mode only) |

**Parameter Effects**:

| Parameter | Increase Effect | Decrease Effect |
|-----------|----------------|-----------------|
| `alpha` | More exploration when heated | Less exploration when heated |
| `beta` | More exploration when cooled | Less exploration when cooled |
| `alpha/beta` ratio | Larger CATS multiplier range | Smaller CATS multiplier range |
| `exploration_phase_end` | Longer pure exploration phase | Shorter pure exploration phase |
| `transition_phase_end` | Longer transition to full CATS | Faster transition to full CATS |
| `min_observations` | More conservative criticality | Earlier criticality activation |

**Example**:
```python
from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig

strategy_config = RouletteWheelConfig(
    mode="minimize",              # For docking scores
    alpha=0.15,                   # More aggressive exploration when heated
    beta=0.03,                    # Stronger exploitation when cooled
    exploration_phase_end=0.15,   # Start CATS earlier
    transition_phase_end=0.50,    # Reach full CATS sooner
    min_observations=3,           # Trust criticality with fewer samples
)
```

**Recommended Presets**:

| Scenario | alpha | beta | exploration_phase_end | transition_phase_end |
|----------|-------|------|----------------------|---------------------|
| Default (balanced) | 0.10 | 0.05 | 0.20 | 0.60 |
| Aggressive exploration | 0.20 | 0.05 | 0.30 | 0.70 |
| Fast convergence | 0.08 | 0.02 | 0.10 | 0.40 |
| Large libraries | 0.15 | 0.05 | 0.25 | 0.65 |

### 10.3 Warmup Strategy Parameters

#### 10.3.1 Balanced Warmup (Recommended)

**Class**: `BalancedWarmupConfig`

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `observations_per_reagent` | int | 5 | [3, 50] | Guaranteed observations per reagent ($K$) |
| `seed` | int | None | Any int or None | Random seed for reproducibility |
| `use_per_reagent_variance` | bool | True | True/False | Use per-reagent variance with James-Stein shrinkage |
| `shrinkage_strength` | float | 3.0 | > 0 | Shrinkage parameter ($\lambda$) for variance estimation |

**Shrinkage formula**:
$$\text{weight} = \frac{n}{n + \lambda}$$

With default $\lambda = 3.0$:
- 3 observations: weight = 0.5 (50% per-reagent, 50% global)
- 5 observations: weight = 0.625 (62.5% per-reagent)
- 10 observations: weight = 0.77 (77% per-reagent)

**Example**:
```python
from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig

warmup_config = BalancedWarmupConfig(
    observations_per_reagent=5,
    seed=42,                      # Reproducible
    use_per_reagent_variance=True,
    shrinkage_strength=3.0,
)
```

#### 10.3.2 Standard Warmup (Legacy)

**Class**: `StandardWarmupConfig`

No additional parameters. Uses `num_warmup_trials` from main config.

#### 10.3.3 Enhanced Warmup (Legacy)

**Class**: `EnhancedWarmupConfig`

No additional parameters. **Warning**: Creates imbalanced posteriors due to over-sampling of small components.

### 10.4 Alternative Selection Strategies

#### 10.4.1 Greedy Selection

**Class**: `GreedyConfig`

| Parameter | Type | Default | Options | Description |
|-----------|------|---------|---------|-------------|
| `mode` | str | "maximize" | "maximize", "minimize" | Optimization direction |

Pure exploitation - always selects highest/lowest posterior mean.

#### 10.4.2 UCB Selection

**Class**: `UCBConfig`

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `mode` | str | "maximize" | "maximize", "minimize" | Optimization direction |
| `c` | float | 2.0 | > 0 | Exploration parameter (higher = more exploration) |

Uses classical UCB formula: $\text{UCB}_i = \mu_i + c\sqrt{\frac{\ln(t)}{n_i}}$

#### 10.4.3 Epsilon-Greedy Selection

**Class**: `EpsilonGreedyConfig`

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `mode` | str | "maximize" | "maximize", "minimize" | Optimization direction |
| `epsilon` | float | 0.1 | [0, 1] | Initial exploration probability |
| `decay` | float | 0.995 | (0, 1] | Decay rate per iteration |

After $t$ iterations: $\epsilon_t = \epsilon_0 \times \text{decay}^t$

### 10.5 Parameter Tuning Guidelines

#### When to Adjust Parameters

| Symptom | Likely Cause | Parameter Adjustment |
|---------|--------------|---------------------|
| Missing top compounds | Under-exploration | Increase `alpha`, decrease `beta`, or increase `exploration_phase_end` |
| Slow convergence | Over-exploration | Decrease `alpha/beta` ratio, decrease `transition_phase_end` |
| Unstable criticality | Insufficient warmup | Increase `num_warmup_trials` or `min_observations` |
| Poor early performance | Inadequate warmup | Increase `observations_per_reagent` |
| All components critical too early | Low variance in scores | Increase `alpha`, or check if score range is appropriate |
| All components flexible too late | High variance in scores | Decrease `alpha`, or normalize scores |

#### Library Size Recommendations

| Library Size | Recommended `num_ts_iterations` | `num_warmup_trials` | `alpha` |
|--------------|--------------------------------|---------------------|---------|
| < 100K compounds | 1,000 - 5,000 | 3-5 | 0.08 - 0.10 |
| 100K - 1M compounds | 5,000 - 20,000 | 5 | 0.10 - 0.12 |
| 1M - 10M compounds | 20,000 - 50,000 | 5-7 | 0.12 - 0.15 |
| > 10M compounds | 50,000+ | 7-10 | 0.15 - 0.20 |

---

## 11. References

1. Thompson, W. R. (1933). On the likelihood that one unknown probability exceeds another in view of the evidence of two samples. *Biometrika*, 25(3-4), 285-294.

2. Zhao, H., Nittinger, E., & Tyrchan, C. (2024). Enhanced Thompson Sampling by Roulette Wheel Selection for Screening Ultra-Large Combinatorial Libraries. *bioRxiv* 2024.05.16.594622.

3. Russo, D. J., Van Roy, B., Kazerouni, A., Osband, I., & Wen, Z. (2018). A Tutorial on Thompson Sampling. *Foundations and Trends in Machine Learning*, 11(1), 1-96.

4. Shannon, C. E. (1948). A Mathematical Theory of Communication. *Bell System Technical Journal*, 27(3), 379-423.

5. James, W., & Stein, C. (1961). Estimation with Quadratic Loss. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 361-379.

6. Cover, T. M., & Thomas, J. A. (2006). *Elements of Information Theory* (2nd ed.). Wiley-Interscience.

---

*Document Version: 3.0*
*Last Updated: March 2026*
*Authors: TACTICS Development Team*
