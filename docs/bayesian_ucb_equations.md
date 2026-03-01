# Bayesian Upper Confidence Bound Selection in TACTICS

## A Scientific Treatment with Mathematical Derivations

---

## Abstract

This document provides rigorous mathematical foundations for the Bayesian Upper Confidence Bound (Bayes-UCB) selection strategy as implemented in TACTICS. Bayes-UCB offers a principled alternative to Roulette Wheel Selection by using Student-t quantiles for proper Bayesian treatment of uncertainty. We derive the theoretical underpinnings, demonstrate the integration with Component-Aware Thompson Sampling (CATS), and compare Bayes-UCB systematically with classical UCB and RWS approaches. The key advantages of Bayes-UCB include deterministic selection given posteriors, explicit control over exploration via percentiles, and proper handling of small-sample uncertainty through heavy-tailed distributions.

---

## Table of Contents

1. [Theoretical Background](#1-theoretical-background)
2. [Derivation of Bayes-UCB](#2-derivation-of-bayes-ucb)
3. [Student-t Distribution in Bayesian Inference](#3-student-t-distribution-in-bayesian-inference)
4. [Integration with CATS](#4-integration-with-cats)
5. [Percentile-Based Thermal Cycling](#5-percentile-based-thermal-cycling)
6. [Complete Algorithm Specification](#6-complete-algorithm-specification)
7. [Regret Analysis](#7-regret-analysis)
8. [Comparison with Alternative Strategies](#8-comparison-with-alternative-strategies)
9. [Practical Considerations](#9-practical-considerations)
10. [References](#10-references)

---

## 1. Theoretical Background

### 1.1 The Multi-Armed Bandit Framework

In combinatorial library screening, each reagent selection problem can be formulated as a multi-armed bandit:

**Definition 1.1 (Stochastic Multi-Armed Bandit)**: A stochastic MAB consists of:
- $K$ arms (reagents), indexed $i \in \{1, \ldots, K\}$
- Unknown reward distributions $\nu_i$ with means $\mu_i^*$
- At each round $t$, the agent selects arm $I_t$ and receives reward $X_t \sim \nu_{I_t}$

**Objective**: Minimize cumulative regret:
$$R_T = T \cdot \mu^* - \sum_{t=1}^{T} \mu_{I_t}^*$$

where $\mu^* = \max_i \mu_i^*$.

### 1.2 From Frequentist to Bayesian UCB

#### Classical UCB1 (Auer et al., 2002)

**Definition 1.2 (UCB1 Index)**: The UCB1 index for arm $i$ at time $t$ is:
$$\text{UCB}_i(t) = \bar{X}_i(t) + c\sqrt{\frac{\ln(t)}{N_i(t)}}$$

where:
- $\bar{X}_i(t) = \frac{1}{N_i(t)}\sum_{s \leq t: I_s = i} X_s$ is the sample mean
- $N_i(t)$ is the number of times arm $i$ was pulled
- $c$ is an exploration constant (typically $c = \sqrt{2}$ for $[0,1]$ rewards)

**Derivation of Exploration Term**: The $\sqrt{\ln(t)/N_i}$ term comes from Hoeffding's inequality and the union bound over time.

**Theorem 1.1 (UCB1 Regret Bound)**: For rewards in $[0, 1]$, UCB1 achieves:
$$R_T \leq 8 \sum_{i: \mu_i < \mu^*} \frac{\ln T}{\Delta_i} + (1 + \frac{\pi^2}{3}) \sum_{i=1}^{K} \Delta_i$$

where $\Delta_i = \mu^* - \mu_i$ is the suboptimality gap.

#### Limitations of Classical UCB

1. **Ignores prior information**: Treats all arms as initially equivalent
2. **Bounded reward assumption**: Requires $[0, 1]$ rewards for standard analysis
3. **No posterior uncertainty**: Uses point estimates rather than distributions
4. **Fixed exploration formula**: Cannot incorporate domain knowledge

### 1.3 Bayesian Perspective

**Bayesian Motivation**: Instead of concentration inequalities, use the posterior distribution directly.

**Key Insight**: The upper confidence bound can be interpreted as a quantile of the posterior:
$$\text{UCB}_i = \text{Quantile}_{1-\alpha}(\theta_i | \mathcal{D})$$

where $\theta_i$ is the unknown mean of arm $i$ and $\mathcal{D}$ is the observed data.

---

## 2. Derivation of Bayes-UCB

### 2.1 Normal-Normal Model with Unknown Variance

Consider the model:
- **Prior on mean**: $\theta \sim \mathcal{N}(\mu_0, \sigma_0^2)$
- **Likelihood**: $X | \theta, \sigma^2 \sim \mathcal{N}(\theta, \sigma^2)$
- **Prior on variance**: $\sigma^2$ unknown, estimated from data

After observing $n$ samples $\{X_1, \ldots, X_n\}$ with sample mean $\bar{X}$ and sample variance $s^2$:

### 2.2 Posterior Distribution

**Theorem 2.1 (Posterior for Unknown Variance)**: Under the normal-normal model with unknown variance, the posterior predictive distribution for a new observation is Student-t:
$$\frac{\bar{X} - \theta}{s/\sqrt{n}} \sim t_{n-1}$$

where $t_{n-1}$ denotes the Student-t distribution with $n-1$ degrees of freedom.

**Proof**:

Let $\bar{X} = \frac{1}{n}\sum_{i=1}^n X_i$ and $s^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2$.

Under the normal model:
1. $\bar{X} \sim \mathcal{N}(\theta, \sigma^2/n)$
2. $(n-1)s^2/\sigma^2 \sim \chi^2_{n-1}$
3. $\bar{X}$ and $s^2$ are independent

Therefore:
$$T = \frac{\bar{X} - \theta}{s/\sqrt{n}} = \frac{(\bar{X} - \theta)/(\sigma/\sqrt{n})}{s/\sigma}$$

The numerator is $\mathcal{N}(0, 1)$ and the denominator involves $\sqrt{\chi^2_{n-1}/(n-1)}$, giving $T \sim t_{n-1}$. $\square$

### 2.3 Bayes-UCB Index Definition

**Definition 2.1 (Bayes-UCB Index)**: For reagent $i$ with $n_i$ observations, posterior mean $\mu_i$, and posterior standard deviation $\sigma_i$, the Bayes-UCB index at percentile $p$ is:

$$\text{UCB}_i(p) = \mu_i + \frac{\sigma_i \cdot t_{\nu_i}(p)}{\sqrt{n_i}}$$

where:
- $\nu_i = n_i - 1$ is the degrees of freedom
- $t_{\nu}(p)$ is the $p$-th quantile of the Student-t distribution with $\nu$ degrees of freedom

**For minimization mode (Lower Confidence Bound)**:
$$\text{LCB}_i(p) = \mu_i - \frac{\sigma_i \cdot t_{\nu_i}(p)}{\sqrt{n_i}}$$

#### 2.3.1 Note on the $\sqrt{n_i}$ Scaling Factor

A subtle but important detail: the formula includes division by $\sqrt{n_i}$ even though $\sigma_i$ is already the posterior standard deviation (which shrinks as observations accumulate).

**Why this additional scaling?**

In the standard Bayesian posterior, $\sigma_i$ already incorporates the effect of sample size:
$$\sigma_i^2 = \frac{\sigma_0^2}{n_i + 1} \approx \frac{\sigma_0^2}{n_i}$$

However, in TACTICS, the posterior standard deviation $\sigma_i$ is computed via Bayesian updates that may not perfectly track $1/\sqrt{n}$ behavior (especially during warmup or with heterogeneous observation variance). The additional $1/\sqrt{n_i}$ factor serves as:

1. **Conservative shrinkage**: Ensures exploration bonus decreases appropriately with sample size
2. **Finite-sample correction**: Accounts for potential overconfidence in small-sample posteriors
3. **Empirical effectiveness**: Matches the classical UCB form that has proven regret bounds

This design choice prioritizes **empirical robustness** over strict Bayesian purity. In practice, it prevents over-exploitation of reagents that happen to have artificially low posterior variance.

### 2.4 Properties of the Bayes-UCB Index

**Property 2.1 (Monotonicity in $n$)**: For fixed $\mu$, $\sigma$, and $p > 0.5$:
$$\frac{\partial \text{UCB}}{\partial n} < 0$$

The UCB decreases as more observations are collected.

**Property 2.2 (Limiting Behavior)**:
$$\lim_{n \to \infty} t_{\nu}(p) = z_p$$

where $z_p$ is the standard normal quantile. Thus Bayes-UCB approaches classical UCB as $n \to \infty$.

**Property 2.3 (Heavy Tails for Small $n$)**:
$$t_{\nu}(p) > z_p \quad \text{for all } \nu < \infty, \, p > 0.5$$

The Student-t has heavier tails, providing more exploration bonus for under-sampled arms.

### 2.5 Under-Explored Reagent Handling

**Problem**: For $n_i < 2$, the degrees of freedom $\nu_i \leq 0$, making the t-distribution undefined.

**Solution**: Use a conservative exploration bonus:
$$\text{UCB}_i = \mu_i + 3 \cdot \max(\sigma_i, 10^{-6})$$

**Justification**: The factor of 3 corresponds approximately to the 99.9th percentile of a standard normal, encouraging exploration of under-sampled reagents while maintaining numerical stability.

---

## 3. Student-t Distribution in Bayesian Inference

### 3.1 Definition and Properties

**Definition 3.1 (Student-t Distribution)**: A random variable $T$ has a Student-t distribution with $\nu$ degrees of freedom if:
$$f_T(t; \nu) = \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\,\Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{t^2}{\nu}\right)^{-\frac{\nu+1}{2}}$$

**Properties**:
- Symmetric about zero
- $\mathbb{E}[T] = 0$ for $\nu > 1$
- $\text{Var}(T) = \frac{\nu}{\nu - 2}$ for $\nu > 2$
- $T \xrightarrow{d} \mathcal{N}(0, 1)$ as $\nu \to \infty$

### 3.2 Quantile Function

**Definition 3.2 (t-Quantile)**: The $p$-th quantile of $t_\nu$ is the value $t_\nu(p)$ such that:
$$P(T \leq t_\nu(p)) = p$$

**Numerical values** for $p = 0.90$:

| $\nu$ | $t_\nu(0.90)$ | Ratio to $z_{0.90}$ |
|-------|---------------|---------------------|
| 1 | 3.078 | 2.40 |
| 2 | 1.886 | 1.47 |
| 5 | 1.476 | 1.15 |
| 10 | 1.372 | 1.07 |
| 30 | 1.310 | 1.02 |
| $\infty$ | 1.282 | 1.00 |

**Observation**: Small degrees of freedom ($\nu < 5$) significantly increase the exploration bonus.

### 3.3 Why Student-t Matters for Exploration

**Theorem 3.1 (Exploration Bonus Ratio)**: For fixed percentile $p > 0.5$, the ratio of Student-t to normal quantile satisfies:
$$\frac{t_\nu(p)}{z_p} \geq 1$$

with equality only as $\nu \to \infty$.

**Implication**: Under-explored arms (small $n$, hence small $\nu$) receive larger exploration bonuses than the normal approximation would suggest.

**Derivation of Asymptotic Expansion**: For large $\nu$:
$$t_\nu(p) \approx z_p + \frac{z_p^3 + z_p}{4\nu} + O(\nu^{-2})$$

The second term is always positive for $z_p > 0$, confirming the heavier tails.

---

## 4. Integration with CATS

> **📖 Reference**: This section summarizes CATS integration for Bayes-UCB. For complete derivations of criticality, progressive weight, temperature multipliers, and the relationship between thermal cycling and CATS, see **[thompson_sampling_equations.md, Sections 5.0-5.6](thompson_sampling_equations.md#5-component-aware-thompson-sampling-cats)**.

### 4.1 Component Criticality in Bayes-UCB Context

CATS (Component-Aware Thompson Sampling) extends Bayes-UCB by modulating the percentile level based on component criticality.

**Recall**: Component criticality $\kappa_c$ measures posterior concentration using z-score softmax with IPR (see [thompson_sampling_equations.md, Section 5.1](thompson_sampling_equations.md#51-component-criticality-via-z-score-softmax-with-ipr) for full derivation):

$$z_{c,i} = \frac{\mu_{c,i} - \bar{\mu}_c}{\sigma_{\mu,c}} \cdot \lambda_c \cdot \text{sharpening}$$

$$p_{c,i} = \text{softmax}(z_{c,i}), \qquad \text{IPR}_c = \sum_i p_{c,i}^2, \qquad \kappa_c = 1 - \frac{1/\text{IPR}_c}{n_c}$$

where $\lambda_c$ is the SNR dampening factor and sharpening = $\sqrt{\ln n_c}$ (IPR mode with N-adaptive sharpening).

### 4.2 Percentile as Exploration Control

In RWS, temperature controls exploration. In Bayes-UCB, **percentile** serves the analogous role:

| Parameter | Effect on Exploration |
|-----------|----------------------|
| $p \to 1$ | Wide confidence bounds → More exploration |
| $p \to 0.5$ | Tight confidence bounds → More exploitation |

**Mathematical Relationship**:
$$\lim_{p \to 1} t_\nu(p) = +\infty \quad \text{(infinite exploration bonus)}$$
$$t_\nu(0.5) = 0 \quad \text{(no exploration bonus)}$$

### 4.3 CATS Multiplier for Percentiles

**Definition 4.1 (CATS Percentile Multiplier)**: Analogous to temperature multipliers:
$$m_{\max} = \frac{p_{\text{high}}}{p_{\text{low}}}$$
$$m_{\min} = \frac{p_{\text{low}}}{p_{\text{high}}}$$

**Relative Neutral-Point Mapping** (see [thompson_sampling_equations.md, Section 5.2.2](thompson_sampling_equations.md#522-relative-neutral-point-multiplier)):

Let $\bar{\kappa}$ be the mean criticality across components:

$$m_c = \begin{cases}
m_{\max} + \frac{\kappa_c}{\bar{\kappa}} \cdot (1 - m_{\max}) & \text{if } \kappa_c \leq \bar{\kappa} \\[6pt]
1 + \frac{\kappa_c - \bar{\kappa}}{1 - \bar{\kappa}} \cdot (m_{\min} - 1) & \text{if } \kappa_c > \bar{\kappa}
\end{cases}$$

**Behavior**:
- Flexible component ($\kappa_c < \bar{\kappa}$): $m_c > 1$ → Higher percentile → More exploration
- Mean component ($\kappa_c = \bar{\kappa}$): $m_c = 1$ → No adjustment (neutral)
- Critical component ($\kappa_c > \bar{\kappa}$): $m_c < 1$ → Lower percentile → More exploitation

### 4.4 Progressive Blending

**Definition 4.2 (Effective Multiplier)**:
$$m_c^{\text{eff}} = 1 + w(\gamma) \cdot (m_c - 1)$$

where $w(\gamma)$ is the three-phase progressive weight:
$$w(\gamma) = \begin{cases}
0 & \gamma < 0.20 \\
\frac{\gamma - 0.20}{0.40} & 0.20 \leq \gamma < 0.60 \\
1 & \gamma \geq 0.60
\end{cases}$$

**Derivation of Effective Percentile**:
$$p_c^{\text{eff}} = \min\left(p_c^{\text{base}} \cdot m_c^{\text{eff}}, 1\right)$$

The clipping to 1 ensures valid percentile values.

---

## 5. Percentile-Based Thermal Cycling

### 5.1 Base Percentile Assignment

**Definition 5.1 (Thermal Cycling Percentile)**:
$$p_c^{\text{base}} = \begin{cases}
p_{\text{high}} & \text{if } c = c_{\text{hot}} \\
p_{\text{low}} & \text{otherwise}
\end{cases}$$

where:
- $p_{\text{high}} = 0.90$ (default): Heated component gets wider bounds
- $p_{\text{low}} = 0.60$ (default): Cooled components get tighter bounds

### 5.2 Component Rotation

After each selection cycle, the heated component is selected using criticality-weighted probabilities (see [thompson_sampling_equations.md, Section 5.4](thompson_sampling_equations.md#54-criticality-weighted-component-rotation)):

$$f_c = \max(1 - \kappa_c, 0.1), \qquad P(\text{heat } c) = \frac{f_c}{\sum_{c'} f_{c'}}, \qquad c_{\text{hot}} \sim \text{Categorical}(P)$$

Flexible components receive the exploration bonus more frequently than critical ones.

### 5.3 Complete Effective Percentile Formula

Combining all elements:

**Theorem 5.1 (Effective Percentile)**: The effective percentile for component $c$ at iteration $t$ is:
$$p_c^{\text{eff}}(t) = \min\left(p_c^{\text{base}} \cdot \left[1 + w\left(\frac{t}{T}\right) \cdot (m_c - 1)\right], 1\right)$$

where:
- $p_c^{\text{base}}$ comes from thermal cycling
- $w(\gamma)$ is the progressive weight
- $m_c$ is the CATS multiplier based on criticality

### 5.4 Worked Example

**Setup**:
- $p_{\text{high}} = 0.90$, $p_{\text{low}} = 0.60$
- Iteration $t = 4000$, $T = 10000$ ($\gamma = 0.4$)
- Component 1 (acids): $\kappa_1 = 0.8$ (critical), currently heated
- Component 2 (amines): $\kappa_2 = 0.2$ (flexible), cooled

**Derived parameters**:
- $m_{\max} = 0.90/0.60 = 1.5$
- $m_{\min} = 0.60/0.90 \approx 0.667$

**Progressive weight**: $w(0.4) = (0.4 - 0.2)/0.4 = 0.5$

**For acids** (heated, critical):
- $m_1 = 0.667 + (1.5 - 0.667)(1 - 0.8) = 0.667 + 0.167 = 0.833$
- $m_1^{\text{eff}} = 1 + 0.5 \cdot (0.833 - 1) = 0.917$
- $p_1^{\text{eff}} = 0.90 \times 0.917 = 0.825$

**For amines** (cooled, flexible):
- $m_2 = 0.667 + (1.5 - 0.667)(1 - 0.2) = 0.667 + 0.667 = 1.333$
- $m_2^{\text{eff}} = 1 + 0.5 \cdot (1.333 - 1) = 1.167$
- $p_2^{\text{eff}} = 0.60 \times 1.167 = 0.70$

**Result**:
- Acids: $p = 0.825$ (reduced from 0.90 due to criticality)
- Amines: $p = 0.70$ (increased from 0.60 due to flexibility)

---

## 6. Complete Algorithm Specification

### 6.1 Bayes-UCB Selection Algorithm

```
Algorithm: BAYES-UCB WITH CATS
─────────────────────────────────────────────────────────────────────
Input: reagent_list, component_idx, current_cycle, total_cycles
Output: selected_reagent_index

1. BASE PERCENTILE (Thermal Cycling)
   if component_idx == current_heated_component:
       p_base ← p_high
   else:
       p_base ← p_low

2. COMPONENT CRITICALITY (z-score softmax with IPR)
   active ← [r for r in reagent_list if r.n_samples > 0]
   p25_obs ← percentile([r.n_samples for r in active], 25)
   if p25_obs < N_min:
       κ ← 0.5  // neutral
   else:
       means ← [r.mean for r in active]
       if mode == "minimize": means ← -means
       z ← (means - mean(means)) / std(means)
       z ← z * SNR_dampening(active)     // suppress noise-dominated
       z ← z * sqrt(log(len(active)))    // N-adaptive sharpening (IPR)
       probs ← softmax(z)
       IPR ← sum(probs²)
       κ ← 1 - (1/IPR) / len(active)

3. PROGRESSIVE WEIGHT
   γ ← current_cycle / total_cycles
   if γ < 0.20:      w ← 0
   elif γ < 0.60:    w ← (γ - 0.20) / 0.40
   else:             w ← 1

4. CATS MULTIPLIER (relative neutral-point)
   κ̄ ← mean_criticality  // cached from weighted rotation
   if κ ≤ κ̄:
       m ← m_max + (κ / κ̄) * (1 - m_max)
   else:
       m ← 1 + ((κ - κ̄) / (1 - κ̄)) * (m_min - 1)
   m_eff ← 1 + w * (m - 1)

5. EFFECTIVE PERCENTILE
   p_eff ← clip(p_base * m_eff, 0, 1)

6. COMPUTE UCB INDICES
   for each reagent i:
       μ_i, σ_i, n_i ← reagent_list[i].mean, .std, .n_samples
       if n_i < 2:
           UCB_i ← μ_i + 3 * max(σ_i, 1e-6)  // conservative bonus
       else:
           ν_i ← n_i - 1
           t_quantile ← StudentT.ppf(p_eff, ν_i)
           UCB_i ← μ_i + σ_i * t_quantile / sqrt(n_i)

7. APPLY DISALLOW MASK
   if mode == "maximize":
       UCB[disallowed] ← -∞
   else:
       UCB[disallowed] ← +∞

8. SELECT
   if mode == "maximize":
       return argmax(UCB)
   else:
       return argmin(UCB)
─────────────────────────────────────────────────────────────────────
```

### 6.2 Pseudocode for Key Functions

**Criticality Calculation** (z-score softmax with IPR):
```python
def calculate_criticality(reagent_list, mode, min_observations=5,
                          criticality_metric="ipr", n_adaptive_sharpening=True):
    """
    Calculate component criticality using z-score softmax with
    SNR dampening and IPR (or Shannon entropy for legacy).

    Returns criticality κ ∈ [0, 1]:
        κ ≈ 0: Flexible (uniform posterior) → Explore
        κ ≈ 1: Critical (peaked posterior) → Exploit
    """
    active = [r for r in reagent_list if r.n_samples > 0]
    if len(active) < 2:
        return 0.5

    # p25 observation gating
    obs = [r.n_samples for r in active]
    if np.percentile(obs, 25) < min_observations:
        return 0.5

    means = np.array([r.mean for r in active])
    if mode == "minimize":
        means = -means

    # Z-score normalization (scale-invariant)
    mean_std = np.std(means)
    if mean_std < 1e-10:
        return 0.0
    z_scores = (means - np.mean(means)) / mean_std

    # SNR dampening
    se_sq = np.array([r.std**2 / max(r.n_samples, 1) for r in active])
    noise_std = np.sqrt(np.mean(se_sq))
    snr = mean_std / max(noise_std, 1e-10)
    z_scores *= 1.0 - np.exp(-max(snr - 1.0, 0.0))

    # N-adaptive sharpening (IPR mode only)
    N = len(active)
    if criticality_metric == "ipr" and n_adaptive_sharpening and N > 2:
        z_scores *= max(1.0, np.sqrt(np.log(N)))

    # Softmax with numerical stability
    probs = np.exp(z_scores - z_scores.max())
    probs /= probs.sum()

    if criticality_metric == "ipr":
        ipr = np.sum(probs ** 2)
        return 1.0 - (1.0 / ipr) / N
    else:  # shannon
        entropy = -np.sum(probs * np.log(probs + 1e-10))
        return 1.0 - entropy / np.log(N)
```

**UCB Index Computation**:
```python
def compute_ucb_indices(reagent_list, percentile, mode="maximize"):
    """
    Compute Bayes-UCB indices using Student-t quantiles.

    UCB_i = μ_i + σ_i * t_{ν_i}(p) / √n_i
    """
    n = len(reagent_list)
    ucb = np.zeros(n)

    for i, r in enumerate(reagent_list):
        if r.n_samples < 2:
            # Conservative bonus for under-explored
            bonus = 3.0 * max(r.std, 1e-6)
            ucb[i] = r.mean + bonus if mode == "maximize" else r.mean - bonus
        else:
            df = r.n_samples - 1
            t_quantile = scipy.stats.t.ppf(percentile, df)
            bonus = r.std * t_quantile / np.sqrt(r.n_samples)
            ucb[i] = r.mean + bonus if mode == "maximize" else r.mean - bonus

    return ucb
```

---

## 7. Regret Analysis

### 7.1 Bayes-UCB Regret Bounds

**Theorem 7.1 (Kaufmann et al., 2012)**: Under appropriate conditions on the prior, Bayes-UCB achieves asymptotically optimal regret:
$$R_T \sim \sum_{i: \Delta_i > 0} \frac{\Delta_i}{\text{KL}(\mu_i, \mu^*)} \ln T$$

where $\text{KL}(\mu_i, \mu^*)$ is the Kullback-Leibler divergence between the reward distributions.

### 7.2 Comparison of Regret Bounds

| Algorithm | Regret Bound | Optimality |
|-----------|--------------|------------|
| UCB1 | $O\left(\sum_i \frac{\ln T}{\Delta_i}\right)$ | Order-optimal |
| KL-UCB | $\sum_i \frac{\Delta_i}{\text{KL}(\mu_i, \mu^*)} \ln T + o(\ln T)$ | Asymptotically optimal |
| Bayes-UCB | Same as KL-UCB | Asymptotically optimal |
| Thompson Sampling | Same as KL-UCB | Asymptotically optimal |

**Key Result**: Bayes-UCB matches the theoretical optimality of Thompson Sampling while providing deterministic selection.

### 7.3 Finite-Time Analysis

**Theorem 7.2**: For Gaussian rewards with known variance $\sigma^2$, Bayes-UCB with quantile $1 - 1/t$ achieves:
$$\mathbb{E}[R_T] \leq O\left(\sqrt{KT \ln T}\right)$$

This matches the minimax lower bound up to logarithmic factors.

---

## 8. Comparison with Alternative Strategies

### 8.1 Bayes-UCB vs Roulette Wheel Selection

| Aspect | Bayes-UCB | RWS |
|--------|-----------|-----|
| **Selection mechanism** | Argmax on UCB indices | Probabilistic (roulette) |
| **Exploration control** | Percentile level | Temperature |
| **Stochasticity** | Deterministic | Stochastic |
| **Theoretical guarantees** | Regret bounds | Empirical |
| **CATS integration** | Percentile multiplier | Temperature multiplier |
| **Computational cost** | Higher (t-quantiles) | Lower (exponentials) |

**Mathematical Relationship**:

RWS selection probability:
$$P_{\text{RWS}}(i) = \frac{\exp(s_i / T)}{\sum_j \exp(s_j / T)}$$

Bayes-UCB implicit probability (via argmax):
$$P_{\text{UCB}}(i) = \mathbb{1}\left[i = \arg\max_j \text{UCB}_j\right]$$

### 8.2 Bayes-UCB vs Standard UCB

| Aspect | Bayes-UCB | Standard UCB |
|--------|-----------|--------------|
| **Uncertainty model** | Posterior distribution | Concentration inequality |
| **Exploration bonus** | $\sigma_i t_\nu(p) / \sqrt{n_i}$ | $c\sqrt{\ln(t)/n_i}$ |
| **Small sample handling** | Heavy-tailed (Student-t) | May over-exploit |
| **Variance treatment** | Estimated from data | Assumed bounded |
| **Prior incorporation** | Natural | Difficult |

**When Bayes-UCB is Preferred**:
1. Variance differs across arms
2. Prior information is available
3. Small sample sizes are common
4. Deterministic selection is desired

### 8.3 Bayes-UCB vs Thompson Sampling

| Aspect | Bayes-UCB | Thompson Sampling |
|--------|-----------|-------------------|
| **Core idea** | Optimism (upper bound) | Probability matching |
| **Selection** | Argmax on UCB | Argmax on samples |
| **Exploration source** | Explicit (percentile) | Implicit (posterior) |
| **Reproducibility** | Deterministic | Stochastic |
| **Computational cost** | Similar | Similar |
| **Theoretical guarantees** | Both asymptotically optimal | Both asymptotically optimal |

**Theorem 8.1 (Equivalence in Limit)**: As the percentile $p \to 1$ and sample size $n \to \infty$, Bayes-UCB and Thompson Sampling become equivalent in their selection probabilities.

---

## 9. Practical Considerations

### 9.1 Parameter Selection Guidelines

| Parameter | Default | Range | Guidance |
|-----------|---------|-------|----------|
| $p_{\text{high}}$ | 0.90 | [0.7, 0.99] | Higher = more exploration when heated |
| $p_{\text{low}}$ | 0.60 | [0.5, 0.85] | Higher = more exploration when cooled |
| $N_{\min}$ | 5 | [3, 10] | Minimum observations for criticality |
| $\gamma_1$ | 0.20 | [0.10, 0.30] | End of pure exploration |
| $\gamma_2$ | 0.60 | [0.40, 0.80] | End of transition |

### 9.2 Computational Efficiency

**Student-t quantile computation**:
- For unique degrees of freedom, cache ppf values
- Group reagents by $n_i$ to minimize ppf calls

```python
# Efficient batch computation
unique_dfs = np.unique(n_samples - 1)
t_quantiles = {}
for df in unique_dfs:
    t_quantiles[df] = scipy.stats.t.ppf(percentile, df)

# Apply cached values
for i, r in enumerate(reagent_list):
    df = r.n_samples - 1
    ucb[i] = r.mean + r.std * t_quantiles.get(df, 3.0) / np.sqrt(r.n_samples)
```

### 9.3 Numerical Stability

**Issue**: Very small or very large posteriors can cause numerical issues.

**Solutions**:
1. Clamp standard deviation: $\sigma_i \geq 10^{-6}$
2. Clamp percentile: $p \in [0.01, 0.99]$
3. Use log-space for extreme values

### 9.4 Batch Selection

> **⚠️ CRITICAL WARNING: Bayes-UCB Batch Selection is DETERMINISTIC**
>
> Unlike Roulette Wheel Selection (which samples probabilistically), Bayes-UCB uses `argmax` and will return **the same reagent every time** given identical posteriors and CATS state.
>
> **Implications**:
> - Calling `select_reagent()` $B$ times returns the **same reagent** $B$ times
> - For diverse batch selection, you **MUST** use `DisallowTracker` to exclude already-selected reagents
> - Alternatively, consider using RWS for batch mode if diversity without tracking is desired

For batch selection (selecting $B$ reagents simultaneously):

**Correct approach with DisallowTracker**:
```python
selected = []
for _ in range(batch_size):
    # DisallowTracker updates mask after each selection
    idx = strategy.select_reagent(reagent_list, disallow_mask=tracker.get_mask())
    selected.append(idx)
    tracker.update(idx)  # Prevents re-selection
return selected
```

**Simple approach (only for ranking, not diverse selection)**:
```python
return np.argsort(ucb_indices)[-B:]  # Top B for maximize - BUT may have duplicates in meaning
```

**When to use which strategy for batches**:

| Scenario | Recommended Strategy | Reason |
|----------|---------------------|--------|
| Diverse batch needed | RWS or Bayes-UCB + DisallowTracker | Determinism requires explicit exclusion |
| Ranking top reagents | Bayes-UCB without DisallowTracker | Determinism is actually desirable |
| Parallel evaluation | RWS with batch sampling | Natural diversity from stochasticity |

---

## 10. Tunable Parameters Reference

This section provides a complete reference for all tunable parameters specific to the Bayes-UCB selection strategy.

### 10.1 Bayes-UCB Configuration

**Class**: `BayesUCBConfig`

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `mode` | str | "maximize" | "maximize", "minimize" | Optimization direction |
| `initial_p_high` | float | 0.90 | [0.5, 0.999] | Base percentile for **heated** component (wider bounds) |
| `initial_p_low` | float | 0.60 | [0.5, 0.999] | Base percentile for **cooled** components (tighter bounds) |
| `exploration_phase_end` | float | 0.20 | (0, 1] | Fraction of iterations before CATS starts ($\gamma_1$) |
| `transition_phase_end` | float | 0.60 | (0, 1] | Fraction of iterations when CATS is fully active ($\gamma_2$) |
| `min_observations` | int | 5 | > 0 | Minimum observations per reagent before trusting criticality |
| `criticality_metric` | str | "ipr" | "ipr", "shannon" | Criticality metric: IPR (default, recommended) or Shannon entropy (legacy) |
| `n_adaptive_sharpening` | bool | True | True/False | Apply $\sqrt{\ln N}$ sharpening to z-scores before softmax (IPR mode only) |

### 10.2 Percentile vs Temperature Analogy

In Bayes-UCB, **percentile** plays the role that **temperature** plays in RWS:

| RWS Parameter | Bayes-UCB Equivalent | Effect |
|---------------|---------------------|--------|
| `alpha` (high temp) | `initial_p_high` | More exploration when heated |
| `beta` (low temp) | `initial_p_low` | Less exploration when cooled |
| `alpha/beta` ratio | `p_high/p_low` ratio | CATS multiplier range |

**Key Difference**:
- RWS temperature affects the **softmax distribution** (stochastic selection)
- Bayes-UCB percentile affects the **UCB index** (deterministic argmax selection)

### 10.3 Parameter Effects

| Parameter | Increase Effect | Decrease Effect |
|-----------|----------------|-----------------|
| `initial_p_high` | Wider confidence bounds when heated → more exploration | Tighter bounds → less exploration |
| `initial_p_low` | More exploration even when cooled | Stronger exploitation when cooled |
| `p_high/p_low` ratio | Larger CATS adjustment range | Smaller CATS adjustment range |
| `exploration_phase_end` | Longer pure exploration (no CATS) | Earlier CATS activation |
| `transition_phase_end` | Slower transition to full CATS | Faster full CATS effect |
| `min_observations` | More conservative criticality estimates | Earlier (potentially noisier) criticality |

### 10.4 Recommended Configurations

**Default (Balanced)**:
```python
from TACTICS.thompson_sampling.strategies.config import BayesUCBConfig

config = BayesUCBConfig(
    mode="minimize",               # For docking scores
    initial_p_high=0.90,           # 90th percentile when heated
    initial_p_low=0.60,            # 60th percentile when cooled
    exploration_phase_end=0.20,    # CATS starts at 20%
    transition_phase_end=0.60,     # Full CATS at 60%
    min_observations=5,
)
```

**Aggressive Exploration**:
```python
config = BayesUCBConfig(
    mode="minimize",
    initial_p_high=0.95,           # Very wide bounds when heated
    initial_p_low=0.55,            # Still some exploration when cooled
    exploration_phase_end=0.30,    # Longer exploration phase
    transition_phase_end=0.70,
    min_observations=3,
)
```

**Fast Convergence**:
```python
config = BayesUCBConfig(
    mode="minimize",
    initial_p_high=0.85,           # Moderate exploration
    initial_p_low=0.55,            # Exploitation-focused
    exploration_phase_end=0.10,    # Quick CATS start
    transition_phase_end=0.40,     # Fast transition
    min_observations=5,
)
```

### 10.5 CATS Multiplier Range Computation

The CATS multiplier range is derived from the percentile ratio:

$$m_{\max} = \frac{p_{\text{high}}}{p_{\text{low}}}$$
$$m_{\min} = \frac{p_{\text{low}}}{p_{\text{high}}}$$

| p_high | p_low | m_max | m_min | Range |
|--------|-------|-------|-------|-------|
| 0.90 | 0.60 | 1.50 | 0.67 | 0.83 |
| 0.95 | 0.55 | 1.73 | 0.58 | 1.15 |
| 0.85 | 0.65 | 1.31 | 0.76 | 0.55 |
| 0.99 | 0.50 | 1.98 | 0.51 | 1.47 |

**Interpretation**: Larger range = CATS has more ability to adjust exploration based on criticality.

### 10.6 Comparison with RWS Parameters

If you're familiar with RWS and want equivalent Bayes-UCB settings:

| RWS Setting | Equivalent Bayes-UCB |
|-------------|---------------------|
| `alpha=0.10, beta=0.05` | `p_high=0.90, p_low=0.60` (default) |
| `alpha=0.15, beta=0.05` | `p_high=0.95, p_low=0.55` (more exploration) |
| `alpha=0.08, beta=0.04` | `p_high=0.85, p_low=0.65` (less exploration) |

### 10.7 When to Choose Bayes-UCB vs RWS

| Scenario | Recommended Strategy | Reason |
|----------|---------------------|--------|
| Deterministic results needed | Bayes-UCB | Argmax is reproducible |
| Parallel batch evaluation | RWS | Natural diversity from stochasticity |
| Small sample sizes | Bayes-UCB | Student-t handles uncertainty better |
| Fast screening | RWS | Slightly faster computation |
| Theoretical guarantees needed | Bayes-UCB | Has regret bounds |
| Large batches without DisallowTracker | RWS | Bayes-UCB would pick same reagent repeatedly |

---

## 11. References

1. Kaufmann, E., Cappé, O., & Garivier, A. (2012). On Bayesian Upper Confidence Bounds for Bandit Problems. *AISTATS 2012*.

2. Auer, P., Cesa-Bianchi, N., & Fischer, P. (2002). Finite-time Analysis of the Multiarmed Bandit Problem. *Machine Learning*, 47(2-3), 235-256.

3. Russo, D. J., Van Roy, B., Kazerouni, A., Osband, I., & Wen, Z. (2018). A Tutorial on Thompson Sampling. *Foundations and Trends in Machine Learning*, 11(1), 1-96.

4. Agrawal, S., & Goyal, N. (2012). Analysis of Thompson Sampling for the Multi-armed Bandit Problem. *COLT 2012*.

5. Chapelle, O., & Li, L. (2011). An Empirical Evaluation of Thompson Sampling. *NeurIPS 2011*.

6. Zhao, H., Nittinger, E., & Tyrchan, C. (2024). Enhanced Thompson Sampling by Roulette Wheel Selection for Screening Ultra-Large Combinatorial Libraries. *bioRxiv* 2024.05.16.594622.

7. Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

---

*Document Version: 3.0*
*Last Updated: March 2026*
*Authors: TACTICS Development Team*
