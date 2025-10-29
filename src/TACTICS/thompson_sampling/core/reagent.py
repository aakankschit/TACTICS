"""
Unified Reagent class for Thompson Sampling.

This module provides a single Reagent class that works with all selection strategies.
"""

import numpy as np
from rdkit import Chem


class Reagent:
    """
    Unified reagent class for Thompson Sampling.

    Handles both warmup and search phases with Bayesian updating of
    posterior distributions.

    Attributes:
        reagent_name (str): Name/ID of the reagent
        smiles (str): SMILES string representation
        mol (Mol): RDKit molecule object
        mean (float): Current posterior mean
        std (float): Current posterior standard deviation
        n_samples (int): Number of times this reagent has been sampled
        known_var (float): Known variance (set during initialization)
        initial_scores (list): Scores collected during warmup
        current_phase (str): Either "warmup" or "search"
    """

    __slots__ = [
        "reagent_name",
        "smiles",
        "mol",
        "mean",
        "std",
        "n_samples",
        "known_var",
        "initial_scores",
        "current_phase"
    ]

    def __init__(self, reagent_name: str, smiles: str):
        """
        Initialize a reagent.

        Parameters:
            reagent_name: Unique identifier for this reagent
            smiles: SMILES string representation of the molecule
        """
        self.reagent_name = reagent_name
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(self.smiles)
        self.mean = 0.0
        self.std = 0.0
        self.n_samples = 0
        self.known_var = None
        self.initial_scores = []
        self.current_phase = "warmup"

    def add_score(self, score: float) -> None:
        """
        Add an observed score for this reagent.

        During warmup, scores are collected. During search, Bayesian
        updating is performed on the posterior distribution.

        Parameters:
            score: Observed score value
        """
        self.n_samples += 1

        if self.current_phase == "warmup":
            self.initial_scores.append(score)
        elif self.current_phase == "search":
            # Bayesian update
            current_var = self.std ** 2
            self.mean = self._update_mean(current_var, score)
            self.std = self._update_std(current_var)
        else:
            raise ValueError(f"Invalid phase: {self.current_phase}")

    def init_prior(self, prior_mean: float, prior_std: float) -> None:
        """
        Initialize the prior distribution from warmup statistics.

        This is called after warmup to set the prior mean and std,
        then replays all warmup scores as Bayesian updates.

        Parameters:
            prior_mean: Mean of the prior distribution (from warmup)
            prior_std: Standard deviation of the prior (from warmup)
        """
        if self.current_phase != "warmup":
            raise ValueError(f"Reagent {self.reagent_name} has already been initialized")

        if not self.initial_scores:
            raise ValueError(f"No warmup scores for reagent {self.reagent_name}")

        # Set prior parameters
        self.mean = prior_mean
        self.std = prior_std
        self.known_var = prior_std ** 2

        # Transition to search phase
        self.current_phase = "search"

        # Replay warmup scores as Bayesian updates
        warmup_scores = self.initial_scores.copy()
        self.initial_scores = []  # Clear to avoid double-counting
        for score in warmup_scores:
            self.add_score(score)

    def sample(self) -> float:
        """
        Sample a value from the current posterior distribution.

        Returns:
            float: Random sample from N(mean, std^2)
        """
        if self.current_phase != "search":
            raise ValueError(f"Must initialize prior before sampling")

        return np.random.normal(loc=self.mean, scale=self.std)

    def _update_mean(self, current_var: float, observed_value: float) -> float:
        """
        Bayesian update for the posterior mean.

        Parameters:
            current_var: Current variance
            observed_value: New observed score

        Returns:
            float: Updated posterior mean
        """
        numerator = current_var * observed_value + self.known_var * self.mean
        denominator = current_var + self.known_var
        return numerator / denominator

    def _update_std(self, current_var: float) -> float:
        """
        Bayesian update for the posterior standard deviation.

        Parameters:
            current_var: Current variance

        Returns:
            float: Updated posterior standard deviation
        """
        numerator = current_var * self.known_var
        denominator = current_var + self.known_var
        return np.sqrt(numerator / denominator)

    def __repr__(self) -> str:
        return f"Reagent('{self.reagent_name}', mean={self.mean:.3f}, std={self.std:.3f}, n={self.n_samples})"
