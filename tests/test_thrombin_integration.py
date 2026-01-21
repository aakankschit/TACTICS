"""
Integration tests using bundled Thrombin dataset.

These tests verify that:
1. Bundled data files are accessible via importlib.resources
2. LookupEvaluator works with real product_scores.csv
3. SynthesisPipeline works with real reagent files
4. Full Thompson Sampling workflow (warmup → search) works with real data
5. Product naming matches scores file
"""

import pytest
import importlib.resources
import polars as pl

from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator
from TACTICS.thompson_sampling.strategies.greedy_selection import GreedySelection
from TACTICS.thompson_sampling.strategies.epsilon_greedy import EpsilonGreedySelection
from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef


# Amide coupling reaction SMARTS from tutorial
AMIDE_COUPLING_SMARTS = (
    "[#6:1](=[O:2])[OH]."
    "[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]"
    ">>[#6:1](=[O:2])[#7:3]"
)


class TestThrombinDataAccess:
    """Test bundled Thrombin data is accessible."""

    def test_bundled_files_exist(self):
        """Verify all bundled data files are accessible."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")

        # Check each file exists and is a file
        acids = data_files / "acids.smi"
        coupled_aa = data_files / "coupled_aa_sub.smi"
        amino_acids = data_files / "amino_acids_no_fmoc.smi"
        scores = data_files / "product_scores.csv"

        assert acids.is_file(), "acids.smi not found"
        assert coupled_aa.is_file(), "coupled_aa_sub.smi not found"
        assert amino_acids.is_file(), "amino_acids_no_fmoc.smi not found"
        assert scores.is_file(), "product_scores.csv not found"

    def test_acids_file_format(self):
        """Verify acids.smi has correct format."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")

        with open(acids_path, "r") as f:
            first_line = f.readline().strip()

        # Format: SMILES SPACE NAME
        parts = first_line.split()
        assert len(parts) >= 2, "Expected SMILES NAME format"
        assert "CA" in parts[-1], "Acid name should contain 'CA'"

    def test_coupled_aa_file_format(self):
        """Verify coupled_aa_sub.smi has correct format."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        aa_path = str(data_files / "coupled_aa_sub.smi")

        with open(aa_path, "r") as f:
            first_line = f.readline().strip()

        # Format: SMILES SPACE NAME (AA#_AA#)
        parts = first_line.split()
        assert len(parts) >= 2, "Expected SMILES NAME format"
        assert "AA" in parts[-1], "Name should contain 'AA'"
        assert "_" in parts[-1], "Name should have underscore (AA#_AA#)"

    def test_scores_file_format(self):
        """Verify product_scores.csv has correct format."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        scores_path = str(data_files / "product_scores.csv")

        df = pl.read_csv(scores_path, n_rows=10)

        assert "Product_Code" in df.columns, "Missing Product_Code column"
        assert "Scores" in df.columns, "Missing Scores column"

        # Check product code format (CA#_AA#_AA#)
        first_code = df["Product_Code"][0]
        assert first_code.startswith("CA"), "Product code should start with 'CA'"
        assert "_AA" in first_code, "Product code should contain '_AA'"

    def test_scores_count(self):
        """Verify product_scores.csv has expected row count."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        scores_path = str(data_files / "product_scores.csv")

        df = pl.read_csv(scores_path)

        # Should have ~500K products
        assert len(df) > 400000, f"Expected ~500K products, got {len(df)}"


class TestThrombinLookupEvaluator:
    """Test LookupEvaluator with real Thrombin scores."""

    @pytest.fixture
    def evaluator(self):
        """Create LookupEvaluator with Thrombin scores."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        scores_path = str(data_files / "product_scores.csv")

        return LookupEvaluator({"ref_filename": scores_path})

    def test_evaluator_loads_scores(self, evaluator):
        """Test that evaluator loads the scores file."""
        assert evaluator is not None
        assert evaluator.counter == 0

    def test_lookup_returns_score(self, evaluator):
        """Test that lookup returns a valid score for known product."""
        # Get a known product code from the file
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        scores_path = str(data_files / "product_scores.csv")
        df = pl.read_csv(scores_path, n_rows=1)
        known_product = df["Product_Code"][0]
        expected_score = df["Scores"][0]

        # LookupEvaluator.evaluate uses name to look up score
        score = evaluator.evaluate(known_product)

        assert score is not None
        assert abs(score - expected_score) < 0.001
        assert evaluator.counter == 1


class TestThrombinSynthesisPipeline:
    """Test SynthesisPipeline with real Thrombin reagent files."""

    @pytest.fixture
    def pipeline(self):
        """Create SynthesisPipeline with Thrombin data."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")
        amines_path = str(data_files / "coupled_aa_sub.smi")

        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts=AMIDE_COUPLING_SMARTS,
                    step_index=0,
                    description="Amide coupling",
                )
            ],
            reagent_file_list=[acids_path, amines_path],
        )
        return SynthesisPipeline(config)

    def test_pipeline_creation(self, pipeline):
        """Test that pipeline is created successfully."""
        assert pipeline is not None
        assert pipeline.num_components == 2
        assert pipeline.num_steps == 1

    def test_pipeline_enumeration(self, pipeline):
        """Test that pipeline can enumerate products."""
        from rdkit import Chem

        # Read first acid and amine
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")
        amines_path = str(data_files / "coupled_aa_sub.smi")

        with open(acids_path, "r") as f:
            acid_line = f.readline().strip()
        with open(amines_path, "r") as f:
            amine_line = f.readline().strip()

        acid_smiles = acid_line.split("\t")[0]
        amine_smiles = amine_line.split()[0]

        acid_mol = Chem.MolFromSmiles(acid_smiles)
        amine_mol = Chem.MolFromSmiles(amine_smiles)

        assert acid_mol is not None, f"Failed to parse acid SMILES: {acid_smiles}"
        assert amine_mol is not None, f"Failed to parse amine SMILES: {amine_smiles}"

        result = pipeline.enumerate([acid_mol, amine_mol])

        assert result.success, f"Enumeration failed: {result.error}"
        assert result.product is not None
        assert result.product_smiles is not None


class TestThrombinWorkflow:
    """Test full Thompson Sampling workflow with Thrombin data."""

    @pytest.fixture
    def sampler(self):
        """Create ThompsonSampler with Thrombin data."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")
        amines_path = str(data_files / "coupled_aa_sub.smi")
        scores_path = str(data_files / "product_scores.csv")

        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts=AMIDE_COUPLING_SMARTS,
                    step_index=0,
                    description="Amide coupling",
                )
            ],
            reagent_file_list=[acids_path, amines_path],
        )
        pipeline = SynthesisPipeline(config)

        # Use minimize mode for docking scores (lower is better)
        strategy = GreedySelection(mode="minimize")
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)
        sampler.read_reagents([acids_path, amines_path])

        evaluator = LookupEvaluator({"ref_filename": scores_path})
        sampler.set_evaluator(evaluator)

        yield sampler
        sampler.close()

    def test_sampler_setup(self, sampler):
        """Test that sampler is set up correctly."""
        assert sampler is not None
        assert sampler.evaluator is not None
        assert len(sampler.reagent_lists) == 2

        # Check reagent counts
        acids_count = len(sampler.reagent_lists[0])
        amines_count = len(sampler.reagent_lists[1])

        assert acids_count == 130, f"Expected 130 acids, got {acids_count}"
        assert amines_count == 3844, f"Expected 3844 amines, got {amines_count}"

    def test_single_evaluation(self, sampler):
        """Test evaluating a single compound."""
        smiles, name, score = sampler.evaluate([0, 0])

        assert smiles is not None
        assert name is not None
        assert isinstance(score, (int, float))
        assert sampler.evaluator.counter == 1

    def test_batch_evaluation(self, sampler):
        """Test evaluating a batch of compounds."""
        indices = [[0, 0], [0, 1], [1, 0], [1, 1]]
        results = sampler.evaluate_batch(indices)

        assert len(results) == 4
        for smiles, name, score in results:
            assert smiles is not None
            assert name is not None
            assert isinstance(score, (int, float))

    def test_warmup_phase(self, sampler):
        """Test warmup phase with small number of trials."""
        warmup_df = sampler.warm_up(num_warmup_trials=5)

        assert len(warmup_df) > 0
        assert "SMILES" in warmup_df.columns
        assert "Name" in warmup_df.columns
        assert "score" in warmup_df.columns

    def test_search_phase(self, sampler):
        """Test search phase with small number of cycles."""
        # Run brief warmup first
        sampler.warm_up(num_warmup_trials=3)

        # Run brief search
        search_df = sampler.search(num_cycles=10)

        assert len(search_df) > 0
        assert "SMILES" in search_df.columns
        assert "Name" in search_df.columns
        assert "score" in search_df.columns


class TestThrombinProductNaming:
    """Test that product naming matches scores file convention."""

    def test_product_name_format_matches_scores(self):
        """Verify product names match scores file format (CA#_AA#_AA#)."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")
        amines_path = str(data_files / "coupled_aa_sub.smi")
        scores_path = str(data_files / "product_scores.csv")

        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts=AMIDE_COUPLING_SMARTS,
                    step_index=0,
                    description="Amide coupling",
                )
            ],
            reagent_file_list=[acids_path, amines_path],
        )
        pipeline = SynthesisPipeline(config)

        strategy = GreedySelection(mode="minimize")
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)
        sampler.read_reagents([acids_path, amines_path])

        evaluator = LookupEvaluator({"ref_filename": scores_path})
        sampler.set_evaluator(evaluator)

        # Evaluate a compound and check naming
        smiles, name, score = sampler.evaluate([0, 0])

        # Name should match format CA#_AA#_AA#
        parts = name.split("_")
        assert len(parts) == 3, f"Expected 3 parts in name, got {len(parts)}: {name}"
        assert parts[0].startswith("CA"), f"First part should start with 'CA': {name}"
        assert parts[1].startswith("AA"), f"Second part should start with 'AA': {name}"
        assert parts[2].startswith("AA"), f"Third part should start with 'AA': {name}"

        # Score should be valid (lookup succeeded)
        assert score is not None and score != 0, f"Score lookup may have failed for {name}"

        sampler.close()


class TestThrombinEpsilonGreedy:
    """Test Epsilon-Greedy strategy with Thrombin data."""

    def test_epsilon_greedy_workflow(self):
        """Test epsilon-greedy strategy finds good compounds."""
        data_files = importlib.resources.files("TACTICS.data.thrombin")
        acids_path = str(data_files / "acids.smi")
        amines_path = str(data_files / "coupled_aa_sub.smi")
        scores_path = str(data_files / "product_scores.csv")

        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts=AMIDE_COUPLING_SMARTS,
                    step_index=0,
                    description="Amide coupling",
                )
            ],
            reagent_file_list=[acids_path, amines_path],
        )
        pipeline = SynthesisPipeline(config)

        # Epsilon-Greedy with minimize mode for docking scores
        strategy = EpsilonGreedySelection(mode="minimize", epsilon=0.2)
        sampler = ThompsonSampler(pipeline, selection_strategy=strategy)
        sampler.read_reagents([acids_path, amines_path])

        evaluator = LookupEvaluator({"ref_filename": scores_path})
        sampler.set_evaluator(evaluator)

        # Run brief warmup and search
        warmup_df = sampler.warm_up(num_warmup_trials=5)
        search_df = sampler.search(num_cycles=20)

        assert len(warmup_df) > 0
        assert len(search_df) > 0

        # Check that scores are reasonable (negative for docking scores)
        scores = search_df["score"].to_list()
        assert all(isinstance(s, (int, float)) for s in scores)

        sampler.close()
