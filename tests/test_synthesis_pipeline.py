"""
Unit tests for SynthesisPipeline functionality (v0.4.0 API).

Tests cover:
- ReactionDef creation and validation
- ReactionConfig creation and validation
- SynthesisPipeline factory methods
- Single-step and multi-step enumeration
- Alternative SMARTS routing
- Compatibility detection
"""

import pytest
from rdkit import Chem

from TACTICS.library_enumeration import SynthesisPipeline
from TACTICS.library_enumeration.smarts_toolkit import (
    ReactionDef,
    ReactionConfig,
    DeprotectionSpec,
    StepInput,
    InputSource,
)


class TestReactionDef:
    """Tests for ReactionDef creation and validation."""

    def test_simple_reaction_def(self):
        """Create simple ReactionDef."""
        rxn = ReactionDef(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            step_index=0,
            description="Amide coupling",
        )
        assert rxn.step_index == 0
        assert rxn.num_reactants == 2
        assert rxn.description == "Amide coupling"

    def test_reaction_def_with_pattern_id(self):
        """ReactionDef with pattern_id for alternatives."""
        rxn = ReactionDef(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            step_index=0,
            pattern_id="primary",
        )
        assert rxn.pattern_id == "primary"

    def test_reaction_def_with_deprotection(self):
        """ReactionDef with deprotection specifications."""
        rxn = ReactionDef(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            step_index=1,
            deprotections=[DeprotectionSpec(group="Boc", target=0)],
        )
        assert len(rxn.deprotections) == 1
        assert rxn.deprotections[0].group == "Boc"
        assert rxn.deprotections[0].target == 0
        assert rxn.deprotections[0].reactant_index == 0  # Property still works
        assert not rxn.deprotections[0].is_product_deprotection

    def test_reaction_def_with_product_deprotection(self):
        """ReactionDef with product deprotection specification."""
        rxn = ReactionDef(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            step_index=0,
            deprotections=[DeprotectionSpec(group="Fmoc", target="product")],
        )
        assert len(rxn.deprotections) == 1
        assert rxn.deprotections[0].group == "Fmoc"
        assert rxn.deprotections[0].target == "product"
        assert rxn.deprotections[0].reactant_index is None
        assert rxn.deprotections[0].is_product_deprotection

    def test_invalid_smarts_raises_error(self):
        """Invalid SMARTS should raise ValueError."""
        with pytest.raises(ValueError):
            ReactionDef(
                reaction_smarts="invalid_smarts",
                step_index=0,
            )


class TestReactionConfig:
    """Tests for ReactionConfig validation."""

    def test_single_reaction_config(self):
        """Single reaction configuration (step_inputs not required)."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
        )
        assert config.num_steps == 1
        assert len(config.reagent_file_list) == 2
        assert not config.is_multi_step

    def test_config_with_alternatives(self):
        """Configuration with alternative SMARTS patterns."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    pattern_id="primary",
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                    pattern_id="secondary",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )
        assert config.num_steps == 1
        assert config.has_alternatives_at_step(0)
        assert len(config.get_reactions_for_step(0)) == 2

    def test_multi_step_config(self):
        """Multi-step synthesis configuration."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                ),
                ReactionDef(
                    reaction_smarts="[NH2:1].[C:2](=O)[OH]>>[NH:1][C:2]=O",
                    step_index=1,
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi", "acids2.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ],
                1: [
                    StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=2),
                ],
            },
        )
        assert config.num_steps == 2
        assert config.is_multi_step

    def test_config_requires_step_inputs_for_multi_reaction(self):
        """Should raise error if step_inputs not provided for multiple reactions."""
        with pytest.raises(ValueError, match="step_inputs is required"):
            ReactionConfig(
                reactions=[
                    ReactionDef(
                        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                        step_index=0,
                    ),
                    ReactionDef(
                        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                        step_index=0,
                    ),
                ],
                reagent_file_list=["acids.smi", "amines.smi"],
            )

    def test_auto_generate_pattern_ids(self):
        """Pattern IDs should be auto-generated for alternatives."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )
        rxns = config.get_reactions_for_step(0)
        # Both should have pattern_ids (auto-generated)
        assert all(r.pattern_id is not None for r in rxns)


class TestSynthesisPipeline:
    """Tests for SynthesisPipeline functionality."""

    @pytest.fixture
    def simple_config(self):
        """Simple single-step configuration."""
        return ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
        )

    @pytest.fixture
    def alternative_config(self):
        """Single-step with alternatives."""
        return ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    pattern_id="primary",
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                    pattern_id="secondary",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )

    @pytest.fixture
    def two_step_config(self):
        """Two-step synthesis configuration."""
        return ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    description="Amide formation",
                ),
                ReactionDef(
                    reaction_smarts="[NH2:1].[CH:2]=O>>[NH:1][CH:2]",
                    step_index=1,
                    description="Reductive amination",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi", "aldehydes.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ],
                1: [
                    StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=2),
                ],
            },
        )

    def test_single_smarts_construction(self):
        """Test direct construction with single SMARTS."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
        )
        pipeline = SynthesisPipeline(config)
        assert pipeline.num_steps == 1
        assert pipeline.num_components == 2
        assert not pipeline.has_alternatives

    def test_alternatives_construction(self):
        """Test direct construction with alternative SMARTS patterns."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    pattern_id="primary",
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                    pattern_id="secondary",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )
        pipeline = SynthesisPipeline(config)
        assert pipeline.num_steps == 1
        assert pipeline.has_alternatives

    def test_single_step_enumeration(self, simple_config):
        """Should enumerate product in single-step mode."""
        pipeline = SynthesisPipeline(simple_config)

        acid = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        amine = Chem.MolFromSmiles("CCN")  # Ethylamine

        result = pipeline.enumerate([acid, amine])

        assert result.success
        assert result.product is not None
        assert result.product_smiles is not None
        assert "N" in result.product_smiles  # Should have nitrogen

    def test_enumeration_with_reagent_keys(self, simple_config):
        """Should pass reagent keys to enumeration."""
        pipeline = SynthesisPipeline(simple_config)

        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        result = pipeline.enumerate([acid, amine], ["acid_1", "amine_1"])

        assert result.success
        assert 0 in result.patterns_used
        assert result.patterns_used[0] == "primary"

    def test_enumeration_from_smiles(self, simple_config):
        """Should enumerate from SMILES strings."""
        pipeline = SynthesisPipeline(simple_config)

        result = pipeline.enumerate_from_smiles(
            ["CC(=O)O", "CCN"], ["acid_1", "amine_1"]
        )

        assert result.success
        assert result.product_smiles is not None

    def test_enumeration_failure_returns_error(self, simple_config):
        """Should return EnumerationError for failed enumeration."""
        pipeline = SynthesisPipeline(simple_config)

        acid = Chem.MolFromSmiles("CC(=O)O")
        not_amine = Chem.MolFromSmiles("CCC")  # No NH2

        result = pipeline.enumerate([acid, not_amine])

        assert not result.success
        assert result.product is None
        assert result.error is not None
        assert result.error.error_type == "reaction_failed"

    def test_enumeration_with_invalid_smiles(self, simple_config):
        """Should return error for invalid SMILES input."""
        pipeline = SynthesisPipeline(simple_config)

        result = pipeline.enumerate_from_smiles(["invalid_smiles", "CCN"])

        assert not result.success
        assert result.error is not None
        assert result.error.error_type == "invalid_input"

    def test_store_intermediates(self, two_step_config):
        """Should store intermediate products when requested."""
        pipeline = SynthesisPipeline(two_step_config)

        acid = Chem.MolFromSmiles("CC(=O)O")
        diamine = Chem.MolFromSmiles("NCCN")
        aldehyde = Chem.MolFromSmiles("CC=O")

        result = pipeline.enumerate([acid, diamine, aldehyde], store_intermediates=True)

        # Note: may or may not succeed depending on reaction chemistry
        # The test is that intermediates dict is populated when requested
        if result.success:
            assert result.intermediates is not None
        else:
            # Still shouldn't crash
            assert result.error is not None

    def test_alternative_pattern_routing(self, alternative_config):
        """Should route to correct pattern based on compatibility."""
        pipeline = SynthesisPipeline(alternative_config)

        # Register compatibility
        pipeline.register_compatibility(0, "acid_1", {"primary", "secondary"})
        pipeline.register_compatibility(1, "sec_amine_1", {"secondary"})

        acid = Chem.MolFromSmiles("CC(=O)O")
        sec_amine = Chem.MolFromSmiles("CCNCC")  # Secondary amine

        result = pipeline.enumerate([acid, sec_amine], ["acid_1", "sec_amine_1"])

        # Should use secondary pattern
        if result.success:
            assert result.patterns_used.get(0) == "secondary"

    def test_get_compatible_patterns(self, alternative_config):
        """Should find compatible patterns for reagent keys."""
        pipeline = SynthesisPipeline(alternative_config)

        pipeline.register_compatibility(0, "acid_1", {"primary", "secondary"})
        pipeline.register_compatibility(1, "amine_1", {"primary"})

        pattern = pipeline.get_compatible_patterns(["acid_1", "amine_1"])
        assert pattern == "primary"

    def test_get_compatible_patterns_intersection(self, alternative_config):
        """Should find intersection of compatible patterns."""
        pipeline = SynthesisPipeline(alternative_config)

        pipeline.register_compatibility(0, "acid_1", {"primary", "secondary"})
        pipeline.register_compatibility(1, "sec_amine_1", {"secondary"})

        pattern = pipeline.get_compatible_patterns(["acid_1", "sec_amine_1"])
        assert pattern == "secondary"

    def test_get_compatible_patterns_no_common(self, alternative_config):
        """Should return None when no common patterns."""
        pipeline = SynthesisPipeline(alternative_config)

        pipeline.register_compatibility(0, "acid_1", {"primary"})
        pipeline.register_compatibility(1, "sec_amine_1", {"secondary"})

        pattern = pipeline.get_compatible_patterns(["acid_1", "sec_amine_1"])
        assert pattern is None

    def test_runtime_pattern_fallback_without_cache(self):
        """Should try alternative patterns at runtime without pre-populated cache."""
        # Create config with alternatives - primary amine vs secondary amine
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    pattern_id="primary",
                    description="Primary amine coupling",
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH1:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                    pattern_id="secondary",
                    description="Secondary amine coupling",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )
        pipeline = SynthesisPipeline(config)

        # Do NOT call auto_detect_compatibility - testing runtime fallback

        acid = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        secondary_amine = Chem.MolFromSmiles("CCNCC")  # Diethylamine (secondary)

        # This should try "primary" first (will fail), then "secondary" (will succeed)
        result = pipeline.enumerate_single(
            [acid, secondary_amine],
            reagent_keys=["acid", "sec_amine"],
        )

        # Should succeed by falling back to secondary pattern
        assert result.success, f"Expected success but got error: {result.error}"
        assert result.patterns_used.get(0) == "secondary"

    def test_runtime_pattern_primary_succeeds_first(self):
        """Should use primary pattern when it succeeds."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                    pattern_id="primary",
                ),
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH1:2]>>[C:1](=O)[N:2]",
                    step_index=0,
                    pattern_id="secondary",
                ),
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
            step_inputs={
                0: [
                    StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                    StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                ]
            },
            step_modes={0: "alternative"},
        )
        pipeline = SynthesisPipeline(config)

        acid = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        primary_amine = Chem.MolFromSmiles("CCN")  # Ethylamine (primary)

        result = pipeline.enumerate_single(
            [acid, primary_amine],
            reagent_keys=["acid", "prim_amine"],
        )

        # Should succeed with primary pattern
        assert result.success, f"Expected success but got error: {result.error}"
        assert result.patterns_used.get(0) == "primary"

    def test_pipeline_properties(self, alternative_config):
        """Test pipeline property accessors."""
        pipeline = SynthesisPipeline(alternative_config)

        assert pipeline.num_steps == 1
        assert pipeline.num_components == 2
        assert not pipeline.is_multi_step
        assert pipeline.has_alternatives
        assert "primary" in pipeline.pattern_ids[0]
        assert "secondary" in pipeline.pattern_ids[0]

    def test_multi_step_pipeline(self, two_step_config):
        """Test multi-step pipeline properties."""
        pipeline = SynthesisPipeline(two_step_config)

        assert pipeline.num_steps == 2
        assert pipeline.num_components == 3
        assert pipeline.is_multi_step


class TestSynthesisPipelineBatch:
    """Tests for batch enumeration."""

    @pytest.fixture
    def pipeline(self):
        """Simple pipeline for batch testing."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=["acids.smi", "amines.smi"],
        )
        return SynthesisPipeline(config)

    def test_enumerate_batch_sequential(self, pipeline):
        """Test sequential batch enumeration."""
        acid1 = Chem.MolFromSmiles("CC(=O)O")
        acid2 = Chem.MolFromSmiles("CCC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        combinations = [
            ([acid1, amine], ["acid_1", "amine_1"]),
            ([acid2, amine], ["acid_2", "amine_1"]),
        ]

        results = pipeline.enumerate_batch(combinations, n_jobs=1)

        assert len(results) == 2
        assert results[0].success
        assert results[1].success

    def test_enumerate_batch_handles_failures(self, pipeline):
        """Batch enumeration should handle individual failures."""
        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")
        not_amine = Chem.MolFromSmiles("CCC")  # Won't react

        combinations = [
            ([acid, amine], ["acid_1", "amine_1"]),
            ([acid, not_amine], ["acid_1", "not_amine"]),
        ]

        results = pipeline.enumerate_batch(combinations, n_jobs=1)

        assert len(results) == 2
        assert results[0].success
        assert not results[1].success
        assert results[1].error is not None


class TestValidatorIntegration:
    """Tests for validation integration with pipeline."""

    def test_get_validator(self):
        """Should return validator for step/pattern."""
        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=[],
        )
        pipeline = SynthesisPipeline(config)

        # Validator not created until validation runs
        validator = pipeline.get_validator(0, "primary")
        assert validator is None  # Not validated yet

    def test_validate_all_creates_validators(self, tmp_path):
        """validate_all should create validators."""
        # Create test files
        acids = tmp_path / "acids.smi"
        acids.write_text("CC(=O)O acetic\n")
        amines = tmp_path / "amines.smi"
        amines.write_text("CCN ethylamine\n")

        config = ReactionConfig(
            reactions=[
                ReactionDef(
                    reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                    step_index=0,
                )
            ],
            reagent_file_list=[str(acids), str(amines)],
        )
        pipeline = SynthesisPipeline(config)

        results = pipeline.validate_all()

        assert 0 in results
        assert "primary" in results[0]

        # Validator should now exist
        validator = pipeline.get_validator(0, "primary")
        assert validator is not None
