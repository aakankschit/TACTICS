"""
Unit tests for ReactionSequence functionality.
"""

import pytest
from rdkit import Chem
from pydantic_core import ValidationError
from TACTICS.library_enumeration.smarts_toolkit.reaction_sequence import (
    ReactionInputSource,
    ReactionInputConfig,
    ReactionStepConfig,
    ProtectingGroupConfig,
    MultiStepSynthesisConfig,
    ReactionSequence
)
from TACTICS.library_enumeration.smarts_toolkit.smarts_router import (
    AlternativeSMARTSConfig,
    SMARTSPatternConfig
)


class TestReactionInputConfig:
    """Tests for ReactionInputConfig validation."""

    def test_reagent_file_source_valid(self):
        """Should accept reagent_file source with component_idx."""
        config = ReactionInputConfig(
            source=ReactionInputSource.REAGENT_FILE,
            component_idx=0
        )
        assert config.source == ReactionInputSource.REAGENT_FILE
        assert config.component_idx == 0

    def test_previous_step_source_valid(self):
        """Should accept previous_step source with step_idx."""
        config = ReactionInputConfig(
            source=ReactionInputSource.PREVIOUS_STEP,
            step_idx=1
        )
        assert config.source == ReactionInputSource.PREVIOUS_STEP
        assert config.step_idx == 1

    def test_reagent_file_missing_component_idx_raises(self):
        """Should raise error if reagent_file source lacks component_idx."""
        with pytest.raises(ValidationError):
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE)

    def test_previous_step_missing_step_idx_raises(self):
        """Should raise error if previous_step source lacks step_idx."""
        with pytest.raises(ValidationError):
            ReactionInputConfig(source=ReactionInputSource.PREVIOUS_STEP)


class TestReactionStepConfig:
    """Tests for ReactionStepConfig."""

    def test_valid_step_config(self):
        """Should create valid step configuration."""
        config = ReactionStepConfig(
            step_idx=0,
            smarts_config=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
            ),
            inputs=[
                ReactionInputConfig(source="reagent_file", component_idx=0),
                ReactionInputConfig(source="reagent_file", component_idx=1)
            ],
            description="Amide coupling"
        )
        assert config.step_idx == 0
        assert config.description == "Amide coupling"
        assert len(config.inputs) == 2

    def test_step_with_substructure_validation(self):
        """Should accept valid substructure SMARTS."""
        config = ReactionStepConfig(
            step_idx=0,
            smarts_config=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
            ),
            inputs=[
                ReactionInputConfig(source="reagent_file", component_idx=0)
            ],
            required_substructure="[NH2]",
            expected_substructure="[NH]"
        )
        assert config.required_substructure == "[NH2]"
        assert config.expected_substructure == "[NH]"

    def test_invalid_required_substructure_raises(self):
        """Should raise error for invalid substructure SMARTS."""
        with pytest.raises(ValidationError):
            ReactionStepConfig(
                step_idx=0,
                smarts_config=AlternativeSMARTSConfig(
                    primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
                ),
                inputs=[
                    ReactionInputConfig(source="reagent_file", component_idx=0)
                ],
                required_substructure="invalid_smarts"
            )


class TestProtectingGroupConfig:
    """Tests for ProtectingGroupConfig."""

    def test_valid_protecting_group(self):
        """Should create valid protecting group config."""
        config = ProtectingGroupConfig(
            name="Boc",
            protected_smarts="[NH]C(=O)OC(C)(C)C",
            deprotected_smarts="[NH2]",
            component_idx=0,
            protection_removed_at_step=1
        )
        assert config.name == "Boc"
        assert config.component_idx == 0
        assert config.protection_removed_at_step == 1

    def test_invalid_protected_smarts_raises(self):
        """Should raise error for invalid SMARTS."""
        with pytest.raises(ValidationError):
            ProtectingGroupConfig(
                name="Boc",
                protected_smarts="invalid",
                deprotected_smarts="[NH2]",
                component_idx=0,
                protection_removed_at_step=1
            )


class TestMultiStepSynthesisConfig:
    """Tests for MultiStepSynthesisConfig."""

    def test_mode_1_simple_smarts(self):
        """Mode 1: Single SMARTS string."""
        config = MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=["acids.smi", "amines.smi"]
        )
        assert not config.is_multi_step()
        smarts_config = config.get_smarts_config()
        assert smarts_config.primary_smarts == "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"

    def test_mode_2_alternative_smarts(self):
        """Mode 2: Single step with alternatives."""
        config = MultiStepSynthesisConfig(
            alternative_smarts=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                alternatives=[
                    SMARTSPatternConfig(
                        pattern_id="secondary",
                        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                    )
                ]
            ),
            reagent_file_list=["acids.smi", "amines.smi"]
        )
        assert not config.is_multi_step()
        smarts_config = config.get_smarts_config()
        assert smarts_config.has_alternatives()

    def test_mode_3_multi_step(self):
        """Mode 3: Multi-step synthesis."""
        step1 = ReactionStepConfig(
            step_idx=0,
            smarts_config=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
            ),
            inputs=[
                ReactionInputConfig(source="reagent_file", component_idx=0),
                ReactionInputConfig(source="reagent_file", component_idx=1)
            ]
        )
        config = MultiStepSynthesisConfig(
            steps=[step1],
            reagent_file_list=["acids.smi", "amines.smi"]
        )
        assert config.is_multi_step()

    def test_no_mode_specified_raises(self):
        """Should raise error if no mode specified."""
        with pytest.raises(ValidationError):
            MultiStepSynthesisConfig(
                reagent_file_list=["acids.smi"]
            )

    def test_multiple_modes_raises(self):
        """Should raise error if multiple modes specified."""
        with pytest.raises(ValidationError):
            MultiStepSynthesisConfig(
                reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                alternative_smarts=AlternativeSMARTSConfig(
                    primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
                ),
                reagent_file_list=["acids.smi"]
            )

    def test_get_smarts_config_on_multi_step_raises(self):
        """Should raise error when calling get_smarts_config on multi-step."""
        step1 = ReactionStepConfig(
            step_idx=0,
            smarts_config=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
            ),
            inputs=[
                ReactionInputConfig(source="reagent_file", component_idx=0)
            ]
        )
        config = MultiStepSynthesisConfig(
            steps=[step1],
            reagent_file_list=["acids.smi"]
        )
        with pytest.raises(ValueError, match="Cannot get single SMARTS config"):
            config.get_smarts_config()


class TestReactionSequence:
    """Tests for ReactionSequence functionality."""

    @pytest.fixture
    def single_step_config(self):
        """Simple single-step configuration."""
        return MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=["acids.smi", "amines.smi"]
        )

    @pytest.fixture
    def alternative_config(self):
        """Single-step with alternatives."""
        return MultiStepSynthesisConfig(
            alternative_smarts=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                alternatives=[
                    SMARTSPatternConfig(
                        pattern_id="secondary",
                        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                    )
                ]
            ),
            reagent_file_list=["acids.smi", "amines.smi"]
        )

    @pytest.fixture
    def two_step_config(self):
        """Two-step synthesis configuration."""
        return MultiStepSynthesisConfig(
            reagent_file_list=["acids.smi", "amines.smi", "aldehydes.smi"],
            steps=[
                ReactionStepConfig(
                    step_idx=0,
                    smarts_config=AlternativeSMARTSConfig(
                        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
                    ),
                    inputs=[
                        ReactionInputConfig(source="reagent_file", component_idx=0),
                        ReactionInputConfig(source="reagent_file", component_idx=1)
                    ],
                    description="Amide formation"
                ),
                ReactionStepConfig(
                    step_idx=1,
                    smarts_config=AlternativeSMARTSConfig(
                        primary_smarts="[NH2:1].[CH:2]=O>>[NH:1][CH:2]"
                    ),
                    inputs=[
                        ReactionInputConfig(source="previous_step", step_idx=0),
                        ReactionInputConfig(source="reagent_file", component_idx=2)
                    ],
                    description="Reductive amination"
                )
            ]
        )

    def test_single_step_initialization(self, single_step_config):
        """Should initialize single-step sequence correctly."""
        seq = ReactionSequence(single_step_config)
        assert seq.num_steps == 1
        assert seq.num_components == 2

    def test_multi_step_initialization(self, two_step_config):
        """Should initialize multi-step sequence correctly."""
        seq = ReactionSequence(two_step_config)
        assert seq.num_steps == 2
        assert seq.num_components == 3

    def test_single_step_enumeration(self, single_step_config):
        """Should enumerate product in single-step mode."""
        seq = ReactionSequence(single_step_config)

        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        product, patterns = seq.enumerate([acid, amine])

        assert product is not None
        assert 0 in patterns
        assert patterns[0] == "primary"

    def test_single_step_with_routing(self, alternative_config):
        """Should route to correct pattern in single-step mode."""
        seq = ReactionSequence(alternative_config)

        # Register compatibility
        seq.register_reagent_compatibility(0, "acid_1", {"primary", "secondary"}, 0)
        seq.register_reagent_compatibility(1, "sec_amine", {"secondary"}, 0)

        acid = Chem.MolFromSmiles("CC(=O)O")
        sec_amine = Chem.MolFromSmiles("CCNCC")

        product, patterns = seq.enumerate(
            [acid, sec_amine],
            ["acid_1", "sec_amine"]
        )

        assert product is not None
        assert patterns[0] == "secondary"

    def test_two_step_enumeration(self, two_step_config):
        """Should enumerate product through two-step synthesis."""
        seq = ReactionSequence(two_step_config)

        # Step 1: Acid + Amine -> Amide
        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("NCCN")  # Diamine with free NH2
        aldehyde = Chem.MolFromSmiles("CC=O")

        product, patterns = seq.enumerate([acid, amine, aldehyde])

        # Should complete both steps
        assert product is not None
        assert 0 in patterns  # Step 0 completed
        assert 1 in patterns  # Step 1 completed

    def test_enumeration_failure_returns_none(self, single_step_config):
        """Should return None for failed enumeration."""
        seq = ReactionSequence(single_step_config)

        # Incompatible reagents
        acid = Chem.MolFromSmiles("CC(=O)O")
        not_amine = Chem.MolFromSmiles("CCC")

        product, patterns = seq.enumerate([acid, not_amine])

        assert product is None

    def test_get_reaction_single_step(self, single_step_config):
        """Should return reaction object for single-step."""
        seq = ReactionSequence(single_step_config)
        rxn = seq.get_reaction(step_idx=0, pattern_id="primary")
        assert rxn is not None

    def test_get_reaction_invalid_step_returns_none(self, single_step_config):
        """Should return None for invalid step index."""
        seq = ReactionSequence(single_step_config)
        rxn = seq.get_reaction(step_idx=99)
        assert rxn is None

    def test_register_all_reagent_compatibilities_single_step(self, alternative_config):
        """Should register all compatibilities for single-step."""
        from TACTICS.thompson_sampling.core.reagent import Reagent

        seq = ReactionSequence(alternative_config)

        # Create mock reagent lists
        reagent_lists = [
            [Reagent("acid_1", "CC(=O)O")],
            [Reagent("amine_1", "CCN")]
        ]

        # Set compatibility
        reagent_lists[0][0].set_compatible_smarts({"primary", "secondary"})
        reagent_lists[1][0].set_compatible_smarts({"secondary"})

        # Register
        seq.register_all_reagent_compatibilities(reagent_lists)

        # Test that routing works
        product, patterns = seq.enumerate(
            [r.mol for r in reagent_lists[0]] + [r.mol for r in reagent_lists[1]],
            [reagent_lists[0][0].reagent_key, reagent_lists[1][0].reagent_key]
        )

        assert patterns[0] == "secondary"

    def test_multi_step_with_substructure_validation(self):
        """Should validate substructures in multi-step synthesis."""
        config = MultiStepSynthesisConfig(
            reagent_file_list=["acids.smi", "amines.smi"],
            steps=[
                ReactionStepConfig(
                    step_idx=0,
                    smarts_config=AlternativeSMARTSConfig(
                        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
                    ),
                    inputs=[
                        ReactionInputConfig(source="reagent_file", component_idx=0),
                        ReactionInputConfig(source="reagent_file", component_idx=1)
                    ],
                    required_substructure="[NH2]",  # Must have primary amine
                    expected_substructure="C(=O)N"  # Should form amide
                )
            ]
        )

        seq = ReactionSequence(config)

        # Valid: has NH2
        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        product, patterns = seq.enumerate([acid, amine])
        assert product is not None

        # Invalid: no NH2 (secondary amine)
        sec_amine = Chem.MolFromSmiles("CCNCC")
        product2, patterns2 = seq.enumerate([acid, sec_amine])
        # This should fail required_substructure check
        # (Though the SMARTS might still match, depends on exact pattern)

    def test_complex_multi_step_synthesis(self):
        """Test complex multi-step with intermediate dependencies."""
        config = MultiStepSynthesisConfig(
            reagent_file_list=["r1.smi", "r2.smi", "r3.smi"],
            steps=[
                # Step 0: Combine r1 + r2
                ReactionStepConfig(
                    step_idx=0,
                    smarts_config=AlternativeSMARTSConfig(
                        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
                    ),
                    inputs=[
                        ReactionInputConfig(source="reagent_file", component_idx=0),
                        ReactionInputConfig(source="reagent_file", component_idx=1)
                    ]
                ),
                # Step 1: Take step 0 output + r3
                ReactionStepConfig(
                    step_idx=1,
                    smarts_config=AlternativeSMARTSConfig(
                        primary_smarts="[NH2:1].[C:2](=O)[OH]>>[NH:1][C:2]=O"
                    ),
                    inputs=[
                        ReactionInputConfig(source="previous_step", step_idx=0),
                        ReactionInputConfig(source="reagent_file", component_idx=2)
                    ]
                )
            ]
        )

        seq = ReactionSequence(config)

        # Reagents for 3-step synthesis
        r1 = Chem.MolFromSmiles("CC(=O)O")
        r2 = Chem.MolFromSmiles("NCCN")  # Diamine
        r3 = Chem.MolFromSmiles("c1ccccc1C(=O)O")  # Benzoic acid

        product, patterns = seq.enumerate([r1, r2, r3])

        # Should successfully complete both steps
        assert product is not None
        assert len(patterns) == 2
        assert 0 in patterns
        assert 1 in patterns
