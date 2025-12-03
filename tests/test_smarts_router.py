"""
Unit tests for SMARTSRouter functionality.
"""

import pytest
from rdkit import Chem
from TACTICS.library_enumeration.smarts_toolkit.smarts_router import (
    SMARTSPatternConfig,
    AlternativeSMARTSConfig,
    SMARTSRouter
)


class TestSMARTSPatternConfig:
    """Tests for SMARTSPatternConfig validation."""

    def test_valid_smarts(self):
        """Valid SMARTS should be accepted."""
        config = SMARTSPatternConfig(
            pattern_id="test",
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        )
        assert config.pattern_id == "test"
        assert config.reaction_smarts == "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"

    def test_valid_smarts_with_description(self):
        """Valid SMARTS with description should be accepted."""
        config = SMARTSPatternConfig(
            pattern_id="amide",
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            description="Primary amine amide coupling"
        )
        assert config.description == "Primary amine amide coupling"

    def test_invalid_smarts_raises(self):
        """Invalid SMARTS should raise ValidationError."""
        from pydantic_core import ValidationError
        with pytest.raises(ValidationError):
            SMARTSPatternConfig(
                pattern_id="test",
                reaction_smarts="not_valid_smarts"
            )


class TestAlternativeSMARTSConfig:
    """Tests for AlternativeSMARTSConfig."""

    def test_single_smarts_mode(self):
        """Single SMARTS mode should work."""
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        )
        assert not config.has_alternatives()
        patterns = config.get_all_patterns()
        assert len(patterns) == 1
        assert patterns[0].pattern_id == "primary"
        assert patterns[0].reaction_smarts == "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"

    def test_multiple_smarts_mode(self):
        """Multiple SMARTS patterns should be ordered correctly."""
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="alt_1",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                )
            ]
        )
        assert config.has_alternatives()
        patterns = config.get_all_patterns()
        assert len(patterns) == 2
        assert patterns[0].pattern_id == "primary"
        assert patterns[1].pattern_id == "alt_1"

    def test_multiple_alternatives(self):
        """Multiple alternative patterns should all be included."""
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="secondary_amine",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                ),
                SMARTSPatternConfig(
                    pattern_id="tertiary_amine",
                    reaction_smarts="[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
                )
            ]
        )
        patterns = config.get_all_patterns()
        assert len(patterns) == 3
        assert [p.pattern_id for p in patterns] == ["primary", "secondary_amine", "tertiary_amine"]

    def test_invalid_primary_smarts_raises(self):
        """Invalid primary SMARTS should raise ValidationError."""
        from pydantic_core import ValidationError
        with pytest.raises(ValidationError):
            AlternativeSMARTSConfig(
                primary_smarts="invalid_smarts"
            )


class TestSMARTSRouter:
    """Tests for SMARTSRouter functionality."""

    @pytest.fixture
    def simple_router(self):
        """Router with single SMARTS pattern."""
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        )
        return SMARTSRouter(config)

    @pytest.fixture
    def multi_router(self):
        """Router with multiple SMARTS patterns."""
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="secondary_amine",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                    description="For secondary amines"
                )
            ]
        )
        return SMARTSRouter(config)

    def test_initialization_single_pattern(self, simple_router):
        """Router should initialize correctly with single pattern."""
        assert simple_router.num_patterns == 1
        assert "primary" in simple_router.pattern_ids
        assert simple_router.pattern_ids == ["primary"]

    def test_initialization_multiple_patterns(self, multi_router):
        """Router should initialize correctly with multiple patterns."""
        assert multi_router.num_patterns == 2
        assert "primary" in multi_router.pattern_ids
        assert "secondary_amine" in multi_router.pattern_ids
        assert multi_router.pattern_ids == ["primary", "secondary_amine"]

    def test_register_compatibility_valid(self, multi_router):
        """Should register reagent compatibility correctly."""
        multi_router.register_reagent_compatibility(
            "reagent_A", {"primary", "secondary_amine"}
        )
        assert "reagent_A" in multi_router._compatibility_cache
        assert multi_router._compatibility_cache["reagent_A"] == {"primary", "secondary_amine"}

    def test_register_compatibility_single_pattern(self, multi_router):
        """Should register single pattern compatibility."""
        multi_router.register_reagent_compatibility(
            "reagent_B", {"primary"}
        )
        assert multi_router._compatibility_cache["reagent_B"] == {"primary"}

    def test_register_invalid_pattern_raises(self, multi_router):
        """Should raise error for unknown pattern IDs."""
        with pytest.raises(ValueError, match="Unknown pattern IDs"):
            multi_router.register_reagent_compatibility(
                "reagent_A", {"nonexistent_pattern"}
            )

    def test_register_mixed_valid_invalid_raises(self, multi_router):
        """Should raise error if any pattern ID is invalid."""
        with pytest.raises(ValueError, match="Unknown pattern IDs"):
            multi_router.register_reagent_compatibility(
                "reagent_A", {"primary", "nonexistent_pattern"}
            )

    def test_find_compatible_smarts_no_cache(self, multi_router):
        """Should default to primary when no compatibility registered."""
        result = multi_router.find_compatible_smarts(["reagent_A", "reagent_B"])
        assert result == "primary"

    def test_find_compatible_smarts_all_compatible(self, multi_router):
        """Should find common pattern when all reagents compatible."""
        multi_router.register_reagent_compatibility(
            "reagent_A", {"primary", "secondary_amine"}
        )
        multi_router.register_reagent_compatibility(
            "reagent_B", {"primary", "secondary_amine"}
        )

        result = multi_router.find_compatible_smarts(["reagent_A", "reagent_B"])
        # Should return first in priority order (primary)
        assert result == "primary"

    def test_find_compatible_smarts_intersection(self, multi_router):
        """Should find intersection of compatible patterns."""
        multi_router.register_reagent_compatibility(
            "reagent_A", {"primary", "secondary_amine"}
        )
        multi_router.register_reagent_compatibility(
            "reagent_B", {"secondary_amine"}
        )

        result = multi_router.find_compatible_smarts(["reagent_A", "reagent_B"])
        assert result == "secondary_amine"

    def test_find_compatible_smarts_no_common(self, multi_router):
        """Should return None when no common patterns."""
        multi_router.register_reagent_compatibility("reagent_A", {"primary"})
        multi_router.register_reagent_compatibility(
            "reagent_B", {"secondary_amine"}
        )

        result = multi_router.find_compatible_smarts(["reagent_A", "reagent_B"])
        assert result is None

    def test_find_compatible_smarts_unknown_reagent(self, multi_router):
        """Should assume primary for unknown reagents."""
        multi_router.register_reagent_compatibility(
            "reagent_A", {"primary", "secondary_amine"}
        )

        # reagent_B not registered - should assume primary
        result = multi_router.find_compatible_smarts(["reagent_A", "reagent_B"])
        assert result == "primary"

    def test_enumerate_success_simple(self, simple_router):
        """Should successfully enumerate product with simple router."""
        acid = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        amine = Chem.MolFromSmiles("CCN")      # Ethylamine

        product, pattern_id = simple_router.enumerate([acid, amine])

        assert product is not None
        assert pattern_id == "primary"
        # Product should be N-ethylacetamide (CCNC(C)=O or similar)
        product_smiles = Chem.MolToSmiles(product)
        assert "N" in product_smiles
        # Check for amide functional group (carbonyl)
        assert "=O" in product_smiles

    def test_enumerate_failure_returns_none(self, simple_router):
        """Should return None for failed enumeration."""
        # Incompatible reagents (no NH2)
        acid = Chem.MolFromSmiles("CC(=O)O")
        not_amine = Chem.MolFromSmiles("CCC")

        product, pattern_id = simple_router.enumerate([acid, not_amine])

        assert product is None
        assert pattern_id == "primary"  # Pattern was attempted

    def test_enumerate_with_explicit_pattern(self, multi_router):
        """Should use explicit pattern when provided."""
        acid = Chem.MolFromSmiles("CC(=O)O")
        sec_amine = Chem.MolFromSmiles("CCNCC")  # Diethylamine

        product, pattern_id = multi_router.enumerate(
            [acid, sec_amine],
            pattern_id="secondary_amine"
        )

        assert product is not None
        assert pattern_id == "secondary_amine"

    def test_enumerate_with_routing(self, multi_router):
        """Should route to correct SMARTS based on compatibility."""
        multi_router.register_reagent_compatibility(
            "acid_1", {"primary", "secondary_amine"}
        )
        multi_router.register_reagent_compatibility(
            "sec_amine_1", {"secondary_amine"}
        )

        acid = Chem.MolFromSmiles("CC(=O)O")
        sec_amine = Chem.MolFromSmiles("CCNCC")  # Diethylamine

        product, pattern_id = multi_router.enumerate(
            [acid, sec_amine],
            reagent_keys=["acid_1", "sec_amine_1"]
        )

        assert pattern_id == "secondary_amine"
        assert product is not None

    def test_enumerate_routing_failure(self, multi_router):
        """Should return None when routing fails (no common pattern)."""
        multi_router.register_reagent_compatibility("acid_1", {"primary"})
        multi_router.register_reagent_compatibility("amine_1", {"secondary_amine"})

        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        product, pattern_id = multi_router.enumerate(
            [acid, amine],
            reagent_keys=["acid_1", "amine_1"]
        )

        assert product is None
        assert pattern_id is None

    def test_get_reaction_primary(self, simple_router):
        """Should return primary reaction object."""
        rxn = simple_router.get_reaction("primary")
        assert rxn is not None
        assert rxn.GetNumReactantTemplates() == 2

    def test_get_reaction_invalid_pattern(self, simple_router):
        """Should return None for invalid pattern."""
        rxn = simple_router.get_reaction("nonexistent")
        assert rxn is None

    def test_get_compatibility_summary_empty(self, multi_router):
        """Should return zero counts when no reagents registered."""
        summary = multi_router.get_compatibility_summary()
        assert summary == {"primary": 0, "secondary_amine": 0}

    def test_get_compatibility_summary(self, multi_router):
        """Should correctly count compatible reagents."""
        multi_router.register_reagent_compatibility("r1", {"primary"})
        multi_router.register_reagent_compatibility("r2", {"primary", "secondary_amine"})
        multi_router.register_reagent_compatibility("r3", {"secondary_amine"})

        summary = multi_router.get_compatibility_summary()
        assert summary["primary"] == 2
        assert summary["secondary_amine"] == 2

    def test_enumerate_without_reagent_keys_defaults_primary(self, multi_router):
        """Should default to primary when no reagent_keys provided."""
        acid = Chem.MolFromSmiles("CC(=O)O")
        amine = Chem.MolFromSmiles("CCN")

        product, pattern_id = multi_router.enumerate([acid, amine])

        assert product is not None
        assert pattern_id == "primary"

    def test_complex_molecule_enumeration(self, simple_router):
        """Should handle more complex molecules."""
        # Benzoic acid + benzylamine
        acid = Chem.MolFromSmiles("c1ccccc1C(=O)O")
        amine = Chem.MolFromSmiles("NCc1ccccc1")

        product, pattern_id = simple_router.enumerate([acid, amine])

        assert product is not None
        assert pattern_id == "primary"
        product_smiles = Chem.MolToSmiles(product)
        # Should contain amide functional group (N and carbonyl)
        assert "N" in product_smiles
        assert ("=O" in product_smiles or "O=" in product_smiles)
        # Verify it's not just the starting materials
        assert product.GetNumAtoms() > acid.GetNumAtoms()
