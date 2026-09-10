import os
import warnings
from abc import ABC, abstractmethod
import numpy as np
import polars as pl

# OpenEye toolkit modules are loaded lazily by _ensure_openeye() rather than
# at module import time. They stay None until an OpenEye-backed evaluator is
# actually constructed.
oechem = None
oeomega = None
oeshape = None
oedocking = None


def _ensure_openeye():
    """Import the OpenEye toolkits on first use and cache them as module globals.

    Deferred rather than imported at module level for two reasons:

    1. Segfault avoidance. OpenEye's shared libraries load libexpat.so.1 and
       initialize its global C-level XML parser state. tqdm pulls in
       prompt_toolkit, whose progress-bar formatter calls
       xml.dom.minidom.parseString() at *module* level, re-entering libexpat
       through Python's pyexpat extension. When OpenEye initialized expat
       first, that second entry conflicts with the existing global state and
       segfaults the interpreter (exit 139) during `import TACTICS`, before
       any user code runs. Deferring the OpenEye load until an evaluator is
       constructed lets Python own and initialize expat first.

    2. Import cost. Most TACTICS runs use Lookup, FP, MW or DB evaluators and
       never need the OpenEye toolkits at all.

    Safe to call repeatedly; the import happens once.

    Raises:
        ImportError: if the OpenEye toolkits are not installed.
    """
    global oechem, oeomega, oeshape, oedocking

    if oechem is not None:
        return

    try:
        from openeye import oechem as _oechem
        from openeye import oeomega as _oeomega
        from openeye import oeshape as _oeshape
        from openeye import oedocking as _oedocking
    except ImportError as exc:
        raise ImportError(
            "The OpenEye toolkits are required for ROCSEvaluator and "
            "FredEvaluator but are not installed in this environment. "
            "Install them with `pip install openeye-toolkits` (a valid "
            "OpenEye license is required)."
        ) from exc

    oechem, oeomega, oeshape, oedocking = _oechem, _oeomega, _oeshape, _oedocking


from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sqlitedict import SqliteDict

# Shared Morgan (ECFP4-equivalent) generator: radius 2, 2048 bits.
_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

class Evaluator(ABC):
    """Base class for scoring functions.

    An evaluator turns one product into one number. The sampler calls
    :meth:`evaluate` once per product it decides to test and feeds the score
    into the reagent posteriors.

    Most evaluators are constructed by :func:`~TACTICS.thompson_sampling.factories.create_evaluator`
    from their paired Pydantic config (for example
    :class:`~TACTICS.thompson_sampling.core.evaluator_config.LookupEvaluatorConfig`
    builds a :class:`LookupEvaluator`), which is also how parallel workers rebuild
    them. Constructing one directly is fine for single-process use.

    Subclasses implement :meth:`evaluate` and the :attr:`counter` property.
    A score of ``NaN`` means "could not score"; the sampler skips it.
    """

    @abstractmethod
    def evaluate(self, mol):
        """Score one product.

        Args:
            mol: An RDKit ``Mol`` for structure-based evaluators, or the product
                *name* (``str``) for :class:`LookupEvaluator` and
                :class:`DBEvaluator`, which key on the product code.

        Returns:
            float: The score. Higher is better in ``mode="maximize"``; lower is
            better in ``mode="minimize"`` (docking).
        """

    @property
    @abstractmethod
    def counter(self):
        """Number of :meth:`evaluate` calls so far."""


class MWEvaluator(Evaluator):
    """Score = molecular weight. A smoke-test evaluator; it takes no arguments.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.MWEvaluatorConfig`.
    """

    def __init__(self):
        self.num_evaluations = 0

    @property
    def counter(self):
        return self.num_evaluations

    def evaluate(self, mol):
        self.num_evaluations += 1
        return Descriptors.MolWt(mol)


class FPEvaluator(Evaluator):
    """Score = Morgan-fingerprint Tanimoto similarity to a query molecule.

    Fingerprints are radius 2, 2048 bits (ECFP4-equivalent). Fast, needs no
    3D, no licence.

    Args:
        input_dict: ``{"query_smiles": str}`` -- the reference molecule.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.FPEvaluatorConfig`.

    Raises:
        ValueError: if ``query_smiles`` does not parse.
    """

    def __init__(self, input_dict):
        self.ref_smiles = input_dict["query_smiles"]
        ref_mol = Chem.MolFromSmiles(self.ref_smiles)
        if ref_mol is None:
            raise ValueError(f"Could not parse query_smiles: {self.ref_smiles!r}")
        self.ref_fp = _MORGAN_GEN.GetFingerprint(ref_mol)
        self.num_evaluations = 0

    @property
    def counter(self):
        return self.num_evaluations

    def evaluate(self, rd_mol_in):
        self.num_evaluations += 1
        rd_mol_fp = _MORGAN_GEN.GetFingerprint(rd_mol_in)
        return DataStructs.TanimotoSimilarity(self.ref_fp, rd_mol_fp)


class ROCSEvaluator(Evaluator):
    """Score = ROCS shape + colour Tanimoto combo to a 3D query (OpenEye).

    Conformers are generated with Omega on the fly (``max_confs``, default 50;
    change with :meth:`set_max_confs`). Slow: use ``processes > 1``.
    Requires the ``openeye`` extra and a licence.

    Args:
        input_dict: ``{"query_molfile": str}`` -- a 3D query file readable by
            ``oechem.oemolistream`` (SDF, MOL2, OEB).

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.ROCSEvaluatorConfig`.
    """

    def __init__(self, input_dict):
        _ensure_openeye()
        ref_filename = input_dict['query_molfile']
        ref_fs = oechem.oemolistream(ref_filename)
        self.ref_mol = oechem.OEMol()
        oechem.OEReadMolecule(ref_fs, self.ref_mol)
        self.max_confs = 50
        self.score_cache = {}
        self.num_evaluations = 0

    @property
    def counter(self):
        return self.num_evaluations

    def set_max_confs(self, max_confs):
        """Set the maximum number of conformers generated by Omega
        :param max_confs:
        """
        self.max_confs = max_confs

    def evaluate(self, rd_mol_in):
        """Generate conformers with Omega and evaluate the ROCS overlay of conformers to a reference molecule
        :param rd_mol_in: Input RDKit molecule
        :return: ROCS Tanimoto Combo score, returns -1 if conformer generation fails
        """
        self.num_evaluations += 1
        smi = Chem.MolToSmiles(rd_mol_in)
        # Look up to see if we already processed this molecule
        arc_tc = self.score_cache.get(smi)
        if arc_tc is not None:
            tc = arc_tc
        else:
            fit_mol = oechem.OEMol()
            oechem.OEParseSmiles(fit_mol, smi)
            ret_code = generate_confs(fit_mol, self.max_confs)
            if ret_code:
                tc = self.overlay(fit_mol)
            else:
                tc = -1.0
            self.score_cache[smi] = tc
        return tc

    def overlay(self, fit_mol):
        """Use ROCS to overlay two molecules
        :param fit_mol: OEMolecule
        :return: Combo Tanimoto for the overlay
        """
        prep = oeshape.OEOverlapPrep()
        prep.Prep(self.ref_mol)
        overlay = oeshape.OEMultiRefOverlay()
        overlay.SetupRef(self.ref_mol)
        prep.Prep(fit_mol)
        score = oeshape.OEBestOverlayScore()
        overlay.BestOverlay(score, fit_mol, oeshape.OEHighestTanimoto())
        return score.GetTanimotoCombo()


class LookupEvaluator(Evaluator):
    """Score = a value looked up by product code in a precomputed table.

    Used for benchmarking against exhaustive scores and for any workflow
    where scores already exist. Keyed on the product *name*
    (``<reagent1>_<reagent2>_...``), so the sampler skips product synthesis
    entirely when this evaluator is active.

    Args:
        input_dict: ``{"ref_filename": str, "compound_col": str = "Product_Code",
            "score_col": str = "Scores", "default_score": float | None = None}``.
            ``ref_filename`` may be ``.csv`` or ``.parquet``. ``default_score``
            is returned for product codes absent from the table; leave it
            ``None`` (→ ``NaN``, skipped) unless absence has a meaning, e.g.
            ``0.0`` for DEL read counts where an unlisted product is a
            non-binder. A JSON string of the same dict is also accepted.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.LookupEvaluatorConfig`.
    """

    def __init__(self, input_dictionary):
        self.num_evaluations = 0

        # Handle both dictionary and JSON string inputs
        if isinstance(input_dictionary, str):
            import json
            input_dictionary = json.loads(input_dictionary)

        ref_filename = input_dictionary['ref_filename']
        compound_col = input_dictionary.get('compound_col', 'Product_Code')
        score_col = input_dictionary.get('score_col', 'Scores')
        # Score for product codes absent from the table. None preserves the
        # historical np.nan behavior; 0.0 is used for sparse DEL read-count
        # libraries where an unmeasured combination is a true non-binder.
        self.default_score = input_dictionary.get('default_score', None)

        # Determine file type and read accordingly
        if ref_filename.endswith('.parquet'):
            ref_df = pl.read_parquet(ref_filename, columns=[compound_col, score_col])
        elif ref_filename.endswith('.csv'):
            ref_df = pl.read_csv(ref_filename, columns=[compound_col, score_col])
        else:
            raise ValueError(f"Unsupported file format: {ref_filename}. Supported formats: .csv, .parquet")

        # Null score cells become NaN so the sampler's NaN-skip path applies.
        scores = [s if s is not None else np.nan for s in ref_df[score_col].to_list()]
        self.ref_dict = dict(zip(ref_df[compound_col].to_list(), scores))

    @property
    def counter(self):
        return self.num_evaluations

    def evaluate(self, product_name):
        self.num_evaluations += 1
        # Return score from lookup. If the product is absent, return the
        # configured default (e.g. 0.0 for sparse read-count libraries) or
        # np.nan when no default is set (sampler skips NaN evaluations).
        missing = self.default_score if self.default_score is not None else np.nan
        return self.ref_dict.get(product_name, missing)

class DBEvaluator(Evaluator):
    """Score = a value looked up by product code in a ``sqlitedict`` database.

    Like :class:`LookupEvaluator` but backed by SQLite, for tables too large
    to hold in memory. Keyed on the product name, so synthesis is skipped.

    Args:
        input_dict: ``{"db_filename": str, "db_prefix": str}`` -- ``db_prefix``
            is prepended to the product name to form the key.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.DBEvaluatorConfig`.
    """

    def __init__(self, input_dictionary):
        self.num_evaluations = 0
        self.db_prefix = input_dictionary['db_prefix']
        db_filename = input_dictionary['db_filename']
        self.ref_dict = SqliteDict(db_filename)

    def __repr__(self):
        return "DBEvalutor"


    @property
    def counter(self):
        return self.num_evaluations


    def evaluate(self, smiles):
        self.num_evaluations += 1
        res = self.ref_dict.get(f"{self.db_prefix}{smiles}")
        if res is None:
            return np.nan
        else:
            if res == -500:
                return np.nan
            return res
    

class FredEvaluator(Evaluator):
    """Score = FRED docking score into a prepared receptor (OpenEye).

    Lower is better -- run with ``mode="minimize"``. Conformers via Omega
    (``max_confs``, default 50; :meth:`set_max_confs`). Slow: use
    ``processes > 1``. Requires the ``openeye`` extra and a licence.

    Args:
        input_dict: ``{"design_unit_file": str}`` -- an ``.oedu`` design unit.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.FredEvaluatorConfig`.

    Raises:
        FileNotFoundError: if the design unit file does not exist.
    """

    def __init__(self, input_dict):
        _ensure_openeye()
        du_file = input_dict["design_unit_file"]
        if not os.path.isfile(du_file):
            raise FileNotFoundError(f"{du_file} was not found or is a directory")
        self.dock = read_design_unit(du_file)
        self.num_evaluations = 0
        self.max_confs = 50

    @property
    def counter(self):
        return self.num_evaluations

    def set_max_confs(self, max_confs):
        """Set the maximum number of conformers generated by Omega
        :param max_confs:
        """
        self.max_confs = max_confs

    def evaluate(self, mol):
        self.num_evaluations += 1
        smi = Chem.MolToSmiles(mol)
        mc_mol = oechem.OEMol()
        oechem.OEParseSmiles(mc_mol, smi)
        confs_ok = generate_confs(mc_mol, self.max_confs)
        score = 1000.0
        docked_mol = oechem.OEGraphMol()
        if confs_ok:
            ret_code = self.dock.DockMultiConformerMolecule(docked_mol, mc_mol)
        else:
            ret_code = oedocking.OEDockingReturnCode_ConformerGenError
        if ret_code == oedocking.OEDockingReturnCode_Success:
            dock_opts = oedocking.OEDockOptions()
            sd_tag = oedocking.OEDockMethodGetName(dock_opts.GetScoreMethod())
            # this is a stupid hack, I need to figure out how to do this correctly
            oedocking.OESetSDScore(docked_mol, self.dock, sd_tag)
            score = float(oechem.OEGetSDData(docked_mol, sd_tag))
        return score

class CustomEvaluator(Evaluator):
    """Score = whatever your Python function returns.

    The simplest way to plug in your own scoring: pass a callable that takes
    an RDKit ``Mol`` and returns a ``float``. Results are cached by canonical
    SMILES; an exception inside the function yields ``NaN`` (the product is
    skipped, the run continues).

    For ``processes > 1`` the callable must be picklable -- a module-level
    function, not a lambda or closure -- because each worker rebuilds the
    evaluator from its config.

    Args:
        scoring_function: ``Callable[[Mol], float]``.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.CustomEvaluatorConfig`.
    """

    def __init__(self, scoring_function):
        self.scoring_function = scoring_function
        self.num_evaluations = 0
        self.score_cache = {}

    @property
    def counter(self):
        return self.num_evaluations

    def evaluate(self, mol):
        """Evaluate a molecule using the user-defined scoring function
        :param mol: Input RDKit molecule
        :return: Float score from the custom function, np.nan on failure
        """      
        self.num_evaluations += 1      
        
        # Look up to see if we already processed this molecule
        smi = Chem.MolToSmiles(mol)
        cached_score = self.score_cache.get(smi)
        if cached_score is not None:
            return cached_score

        try:
            score = float(self.scoring_function(mol))
        except Exception:
            score = np.nan

        self.score_cache[smi] = score
        return score

def generate_confs(mol, max_confs):
    """Generate conformers with Omega
    :param max_confs: maximum number of conformers to generate
    :param mol: input OEMolecule
    :return: Boolean Omega return code indicating success of conformer generation
    """
    _ensure_openeye()
    rms = 0.5
    strict_stereo = False
    omega = oeomega.OEOmega()
    omega.SetRMSThreshold(rms)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetStrictStereo(strict_stereo)
    omega.SetMaxConfs(max_confs)
    error_level = oechem.OEThrow.GetLevel()
    # Turn off OEChem warnings
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)
    status = omega(mol)
    # Turn OEChem warnings back on
    oechem.OEThrow.SetLevel(error_level)
    return status


def read_design_unit(filename):
    """Read an OpenEye design unit
    :param filename: design unit filename (.oedu)
    :return: a docking grid
    """
    _ensure_openeye()
    du = oechem.OEDesignUnit()
    rfs = oechem.oeifstream()
    if not rfs.open(filename):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % filename)

    du = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(rfs, du):
        oechem.OEThrow.Fatal("Failed to read design unit")
    if not du.HasReceptor():
        oechem.OEThrow.Fatal("Design unit %s does not contain a receptor" % du.GetTitle())
    dock_opts = oedocking.OEDockOptions()
    dock = oedocking.OEDock(dock_opts)
    dock.Initialize(du)
    return dock


class MLClassifierEvaluator(Evaluator):
    """Score = positive-class probability from a pickled scikit-learn classifier.

    The model is loaded with ``joblib`` and fed a 2048-bit Morgan fingerprint
    (radius 2); the score is ``predict_proba(...)[:, 1]``.

    Args:
        input_dict: ``{"model_filename": str}`` -- a joblib/pickle file.

    Config: :class:`~TACTICS.thompson_sampling.core.evaluator_config.MLClassifierEvaluatorConfig`.
    """

    def __init__(self, input_dict):
        # joblib was previously imported inside the OpenEye try/except, so it
        # was undefined whenever OpenEye was absent even though it has nothing
        # to do with OpenEye. Import it where it is actually used.
        import joblib

        self.cls = joblib.load(input_dict["model_filename"])
        self.num_evaluations = 0

    @property
    def counter(self):
        return self.num_evaluations

    def evaluate(self, mol):
        self.num_evaluations += 1
        fp = _MORGAN_GEN.GetFingerprint(mol)
        return self.cls.predict_proba([fp])[:,1][0]

