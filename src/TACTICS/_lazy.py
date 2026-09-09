"""PEP 562 lazy re-exports for package ``__init__`` modules.

The package hubs (``TACTICS``, ``TACTICS.thompson_sampling``, ``.core``,
``TACTICS.library_analysis``) re-export names from their submodules for
convenience. Importing those submodules eagerly makes ``import TACTICS`` --
and any config-only import that passes through a hub -- pay for RDKit,
scipy, sqlitedict and the plotting stack up front. Each hub instead declares
a ``{name: ".submodule"}`` map and calls :func:`install`, which resolves a
name on first attribute access and caches it in the module namespace.

Submodule paths (``TACTICS.thompson_sampling.core.evaluators``) are
unaffected: a lazy hub is still a regular package.
"""

from importlib import import_module
from types import ModuleType
from typing import Dict


def install(module: ModuleType, exports: Dict[str, str]) -> None:
    """Give ``module`` a lazy ``__getattr__``/``__dir__`` over ``exports``.

    Parameters:
        module: The package module (pass ``sys.modules[__name__]``).
        exports: ``{public_name: relative_submodule}`` -- e.g.
            ``{"ThompsonSampler": ".core.sampler"}``.
    """
    package = module.__name__

    def __getattr__(name: str):
        try:
            submodule = exports[name]
        except KeyError:
            raise AttributeError(f"module {package!r} has no attribute {name!r}") from None
        value = getattr(import_module(submodule, package), name)
        setattr(module, name, value)  # cache; __getattr__ is not called again
        return value

    def __dir__():
        return sorted(set(vars(module)) | set(exports))

    module.__getattr__ = __getattr__
    module.__dir__ = __dir__
    module.__all__ = sorted(exports)
