Extending
=========

Base classes and mixins for new strategies, warmups and evaluators. See the
:doc:`guide </extending>`.

.. autoclass:: TACTICS.thompson_sampling.strategies.base_strategy.SelectionStrategy
   :members:

.. autoclass:: TACTICS.thompson_sampling.strategies._thermal.ThermalCyclingMixin
   :members:
   :private-members: _init_thermal_cycling, _rotation_flexibility

.. autoclass:: TACTICS.thompson_sampling.strategies._thermal.GMICCriticalityMixin
   :members:
   :private-members: _init_gmic_state, _calculate_gmic, _calculate_gmic_details
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.warmup.base.WarmupStrategy
   :members:

Package layout
--------------

.. automodule:: TACTICS._lazy
   :members: install
