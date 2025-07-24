Library Enumeration Module
=========================

The Library Enumeration module provides tools for generating and managing chemical libraries.

LibraryEnumerator Class
----------------------

.. autoclass:: PRISMS.library_enumeration.generate_products.LibraryEnumerator
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: enumerate_library
   .. automethod:: get_product_smiles

Utility Functions
----------------

.. autofunction:: PRISMS.library_enumeration.enumeration_utils.find_reactants_from_product_code
.. autofunction:: PRISMS.library_enumeration.enumeration_utils.write_products_to_files 