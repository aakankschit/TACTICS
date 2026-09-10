# Sphinx build for the TACTICS docs.
#
# `make html` mirrors CI: warnings are errors. Requires the package installed
# with its docs extras:  pip install -e ".[docs]"

SPHINXBUILD ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = build
SPHINXOPTS  ?= -W --keep-going

.PHONY: help html clean linkcheck livehtml

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

html:
	$(SPHINXBUILD) $(SPHINXOPTS) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

linkcheck:
	$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck"

clean:
	rm -rf "$(BUILDDIR)"

# Auto-rebuilding local server (pip install sphinx-autobuild)
livehtml:
	sphinx-autobuild $(SPHINXOPTS) "$(SOURCEDIR)" "$(BUILDDIR)/html" --open-browser
