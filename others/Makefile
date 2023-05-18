# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
BUILDDIR      = build
HTMLBUILDDIR  = build/html

ALLSPHINXOPTS   = $(SPHINXOPTS) source
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(SPHINXOPTS) source


.PHONY: view clean html

all: html

view:
	open $(HTMLBUILDDIR)/index.html

clean:
	-rm -rf $(HTMLBUILDDIR) $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(HTMLBUILDDIR)
	@echo
	@echo "Build finished. The HTML pages are in $(HTMLBUILDDIR)."
