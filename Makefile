PKGNAME = mp
PKGVER  = 0.4.0

all: dist install

attributes:
	Rscript -e 'library(Rcpp)' -e 'compileAttributes()'

roxygen: DESCRIPTION
	Rscript -e 'library(roxygen2)' -e 'roxygenize()'

dist: roxygen
	R CMD build .

install: dist
	R CMD INSTALL $(PKGNAME)_$(PKGVER).tar.gz

clean:
	rm -f $(PKGNAME)_$(PKGVER).tar.gz

check:
	R CMD check --as-cran $(PKGNAME)_$(PKGVER).tar.gz

.PHONY: attributes clean roxygen
