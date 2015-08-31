PKGNAME = mp
PKGVER  = 0.3.1

all: dist

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

.PHONY: attributes clean roxygen
