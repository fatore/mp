PKGNAME = mp
PKGVER  = 0.3

all: dist

dist:
	Rscript -e 'library(roxygen2)' -e 'roxygenize()'
	R CMD build .

install: dist
	R CMD INSTALL $(PKGNAME)_$(PKGVER).tar.gz

clean:
	rm -f $(PKGNAME)_$(PKGVER).tar.gz

.PHONY: clean
