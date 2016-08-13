NAME = mp
VERSION = 0.4.0
ARTIFACT=${NAME}_${VERSION}.tar.gz

all: build install check dist

clean:
	find . -name $(ARTIFACT) | xargs rm -f

attributes:
	Rscript -e 'library(Rcpp)' -e 'compileAttributes()'

roxygen: DESCRIPTION
	Rscript -e 'library(roxygen2)' -e 'roxygenize()'

build: clean roxygen attributes
	R CMD build .

install:
	R CMD INSTALL $(ARTIFACT)

check:
	R CMD check --as-cran $(ARTIFACT)

dist:
	cd dist; R CMD build ..

.PHONY: attributes clean roxygen dist
