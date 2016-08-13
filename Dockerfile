FROM rocker/r-devel:latest

WORKDIR /usr/local/src/mp

COPY install-packages.R /usr/local/src/mp

RUN R CMD BATCH install-packages.R

COPY . /usr/local/src/mp

