FROM r-base

RUN apt-get update && \
  apt-get install texlive -y

WORKDIR /usr/local/src/mp

COPY install-packages.R /usr/local/src/mp

RUN R CMD BATCH install-packages.R

COPY . /usr/local/src/mp

