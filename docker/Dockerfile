FROM rocker/drd

RUN apt-get update && apt-get install \
	-y --no-install-recommends \
	libssl-dev
RUN Rscript -e 'utils::install.packages(c("ggplot2"))'

CMD ["bash"]

