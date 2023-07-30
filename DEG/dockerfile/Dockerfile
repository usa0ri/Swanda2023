FROM rocker/rstudio:4.2.1

RUN useradd -m saori
ENV DEBIAN_FRONTEND=noninteractive

#Install all dependencies and necessary R packages for Seurat.
RUN R -e "install.packages(c('BiocManager'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('edgeR')"

# additional R packages
RUN R -e "install.packages(c('tidyverse','tictoc','R.utils'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

# ####################################
# COPY entrypoint.sh /
# RUN chmod +x entrypoint.sh
# WORKDIR /home/rstudio
# ENTRYPOINT ["/entrypoint.sh"]
