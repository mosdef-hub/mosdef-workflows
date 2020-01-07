FROM continuumio/miniconda3:4.7.12-alpine

#ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV MOSDEF https://github.com/mosdef-hub
ENV PY_VER 3.7

USER root

# Copy files
COPY . / /tmp/
#ADD . / /tmp/

# Set workdir
WORKDIR /tmp

RUN apk add --no-cache git

#RUN cd gromacs_files \
#    && tar -xvf gromacs-2018.5.tar \
#    && cd gromacs-2018.5 \
#    && mkdir build \
#    && cd build \
#    && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=/tmp \
#    && make \
#    && make check \
#    && make install \
#    && source /tmp/gromacs/bin/GMXRC

RUN conda config --set always_yes yes --set changeps1 no \
    && conda config --add channels omnia \
    && conda config --add channels mosdef \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda update --all -y \
    && conda install -y python=3 conda conda-build pip mbuild foyer matplotlib gromacs\
    && git clone https://github.com/rmatsum836/Pore-Builder.git \
    && cd Pore-Builder \
    && pip install -e . \
    && cd .. \ 
    && conda clean -ya

CMD [ "python", "./graphene-pore.py" ]
