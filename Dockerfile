FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN apt-get update --fix-missing
RUN apt-get install -y gcc wget ca-certificates samtools && rm -rf /var/lib/apt/lists/*

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip \
    && python --version

RUN pip install coverage==5.5 && \
    pip install in_place==1.0.1 && \
    pip install pathos==0.3.4

# download StringTie software and untar it
RUN STRINGTIE_V='2.1.4' && \
    mkdir -p /kb/deployment/bin/modules && \
    cd /kb/deployment/bin/modules && \
    mkdir -p StringTie && cd StringTie && \
    wget https://ccb.jhu.edu/software/stringtie/dl/stringtie-${STRINGTIE_V}.Linux_x86_64.tar.gz && \
    tar xvfz stringtie-${STRINGTIE_V}.Linux_x86_64.tar.gz && \
    cd stringtie-${STRINGTIE_V}.Linux_x86_64 && \
    mkdir -p /kb/deployment/bin/StringTie && \
    cp -R stringtie /kb/deployment/bin/StringTie/stringtie && \
    /kb/deployment/bin/StringTie/stringtie --version

# -----------------------------------------

# download prepDE script
RUN cd /kb/deployment/bin/modules && \
    mkdir -p prepDE && cd prepDE && \
    wget https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py3 && \
    mkdir -p /kb/deployment/bin/prepDE && \
    cp -R prepDE.py3 /kb/deployment/bin/prepDE/prepDE.py && \
    chmod 755 /kb/deployment/bin/prepDE/prepDE.py

# -----------------------------------------

# download gffread script
RUN GFFREAD_V='0.12.7' && \
    cd /kb/deployment/bin/modules && \
    mkdir -p gffread && cd gffread && \
    wget https://ccb.jhu.edu/software/stringtie/dl/gffread-${GFFREAD_V}.Linux_x86_64.tar.gz && \
    tar xvfz gffread-${GFFREAD_V}.Linux_x86_64.tar.gz && \
    cd gffread-${GFFREAD_V}.Linux_x86_64 && \
    mkdir -p /kb/deployment/bin/gffread && \
    cp -R gffread /kb/deployment/bin/gffread/gffread && \
    chmod 755 /kb/deployment/bin/gffread/gffread

# -----------------------------------------

# download gffcompare script
RUN GFFCOMP_V='0.12.10' && \
    cd /kb/deployment/bin/modules && \
    mkdir -p gffcompare && cd gffcompare && \
    wget https://ccb.jhu.edu/software/stringtie/dl/gffcompare-${GFFCOMP_V}.Linux_x86_64.tar.gz && \
    tar xvfz gffcompare-${GFFCOMP_V}.Linux_x86_64.tar.gz && \
    cd gffcompare-${GFFCOMP_V}.Linux_x86_64 && \
    mkdir -p /kb/deployment/bin/gffcompare && \
    cp -R gffcompare /kb/deployment/bin/gffcompare/gffcompare && \
    chmod 755 /kb/deployment/bin/gffcompare/gffcompare

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
