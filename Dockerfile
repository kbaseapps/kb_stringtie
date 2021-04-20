FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN apt-get update --fix-missing
RUN apt-get install -y gcc wget

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip \
    && python --version

RUN pip install coverage==5.5 && \
    pip install in_place==0.5.0 && \
    pip install pathos==0.2.7

# download StringTie software and untar it
RUN mkdir /kb/deployment/bin/modules && \
    cd /kb/deployment/bin/modules && \
    mkdir StringTie && cd StringTie && \
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz &&\
    tar xvfz stringtie-1.3.3b.Linux_x86_64.tar.gz && \
    cd stringtie-1.3.3b.Linux_x86_64 && \
    mkdir /kb/deployment/bin/StringTie && \
    cp -R stringtie /kb/deployment/bin/StringTie/stringtie

# -----------------------------------------

# download prepDE script
RUN cd /kb/deployment/bin/modules && \
    mkdir prepDE && cd prepDE && \
    wget http://ccb.jhu.edu/software/stringtie/dl/prepDE.py &&\
    mkdir /kb/deployment/bin/prepDE && \
    cp -R prepDE.py /kb/deployment/bin/prepDE/prepDE.py && \
    chmod 777 /kb/deployment/bin/prepDE/prepDE.py

# -----------------------------------------

# download gffread script
RUN cd /kb/deployment/bin/modules && \
    mkdir gffread && cd gffread && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.6.Linux_x86_64.tar.gz &&\
    tar xvfz gffread-0.12.6.Linux_x86_64.tar.gz && \
    cd gffread-0.12.6.Linux_x86_64 && \
    mkdir /kb/deployment/bin/gffread && \
    cp -R gffread /kb/deployment/bin/gffread/gffread

# -----------------------------------------

# download gffcompare script
RUN cd /kb/deployment/bin/modules && \
    mkdir gffcompare && cd gffcompare && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.4.Linux_x86_64.tar.gz	 &&\
    tar xvfz gffcompare-0.10.4.Linux_x86_64.tar.gz	 && \
    cd gffcompare-0.10.4.Linux_x86_64	 && \
    mkdir /kb/deployment/bin/gffcompare && \
    cp -R gffcompare /kb/deployment/bin/gffcompare/gffcompare

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
