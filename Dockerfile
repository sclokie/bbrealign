
FROM openjdk:11
FROM ubuntu:22.04

RUN apt-get update
#RUN apt-get install -y bash
RUN apt-get install -y mysql-client
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y python3-pip
RUN apt-get install -y libcurl4
RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev # for samtools-1.11 (needed for bbmap)
RUN apt-get install bedtools -y
RUN ln -s /usr/bin/python3 /usr/bin/python

#Install updated crypto library for samtools
RUN wget http://nz2.archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2.19_amd64.deb
RUN dpkg -i libssl1.1_1.1.1f-1ubuntu2.19_amd64.deb

#Install java
RUN apt-get install -y openjdk-11-jdk
RUN apt-get install -y openjdk-11-jre

#Install Miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version
RUN conda install -c bioconda bbmap
RUN conda install -c bioconda picard
###RUN conda install -c bioconda bedtools # Do not use conda - bedtools needs to be >=2.30.0

#Set up App

RUN mkdir /data
RUN mkdir /app
RUN mkdir /genomes
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN pip install "setuptools<58" --upgrade # A way to circumvent a build bug in pyvcf (https://github.com/KarchinLab/open-cravat/issues/98)
RUN pip install pyvcf

#Copy the contents of the current directory to the now working directory: /app
COPY . .

#Gather neded programs from the internet, rather than keep on github
RUN python /app/download_resources.py
# sort sambamba permissions
RUN chmod +x sambamba-0.8.2-linux-amd64-static

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
RUN tar -jvxf samtools-1.11.tar.bz2
WORKDIR samtools-1.11
RUN ./configure --prefix=/app/samtools-1.11
RUN make
RUN make install
RUN ln -s /app/samtools-1.11/samtools /usr/bin/samtools

RUN rm -rf /var/lib/apt/lists/*

WORKDIR /app
ENTRYPOINT ["python", "bbrealign.py"]




