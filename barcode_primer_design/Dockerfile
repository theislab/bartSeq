
# start from plain ubuntu machine
FROM ubuntu
RUN apt-get update; apt-get install -y wget; apt-get install -y bzip2

# install conda
RUN wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda
ENV PATH='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/miniconda/bin'


# get to bioconda
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda


# install primer3 from bioconda
RUN conda install primer3

# install blast (legacy mode) from bioconda
RUN conda install blast-legacy

# install vienna for RNACofold tool from bioconda
RUN conda install viennarna


# install python 2 as last thing as conda thingies override upfront python 2 installation
RUN apt-get install -y python
RUN apt-get install -y python-pip

# pip install required packages for bartender/Primerselect software
RUN pip2 install flask
RUN pip2 install flask-wtf
RUN pip2 install werkzeug
RUN pip2 install biopython
RUN pip2 install reportlab


# now add folders of bartender software

ADD bartender /barcode_primer_design/bartender
ADD databases /barcode_primer_design/databases
ADD primer3_config /barcode_primer_design/primer3_config

# add sofware to pythonpath
ENV PYTHONPATH='$PYTHONPATH:/barcode_primer_design/bartender/'



# start bartender service
ENTRYPOINT ["python", "/barcode_primer_design/bartender/web_frontend/start_service.py"]
