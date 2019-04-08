FROM continuumio/miniconda3:4.5.12

# get to bioconda
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# install packages from bioconda
RUN conda install primer3 blast-legacy viennarna

# install python 2 as last thing as conda thingies override upfront python 2 installation
RUN apt-get update
RUN apt-get install -y python python-pip

# pip install required packages for bartender/Primerselect software
COPY bartender/requirements.txt /barcode_primer_design/
RUN pip2 install -r /barcode_primer_design/requirements.txt

# now add folders of bartender software
COPY bartender /barcode_primer_design/bartender
COPY databases /barcode_primer_design/databases
COPY primer3_config /barcode_primer_design/primer3_config
# sanity check
RUN bash -c 'ls /barcode_primer_design/databases/{hg19ucsc,{hg38,mrna,refMrna}.fa}.{nhr,nin,nsq} >/dev/null 2>&1'

# add sofware to pythonpath
ENV PYTHONPATH='$PYTHONPATH:/barcode_primer_design/bartender/'

# start bartender service
ENTRYPOINT ["python2", "/barcode_primer_design/bartender/web_frontend/start_service.py"]
