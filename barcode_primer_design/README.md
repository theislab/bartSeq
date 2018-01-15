
# Run Bartender software for optimal primer design for Bart-seq

## INSTALL

### 1.  Get repository

```
git clone https://github.com/theislab/bartSeq.git
```

### 2. Download fasta databases

As files are too big please download fasta files for blasting

Download human genome Hg38 from here, unpack -> move to `databases` folder. Additionally load databases:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
```
unpack and move to database directory

```
makeblastdb -in hg38.fa -dbtype nucl -title hg38
```

Implemented are options to select the following databases with respective [filename].

- UCSC Genome hg38 [hg38.fa]
- UCSC Genome hg19 [hg19ucsc]
- RefSeq mRNA [refMrna.fa]
- GenBank mRNA [mrna.fa]



### 3. Build docker container

Build docker and run it. Running the container will start the web service reachable at `http://0.0.0.0:5000/primerselect`

```
docker build -t bartender .
docker run --rm -it bartender
```

## EXAMPLE

Go to `http://0.0.0.0:5000/primerselect` and paste sequences of genes where primers shall be designed for in FASTA format.
Find example FASTA file and example of primer plus setting files in `example-input-data`.


## DEVELOPMENT

For development we can mount the code as follows by commenting out adding software files to the docker container. Build it as described above and run container. Then start container service and edit files live from local computer.


```
docker run --rm -it -v /Users/nikola/planA/projects/mini-projects/bartseq\ drukker/manuscript/codePublic/bartSeq/barcode_primer_design:/barcode_primer_design -p 5000:5000 bartender

python /barcode_primer_design/bartender/web_frontend/start_service.py

```
