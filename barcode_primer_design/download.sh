#!/bin/bash

rsync -aP heidi:/home/icb/steffen.sass/data/{hg19ucsc,{hg38,mrna,refMrna}.fa}.{nhr,nin,nsq} ./databases/ 
