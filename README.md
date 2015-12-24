# txome-clustering

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/txome-clustering?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
Meaningful and efficient clustering of de novo transcriptome assembly results (better name pending).

Snakemake rules:

get_msa: takes a set of contigs as a **.fa** file and gives multiple sequence alignment of the sequence using tool mafft.
```
snakemake get_msa --config ifile=./contigs.fasta ofile=./msa
```

run_python_dp: takes a multiple sequence aligned file as an input and run dynamic program's python implementation.
```
snakemake run_python_dp --config ifile=./msa ofile=./
```

run_cpp_dp: takes a multiple sequence aligned file as an input and run dynamic program's c++11 implementation.
```
snakemake run_cpp_dp --config ifile=./msa ofile=./
```

generate_cluster: Takes sailfish's **.clust** file and divide Transcriptome according to that in different clusters of **.fa** file
```
snakemake generate_cluster --config ifile=quant_human.clust ofile=Trinity.fasta
```

makeBlastdb : Takes following thing as input:
1. Transcriptome
2. Corset data's mapping from contigs to cufflinks gene name i.e. (C** -> XLOC**)
3. 
