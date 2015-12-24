#Snakemake rules:

* **get_msa:** takes a set of contigs as a **.fa** file and gives multiple sequence alignment of the sequence using tool mafft.
```
snakemake get_msa --config ifile=./contigs.fasta ofile=./msa
```

* **run_python_dp:** takes a multiple sequence aligned file as an input and run dynamic program's python implementation.
```
snakemake run_python_dp --config ifile=./msa ofile=./
```

* **run_cpp_dp:** takes a multiple sequence aligned file as an input and run dynamic program's c++11 implementation.
```
snakemake run_cpp_dp --config ifile=./msa ofile=./
```

* **generate_cluster:** Takes sailfish's **.clust** file and divide set of Contigs according to that in different clusters of **.fa** file
```
snakemake generate_cluster --config ifile=quant_human.clust ofile=Trinity.fasta
```

* **makeBlastdb:** Takes following thing as input:
```
1. Transcriptome
2. Corset data's mapping from contigs to cufflinks gene name i.e. (C** -> XLOC**)
3. Cufflink gene name (from **.knownIsoforms** file) to refseq mapping i.e. (XLOC** -> NM***)
4. refseq to gene mapping from Ensemble i.e. (NM**** -> ENS***)
5. Set of All contigs.
6. Initial sub cluster of contigs on which analysis has to be performed
```
This rule parse for all the mapping of contig to gene in the cluster, get most frequent gene and make a blast database of all the transcript of that specific gene.
```
snakemake makeBlastdb
```

* **makeBlastdb_allGene:** Makes Blast database on the whole Transcriptome based on gene bags. **Deprecated**
```
snakemake makeBlastdb_allGene
```

* **run_blast:** Given blast database and MSA of to be queried contigs of one cluster, this rule run our DP algorithm get exon boudaries and writes a consensus file i.e. **.fa** file with sequence of exons. After that it Blast this set of exons to our given database to make allvsall eval matrix.
```
snakemake run_blast
```

* **run_blast_allGenes** run Blast on all the database created. (i.e. eventually a Transcriptome) **Deprecated**
```
snakemake run_blast_allGenes
```


* **get_similarity** Given allvsall eval Mtrix score, run Hungarian algorithm on that to get Bipartite matching.
```
snakemake get_similarity
```
