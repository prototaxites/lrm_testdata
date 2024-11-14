# generating test data for hifi + hic long read metagenomics

create conda environment
```
conda create -n mg_test_data nextflow pip samtools=1.21 minimap2=2.28 seqkit
pip install git+https://github.com/cerebis/sim3C
sitedir=`python -c 'import site; print(site.getsitepackages())' | tr -d []\'`
sed 's#profile_path = os.path.join(os.path.dirname(args.output_file), args.profile_name)#profile_path = os.path.join(os.path.dirname(args.output_file_1), args.profile_name)#g' $sitedir/sim3C/command_line.py > $sitedir/sim3C/command_line.py
```

first, spit out SRA id to file

```
echo "SRR13128014" > SRA.csv
```

then download using nf-core/fetchngs

```
nextflow run nf-core/fetchngs -r 1.12.0 -profile singularity,sanger --input SRA.csv --outdir zymo_fq
```

download the zymo reference genomes

```
wget https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip
unzip D6331.refseq.zip
```

map reads against the fungal genomes to remove them - wouldn't be a problem but the genomes have chromosomes which means the hi-c test binning data will be problematic

```
mkdir fungal_refs
cat D6331.refseq/genomes/Candida_albican.fasta D6331.refseq/genomes/Saccharomyces_cerevisiae.fasta > fungal_refs/fungal_ref.fa
rm D6331.refseq/genomes/Candida_albican.fasta D6331.refseq/genomes/Saccharomyces_cerevisiae.fasta
minimap2 -ax map-hifi fungal_refs/fungal_ref.fa zymo_fq/fastq/SRX9569057_SRR13128014.fastq.gz -o fungal_refs/map.bam
samtools fasta -f 4 -0 fungal_refs/non_fungal_reads.fa.gz fungal_refs/map.bam
```

subsample reads to get small input fq file

```
seqkit sample -p 0.1 fungal_refs/non_fungal_reads.fa.gz | seqkit head -n 20000 -o subsampled.pacbio.fa.gz
seqkit split subsampled.pacbio.fa.gz -p 2 -O final_data
```

simulate hi-c data

```
cat D6331.refseq/genomes/*.fasta > all.fasta
sim3C --dist uniform -n 250000 -l 150 -e DpnII -e HinfI -v -m hic all.fasta hic_1.fq.gz hic_2.fq.gz
samtools import -1 hic_1.fq.gz -2 hic_2.fq.gz -o final_data/hic.cram
```

## assembly

```
singularity exec --bind /lustre:/lustre /software/team311/jd42/images/metamdbg-1.0.sif metaMDBG asm --out-dir test_asm --threads 12 --in-hifi final_data/subsampled.pacbio.part_001.fa.gz final_data/subsampled.pacbio.part_002.fa.gz
```

```
	Run time:                   3min 20sec
	Peak memory:                0.555568 GB
	Assembly length:            28849185
	Contigs N50:                289133
	Nb contigs:                 586
	Nb Contigs (>1Mb):          6
	Nb circular contigs (>1Mb): 1
```

## binning
```
cd test_asm
mkdir bin
minimap2 -ax map-hifi contigs.fasta.gz ../subsampled.pacbio.fa.gz | samtools sort -o map.bam
singularity exec /software/team311/ng13/local/images/metabat2-2.15--c1941c7.sif jgi_summarize_bam_contig_depths --outputDepth depth.txt map.bam
singularity exec /software/team311/ng13/local/images/metabat2-2.15--c1941c7.sif metabat2 -i contigs.fasta.gz -o bins/bin -a depth.txt
```

```
processed files:  18 / 18 [======================================] ETA: 0s. done
file            format  type  num_seqs    sum_len    min_len    avg_len    max_len
bins/bin.10.fa  FASTA   DNA         59    905,965      3,204   15,355.3     57,634
bins/bin.11.fa  FASTA   DNA          7  2,782,460     16,104  397,494.3  1,976,763
bins/bin.12.fa  FASTA   DNA        110  1,385,011      2,543     12,591     66,404
bins/bin.13.fa  FASTA   DNA         41    942,067      2,890   22,977.2     84,292
bins/bin.14.fa  FASTA   DNA         19    256,095      5,998   13,478.7     23,269
bins/bin.15.fa  FASTA   DNA          1    206,329    206,329    206,329    206,329
bins/bin.16.fa  FASTA   DNA         47  2,111,652      5,320   44,928.8    207,553
bins/bin.17.fa  FASTA   DNA          5  1,524,841     28,989  304,968.2    706,307
bins/bin.18.fa  FASTA   DNA          7    211,469      8,510   30,209.9     71,994
bins/bin.1.fa   FASTA   DNA          5    218,059      8,986   43,611.8    146,796
bins/bin.2.fa   FASTA   DNA          5    261,644     19,457   52,328.8     92,618
bins/bin.3.fa   FASTA   DNA         19    223,814      2,674   11,779.7     47,822
bins/bin.4.fa   FASTA   DNA         21    357,878      8,070   17,041.8     62,966
bins/bin.5.fa   FASTA   DNA          6  1,017,529     23,335  169,588.2    361,104
bins/bin.6.fa   FASTA   DNA         42  4,874,765      5,161  116,065.8  1,155,020
bins/bin.7.fa   FASTA   DNA          1  2,158,040  2,158,040  2,158,040  2,158,040
bins/bin.8.fa   FASTA   DNA          6  5,316,910     61,781  886,151.7  1,554,009
bins/bin.9.fa   FASTA   DNA         68  2,613,812      2,901   38,438.4    161,719
```
## hic binning

```
singularity exec -B /lustre/ /software/tola/images/bwa-mem2-2.2.1.sif bwa-mem2 mem -t 12 contigs.fasta.gz ../hic_1.fq.gz ../hic_2.fq.gz | samtools view -F 0x904 -o hic_map_unsorted.bam
samtools sort -n hic_map_unsorted.bam -o hic_map.bam
singularity exec -B /lustre /software/team311/ng13/local/images/bin3c.sif bin3C mkmap -e DpnII -e HinfI -v contigs.fasta.gz hic_map.bam bin3c_out
zcat contigs.fasta.gz | bgzip -c > contigs.fasta.bgz
singularity exec -B /lustre /software/team311/ng13/local/images/bin3c.sif bin3C cluster -v --fasta contigs.fasta.bgz bin3c_out/contact_map.p.gz bin3c_clust
```

```
seqkit stats bin3c_clust/fasta/*
processed files:  12 / 12 [======================================] ETA: 0s. done
file                        format  type  num_seqs    sum_len    min_len    avg_len    max_len
bin3c_clust/fasta/CL01.fna  FASTA   DNA          8  5,588,595     61,781  698,574.4  1,554,009
bin3c_clust/fasta/CL02.fna  FASTA   DNA         69  5,521,186      3,378   80,017.2  1,155,020
bin3c_clust/fasta/CL03.fna  FASTA   DNA          9  2,916,882     16,104    324,098  1,976,763
bin3c_clust/fasta/CL04.fna  FASTA   DNA         59  2,815,142      5,785   47,714.3    161,719
bin3c_clust/fasta/CL05.fna  FASTA   DNA         10  2,477,014     23,335  247,701.4    706,307
bin3c_clust/fasta/CL06.fna  FASTA   DNA          1  2,158,040  2,158,040  2,158,040  2,158,040
bin3c_clust/fasta/CL07.fna  FASTA   DNA         20  1,411,387     19,517   70,569.4    207,553
bin3c_clust/fasta/CL08.fna  FASTA   DNA         30  1,106,370     12,193     36,879     89,664
bin3c_clust/fasta/CL09.fna  FASTA   DNA          3    157,597     31,765   52,532.3     68,198
bin3c_clust/fasta/CL10.fna  FASTA   DNA          2     79,040     12,636     39,520     66,404
bin3c_clust/fasta/CL11.fna  FASTA   DNA          3     78,627      6,177     26,209     51,822
bin3c_clust/fasta/CL12.fna  FASTA   DNA          3     74,434     19,447   24,811.3     34,960
```
