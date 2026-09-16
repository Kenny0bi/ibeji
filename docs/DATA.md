# Data provenance

Every input to ibeji is public. Files were downloaded on 2026-09-14 into `data/raw/` (not committed to git). MD5 checksums below are of the files as downloaded.

| File | Source | Size (bytes) | MD5 |
|---|---|---|---|
| GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz | GEUVADIS E-GEUV-1, [EBI BioStudies](https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/analysis_results/) | 90,812,592 | 500bffeed8e0f770c157e0189e9e50ae (matches EBI's published MD5) |
| E-GEUV-1.sdrf.txt | GEUVADIS sample sheet, EBI BioStudies | 423,201 | 7187e032d87cb1e91aea1d2dc99e22e7 |
| all_phase3.pgen (decompressed from .zst) | 1000 Genomes phase 3, [PLINK 2 resources](https://www.cog-genomics.org/plink/2.0/resources) | 6,696,871,489 | not computed (6.7 GB; the compressed archive was removed after decompression) |
| all_phase3_noannot.pvar.zst | 1000 Genomes phase 3, PLINK 2 resources | 643,457,044 | |
| phase3_corrected.psam | 1000 Genomes phase 3, PLINK 2 resources | 55,174 | da09232eafc7f1fc4c3bddf7c14e8909 |
| iPSYCH-PGC_ASD_Nov2017.gz | [PGC figshare asd2019](https://figshare.com/articles/dataset/asd2019/14671989) | 183,397,420 | 5ca46780db3b37038bd02bd20c38c85c |
| PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz | [PGC figshare scz2022](https://figshare.com/articles/dataset/scz2022/19426775) | 239,710,564 | 6ebe2376f5cda972d37efa0f214c4df0 |
| pgc-bip2021-all.vcf.tsv.gz | [PGC figshare bip2021](https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594) | 385,014,083 | 02e610aaf630e0c869a22fe179d37067 |
| pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz | [PGC figshare MDD2025](https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255) | 232,895,012 | 2b02d6123ecd7fdb12f82314c19c4952 |
| eur_ptsd_pcs_v4_aug3_2021.vcf.gz | [PGC figshare ptsd2024](https://figshare.com/articles/dataset/ptsd2024) (article 26349322) | 254,036,468 | 033685ac31474f59e5a6ec297eca6b49 |

## GWAS details (from file headers)

| Trait | Reference | Cases / controls | Build | Effect column |
|---|---|---|---|---|
| Autism (ASD) | Grove et al. 2019, iPSYCH-PGC | 18,381 / 27,969 | GRCh37 | OR for A1 |
| Schizophrenia (SCZ) | Trubetskoy et al. 2022, PGC3, European core | 52,017 / 75,889 | GRCh37 | BETA for A1 |
| Bipolar disorder (BIP) | Mullins et al. 2021, PGC | 41,917 / 371,549 | GRCh37 | BETA for A1 |
| Major depression (MDD) | PGC MDD2025, European, excluding 23andMe | 412,305 / 1,588,397 | GRCh37 | BETA for EA |
| PTSD | Nievergelt et al. 2024, PGC-PTSD Freeze 3, European | 137,136 / 1,085,746 | GRCh37 | Z for A1 |

Notes:
- The PTSD file is the PGC-PTSD Freeze 3 European GWAS (137,136 cases, 1,222,882 total), the largest public European PTSD GWAS at the time of analysis.
- An earlier PTSD Freeze 2 file was downloaded, then replaced by the Freeze 3 release.
- The SCZ release's terms of use (Fort Lauderdale agreement) permit analyses of this kind; the PGC3 SCZ manuscript has been published.
