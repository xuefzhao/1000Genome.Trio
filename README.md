# 1000Genome.Trio
We applied SVelter to YRI trio first, then we validated our call sets using pacbio long sequences, and analyzed the common and unique SVs.

###Analysis pipeline:
1. Run SVelter setup, NullModel, BPSearch on YRI trio

2. Integrate breakpoints from each sample by running: 
`<BP.Merge.Multi.Sample.py --chromosome chr1 --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/ --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa>`

3.Run SVelter SVPredict, SVIntegrate on TRI trio
`<Pacbio.Validator.py vcf --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs/NA19240.Common.BPs.vcf --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs/NA19240.vcf.PacbioValidation --pacbio-input /scratch/remills_flux/xuefzhao/SV_discovery_index/smrt.download/alignment/NA19240.XXX.bam --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa>`

#To extract common and unique SVs from trio:
`<SV.Merge.Multi.Sample.py vcf --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs>`
`<SV.Merge.Multi.Sample.py svelter --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs>`





