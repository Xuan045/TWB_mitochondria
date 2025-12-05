import hail as hl

outdir = '/Users/xuanchou/github/Mitochondria_pipeline/custom_annotation'
gnomad_mt = hl.import_vcf("/Users/xuanchou/Documents/gnomad_mt_doc/gnomad.genomes.v3.1.sites.chrM.vcf.bgz", reference_genome='GRCh38')
info = gnomad_mt['info']

# Select AF columns
select_fields = ['AF_hom', 'AF_het', 'pop_AF_hom', 'pop_AF_het']
info = info.select(*select_fields)

# Split pop_AF by |
pop_AF_hom_split = info['pop_AF_hom'].split('\|')
pop_AF_het_split = info['pop_AF_het'].split('\|')

# Extract AF for each population and then store to a new field
info = info.annotate(
    AF_hom_afr = hl.float(pop_AF_hom_split[0]),
    AF_het_afr = hl.float(pop_AF_het_split[0]),
    AF_hom_ami = hl.float(pop_AF_hom_split[1]),
    AF_het_ami = hl.float(pop_AF_het_split[1]),
    AF_hom_amr = hl.float(pop_AF_hom_split[2]),
    AF_het_amr = hl.float(pop_AF_het_split[2]),
    AF_hom_asj = hl.float(pop_AF_hom_split[3]),
    AF_het_asj = hl.float(pop_AF_het_split[3]),
    AF_hom_eas = hl.float(pop_AF_hom_split[4]),
    AF_het_eas = hl.float(pop_AF_het_split[4]),
    AF_hom_sas = hl.float(pop_AF_hom_split[8]),
    AF_het_sas = hl.float(pop_AF_het_split[8]),
    AF_hom_fin = hl.float(pop_AF_hom_split[5]),
    AF_het_fin = hl.float(pop_AF_het_split[5]),
    AF_hom_nfe = hl.float(pop_AF_hom_split[6]),
    AF_het_nfe = hl.float(pop_AF_het_split[6]),
    AF_hom_mid = hl.float(pop_AF_hom_split[9]),
    AF_het_mid = hl.float(pop_AF_het_split[9]),
    AF_hom_oth = hl.float(pop_AF_hom_split[7]),
    AF_het_oth = hl.float(pop_AF_het_split[7])
)

# Drop the intermediate fields
info = info.drop("pop_AF_hom", "pop_AF_het")

# Combine with original table
# Reset 'qual' to "." since the quality in output VCF is weirdly assigned
gnomad_mt = gnomad_mt.annotate_rows(info=info, qual='.')
hl.export_vcf(gnomad_mt, f'{outdir}/gnomad_af.vcf.bgz')
