import hail as hl

twb1465_vaf10_mt = hl.import_vcf('/Users/xuanchou/Documents/TWB_1492/add_annotation/output_age/sample_vcf.vcf.bgz', reference_genome='GRCh38')
twb1465_vaf10 = twb1465_vaf10_mt.select_rows(twb1465_vaf10_mt.info.AF_hom, twb1465_vaf10_mt.info.AF_het)
twb1465_vaf10 = twb1465_vaf10.rename({'AF_het': 'AF_het_vaf10'})

twb1465_vaf05_mt = hl.import_vcf('/Users/xuanchou/Documents/TWB_1492/add_annotation/output_vaf0.05/sample_vcf.vcf.bgz', reference_genome='GRCh38')
twb1465_vaf05 = twb1465_vaf05_mt.select_rows(twb1465_vaf05_mt.info.AF_hom, twb1465_vaf05_mt.info.AF_het)
twb1465_vaf05 = twb1465_vaf05.rename({'AF_het': 'AF_het_vaf05'})

# Combine AF info
twb1465_combined = twb1465_vaf05.annotate_rows(
    AF_het_vaf10 = twb1465_vaf10.index_rows(twb1465_vaf05.locus, twb1465_vaf05.alleles).AF_het_vaf10)

# Get filter info from original dataset
# Aggregate af info to a new row called info
twb1465_combined = twb1465_combined.annotate_rows(
    rsid = twb1465_vaf05_mt.index_rows(twb1465_combined.locus, twb1465_combined.alleles).rsid,
    qual = '.',
    filters = twb1465_vaf10_mt.index_rows(twb1465_combined.locus, twb1465_combined.alleles).filters,
    info = hl.struct(AF_hom=twb1465_combined.AF_hom, 
                    AF_het_vaf05=twb1465_combined.AF_het_vaf05, 
                    AF_het_vaf10=twb1465_combined.AF_het_vaf10)
)
# Drop intermediate rows
twb1465_combined = twb1465_combined.drop("AF_hom", "AF_het_vaf05", "AF_het_vaf10")

# Output custom VCF
vcf_variant_ht = twb1465_combined.rows()
rows_mt = hl.MatrixTable.from_rows_table(vcf_variant_ht).key_cols_by(s="foo")
hl.export_vcf(rows_mt, '/Users/xuanchou/Documents/TWB_1492/custom_annotation/twb1465_af.vcf.bgz')