#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                      SelectVariants: SNPs && PASS                  "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// fasta
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; ref_mutect2_tum_only_mode_channel ; ref_for_create_GenomicsDB_channel ; ref_create_somatic_PoN ; fasta_mutect; fasta_variant_eval ; fasta_filter_mutect_calls ; fasta_vcf2maf ; fasta_select_variants_PASS }
}
// fai
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { fai_mutect; fai_baserecalibrator; fai_haplotypecaller; ref_index_mutect2_tum_only_mode_channel ; ref_index_for_create_GenomicsDB_channel ; ref_index_create_somatic_PoN ; fai_variant_eval ; fai_filter_mutect_calls ; fai_vcf2maf ; fai_select_variants_PASS}
}

// dict
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_interval; dict_mutect; dict_baserecalibrator; dict_haplotypecaller; dict_variant_eval ; ref_dict_mutect2_tum_only_mode_channel ; ref_dict_for_create_GenomicsDB_channel ; ref_dict_create_somatic_PoN ; dict_filter_mutect_calls ; dict_vcf2maf ; dict_select_variants_PASS}
}

Channel
    .fromPath("${params.inputdir}/*.vcf")
    .into {  vcf_filtered_for_select_variants}

Channel
    .fromPath("${params.inputdir}/*.idx")
    .into {  idx_vcf_filtered_for_select_variants}

// SelectVariants only the ones with PASS:
// Found here: 
// https://gatkforums.broadinstitute.org/gatk/discussion/1742/using-selectvariants-to-output-pass-records
// Should have done with grep but let's gatk much, shall we
//     -select 'vc.isNotFiltered()' : only the ones with flag PASS
//     -select-type SNP             : only SNPs

process SelectSNPsPASS {

    tag "${filtered_vcf}"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/SelectedSomaticSNPs_VCF", mode: 'copy'

    input:
    file(filtered_vcf) from vcf_filtered_for_select_variants
    file(unfiltered_vcf_idx) from idx_vcf_filtered_for_select_variants
    each file(fasta) from fasta_select_variants_PASS
    each file(fai) from fai_select_variants_PASS
    each file(dict) from dict_select_variants_PASS

    output:
    file("*vcf") into vcf_SNPs_PASS_for_vcf2maf
    file("*vcf.idx") into idx_vcf_SNPs_PASS_for_vcf2maf
 
    script:
    """
    gatk FilterMutectCalls \
    -R ${fasta} \
    -V $filtered_vcf \
    -O "${filtered_vcf.simpleName}.passed.SNPs.vcf"
    -select 'vc.isNotFiltered()' 
    -select-type SNP
   """
}

process Vcf2maf {

    tag "${passed_SNPs}"
    container 'levim/vcf2maf:1.0'
    publishDir "${params.outdir}/SelectedSomaticSNPs_MAF", mode: 'copy'

    input:
    file(vcf_passed_SNPs) from vcf_SNPs_PASS_for_vcf2maf
    file(idx_vcf_passed_SNPs) from idx_vcf_SNPs_PASS_for_vcf2maf
    each file(fasta) from fasta_vcf2maf
    each file(fai) from fai_vcf2maf
    each file(dict) from dict_vcf2maf

    output:
    file("*") into vcf2maf_annotated_files_channel
 
    script:
    """
    basename=\$(echo ${vcf_passed_SNPs.simpleName})
    tumourID=\$(echo \$basename | cut -f 1 -d '_')
    normalID=\$(echo \$basename | cut -f 4 -d '_')

    `echo "tumourID" $tumourID`
    `echo "normalID" $normalID`

    perl /opt/vcf2maf/vcf2maf.pl \
    --input-vcf $vcf_passed_SNPs \
    --output-maf "\${basename}.maf"  \
    --tumor-id \${tumourID} \
    --normal-id \${normalID} \
    --ref-fasta /vepdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --ncbi-build  GRCh37 \
    --filter-vcf /vepdata/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --vep-path /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vepdata/ \
    --vep-forks 2 \
    --buffer-size 200 \
    --species homo_sapiens \
    --cache-version 89
    """
}