version 1.0

workflow eager {
	input{
		File samplesheet
		String udg_type = "none"
		Boolean? single_stranded
		Boolean? single_end
		Int colour_chemistry = 4
		Boolean? bam
		String? snpcapture_bed
		Boolean? run_convertinputbam
		File? fasta
		String? genome
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		String? bwa_index
		String? bt2_index
		String? fasta_index
		String? seq_dict
		Boolean? large_ref
		Boolean? save_reference
		String outdir = "./results"
		String publish_dir_mode = "copy"
		Boolean? help
		Boolean validate_params = true
		String? email
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean? show_hidden_params
		Boolean? enable_conda
		String schema_ignore_params = "genomes"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		String? awsqueue
		String awsregion = "eu-west-1"
		String? awscli
		Boolean? skip_fastqc
		Boolean? skip_adapterremoval
		Boolean? skip_preseq
		Boolean? skip_deduplication
		Boolean? skip_damage_calculation
		Boolean? skip_qualimap
		Boolean? complexity_filter_poly_g
		Int complexity_filter_poly_g_min = 10
		String clip_forward_adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
		String clip_reverse_adaptor = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
		String? clip_adapters_list
		Int clip_readlength = 30
		Int clip_min_read_quality = 20
		Int min_adap_overlap = 1
		Boolean? skip_collapse
		Boolean? skip_trim
		Boolean? preserve5p
		Boolean? mergedonly
		Int qualitymax = 41
		Boolean? run_post_ar_trimming
		Int post_ar_trim_front = 7
		Int post_ar_trim_tail = 7
		Int post_ar_trim_front2 = 7
		Int post_ar_trim_tail2 = 7
		String mapper = "bwaaln"
		Float bwaalnn = 0.01
		Int bwaalnk = 2
		Int bwaalnl = 1024
		Int bwaalno = 2
		Int circularextension = 500
		String circulartarget = "MT"
		Boolean? circularfilter
		String bt2_alignmode = "local"
		String bt2_sensitivity = "sensitive"
		Int? bt2n
		Int? bt2l
		Int? bt2_trim5
		Int? bt2_trim3
		Int bt2_maxins = 500
		Boolean? hostremoval_input_fastq
		String hostremoval_mode = "remove"
		Boolean? run_bam_filtering
		Int? bam_mapping_quality_threshold
		Int? bam_filter_minreadlength
		String bam_unmapped_type = "discard"
		String dedupper = "markduplicates"
		Boolean? dedup_all_merged
		String preseq_mode = "c_curve"
		Int preseq_step_size = 1000
		Int preseq_maxextrap = 10000000000
		Int preseq_terms = 100
		Int preseq_bootstrap = 100
		Float preseq_cval = 0.95
		Int damageprofiler_length = 100
		Int damageprofiler_threshold = 15
		Float damageprofiler_yaxis = 0.3
		Boolean? run_pmdtools
		Int pmdtools_range = 10
		Int pmdtools_threshold = 3
		String? pmdtools_reference_mask
		Int pmdtools_max_reads = 10000
		Boolean? pmdtools_platypus
		Boolean? run_mapdamage_rescaling
		Int rescale_length_5p = 12
		Int rescale_length_3p = 12
		Boolean? run_bedtools_coverage
		String? anno_file
		Boolean? run_trim_bam
		Int? bamutils_clip_double_stranded_half_udg_left
		Int? bamutils_clip_double_stranded_half_udg_right
		Int? bamutils_clip_double_stranded_none_udg_left
		Int? bamutils_clip_double_stranded_none_udg_right
		Int? bamutils_clip_single_stranded_half_udg_left
		Int? bamutils_clip_single_stranded_half_udg_right
		Int? bamutils_clip_single_stranded_none_udg_left
		Int? bamutils_clip_single_stranded_none_udg_right
		Boolean? bamutils_softclip
		Boolean? run_genotyping
		String? genotyping_tool
		String genotyping_source = "raw"
		Int gatk_call_conf = 30
		Int gatk_ploidy = 2
		Int gatk_downsample = 250
		String? gatk_dbsnp
		String gatk_hc_out_mode = "EMIT_VARIANTS_ONLY"
		String gatk_hc_emitrefconf = "GVCF"
		String gatk_ug_out_mode = "EMIT_VARIANTS_ONLY"
		String gatk_ug_genotype_model = "SNP"
		Boolean? gatk_ug_keep_realign_bam
		String? gatk_ug_defaultbasequalities
		Int freebayes_C = 1
		Int? freebayes_g
		Int freebayes_p = 2
		String? pileupcaller_bedfile
		String? pileupcaller_snpfile
		String pileupcaller_method = "randomHaploid"
		String pileupcaller_transitions_mode = "AllSites"
		Int pileupcaller_min_map_quality = 30
		Int pileupcaller_min_base_quality = 30
		String angsd_glmodel = "samtools"
		String angsd_glformat = "binary"
		Boolean? angsd_createfasta
		String angsd_fastamethod = "random"
		Boolean run_bcftools_stats = true
		Boolean? run_vcf2genome
		String? vcf2genome_outfile
		String? vcf2genome_header
		Int vcf2genome_minc = 5
		Int vcf2genome_minq = 30
		Float vcf2genome_minfreq = 0.8
		Boolean? run_multivcfanalyzer
		Boolean? write_allele_frequencies
		Int min_genotype_quality = 30
		Int min_base_coverage = 5
		Float min_allele_freq_hom = 0.9
		Float min_allele_freq_het = 0.9
		String? additional_vcf_files
		String reference_gff_annotations = "NA"
		String reference_gff_exclude = "NA"
		String snp_eff_results = "NA"
		Boolean? run_mtnucratio
		String mtnucratio_header = "MT"
		Boolean? run_sexdeterrmine
		String? sexdeterrmine_bedfile
		Boolean? run_nuclear_contamination
		String contamination_chrom_name = "X"
		Boolean? metagenomic_complexity_filter
		Float metagenomic_complexity_entropy = 0.3
		Boolean? run_metagenomic_screening
		String? metagenomic_tool
		String? database
		Int metagenomic_min_support_reads = 1
		Int percent_identity = 85
		String malt_mode = "BlastN"
		String malt_alignment_mode = "SemiGlobal"
		Int malt_top_percent = 1
		String malt_min_support_mode = "percent"
		Float malt_min_support_percent = 0.01
		Int malt_max_queries = 100
		String malt_memory_mode = "load"
		Boolean? malt_sam_output
		Boolean? run_maltextract
		String? maltextract_taxon_list
		String? maltextract_ncbifiles
		String maltextract_filter = "def_anc"
		Float maltextract_toppercent = 0.01
		Boolean? maltextract_destackingoff
		Boolean? maltextract_downsamplingoff
		Boolean? maltextract_duplicateremovaloff
		Boolean? maltextract_matches
		Boolean? maltextract_megansummary
		Float maltextract_percentidentity = 85.0
		Boolean? maltextract_topalignment

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			samplesheet = samplesheet,
			udg_type = udg_type,
			single_stranded = single_stranded,
			single_end = single_end,
			colour_chemistry = colour_chemistry,
			bam = bam,
			snpcapture_bed = snpcapture_bed,
			run_convertinputbam = run_convertinputbam,
			fasta = fasta,
			genome = genome,
			igenomes_base = igenomes_base,
			igenomes_ignore = igenomes_ignore,
			bwa_index = bwa_index,
			bt2_index = bt2_index,
			fasta_index = fasta_index,
			seq_dict = seq_dict,
			large_ref = large_ref,
			save_reference = save_reference,
			outdir = outdir,
			publish_dir_mode = publish_dir_mode,
			help = help,
			validate_params = validate_params,
			email = email,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			show_hidden_params = show_hidden_params,
			enable_conda = enable_conda,
			schema_ignore_params = schema_ignore_params,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			hostnames = hostnames,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			awsqueue = awsqueue,
			awsregion = awsregion,
			awscli = awscli,
			skip_fastqc = skip_fastqc,
			skip_adapterremoval = skip_adapterremoval,
			skip_preseq = skip_preseq,
			skip_deduplication = skip_deduplication,
			skip_damage_calculation = skip_damage_calculation,
			skip_qualimap = skip_qualimap,
			complexity_filter_poly_g = complexity_filter_poly_g,
			complexity_filter_poly_g_min = complexity_filter_poly_g_min,
			clip_forward_adaptor = clip_forward_adaptor,
			clip_reverse_adaptor = clip_reverse_adaptor,
			clip_adapters_list = clip_adapters_list,
			clip_readlength = clip_readlength,
			clip_min_read_quality = clip_min_read_quality,
			min_adap_overlap = min_adap_overlap,
			skip_collapse = skip_collapse,
			skip_trim = skip_trim,
			preserve5p = preserve5p,
			mergedonly = mergedonly,
			qualitymax = qualitymax,
			run_post_ar_trimming = run_post_ar_trimming,
			post_ar_trim_front = post_ar_trim_front,
			post_ar_trim_tail = post_ar_trim_tail,
			post_ar_trim_front2 = post_ar_trim_front2,
			post_ar_trim_tail2 = post_ar_trim_tail2,
			mapper = mapper,
			bwaalnn = bwaalnn,
			bwaalnk = bwaalnk,
			bwaalnl = bwaalnl,
			bwaalno = bwaalno,
			circularextension = circularextension,
			circulartarget = circulartarget,
			circularfilter = circularfilter,
			bt2_alignmode = bt2_alignmode,
			bt2_sensitivity = bt2_sensitivity,
			bt2n = bt2n,
			bt2l = bt2l,
			bt2_trim5 = bt2_trim5,
			bt2_trim3 = bt2_trim3,
			bt2_maxins = bt2_maxins,
			hostremoval_input_fastq = hostremoval_input_fastq,
			hostremoval_mode = hostremoval_mode,
			run_bam_filtering = run_bam_filtering,
			bam_mapping_quality_threshold = bam_mapping_quality_threshold,
			bam_filter_minreadlength = bam_filter_minreadlength,
			bam_unmapped_type = bam_unmapped_type,
			dedupper = dedupper,
			dedup_all_merged = dedup_all_merged,
			preseq_mode = preseq_mode,
			preseq_step_size = preseq_step_size,
			preseq_maxextrap = preseq_maxextrap,
			preseq_terms = preseq_terms,
			preseq_bootstrap = preseq_bootstrap,
			preseq_cval = preseq_cval,
			damageprofiler_length = damageprofiler_length,
			damageprofiler_threshold = damageprofiler_threshold,
			damageprofiler_yaxis = damageprofiler_yaxis,
			run_pmdtools = run_pmdtools,
			pmdtools_range = pmdtools_range,
			pmdtools_threshold = pmdtools_threshold,
			pmdtools_reference_mask = pmdtools_reference_mask,
			pmdtools_max_reads = pmdtools_max_reads,
			pmdtools_platypus = pmdtools_platypus,
			run_mapdamage_rescaling = run_mapdamage_rescaling,
			rescale_length_5p = rescale_length_5p,
			rescale_length_3p = rescale_length_3p,
			run_bedtools_coverage = run_bedtools_coverage,
			anno_file = anno_file,
			run_trim_bam = run_trim_bam,
			bamutils_clip_double_stranded_half_udg_left = bamutils_clip_double_stranded_half_udg_left,
			bamutils_clip_double_stranded_half_udg_right = bamutils_clip_double_stranded_half_udg_right,
			bamutils_clip_double_stranded_none_udg_left = bamutils_clip_double_stranded_none_udg_left,
			bamutils_clip_double_stranded_none_udg_right = bamutils_clip_double_stranded_none_udg_right,
			bamutils_clip_single_stranded_half_udg_left = bamutils_clip_single_stranded_half_udg_left,
			bamutils_clip_single_stranded_half_udg_right = bamutils_clip_single_stranded_half_udg_right,
			bamutils_clip_single_stranded_none_udg_left = bamutils_clip_single_stranded_none_udg_left,
			bamutils_clip_single_stranded_none_udg_right = bamutils_clip_single_stranded_none_udg_right,
			bamutils_softclip = bamutils_softclip,
			run_genotyping = run_genotyping,
			genotyping_tool = genotyping_tool,
			genotyping_source = genotyping_source,
			gatk_call_conf = gatk_call_conf,
			gatk_ploidy = gatk_ploidy,
			gatk_downsample = gatk_downsample,
			gatk_dbsnp = gatk_dbsnp,
			gatk_hc_out_mode = gatk_hc_out_mode,
			gatk_hc_emitrefconf = gatk_hc_emitrefconf,
			gatk_ug_out_mode = gatk_ug_out_mode,
			gatk_ug_genotype_model = gatk_ug_genotype_model,
			gatk_ug_keep_realign_bam = gatk_ug_keep_realign_bam,
			gatk_ug_defaultbasequalities = gatk_ug_defaultbasequalities,
			freebayes_C = freebayes_C,
			freebayes_g = freebayes_g,
			freebayes_p = freebayes_p,
			pileupcaller_bedfile = pileupcaller_bedfile,
			pileupcaller_snpfile = pileupcaller_snpfile,
			pileupcaller_method = pileupcaller_method,
			pileupcaller_transitions_mode = pileupcaller_transitions_mode,
			pileupcaller_min_map_quality = pileupcaller_min_map_quality,
			pileupcaller_min_base_quality = pileupcaller_min_base_quality,
			angsd_glmodel = angsd_glmodel,
			angsd_glformat = angsd_glformat,
			angsd_createfasta = angsd_createfasta,
			angsd_fastamethod = angsd_fastamethod,
			run_bcftools_stats = run_bcftools_stats,
			run_vcf2genome = run_vcf2genome,
			vcf2genome_outfile = vcf2genome_outfile,
			vcf2genome_header = vcf2genome_header,
			vcf2genome_minc = vcf2genome_minc,
			vcf2genome_minq = vcf2genome_minq,
			vcf2genome_minfreq = vcf2genome_minfreq,
			run_multivcfanalyzer = run_multivcfanalyzer,
			write_allele_frequencies = write_allele_frequencies,
			min_genotype_quality = min_genotype_quality,
			min_base_coverage = min_base_coverage,
			min_allele_freq_hom = min_allele_freq_hom,
			min_allele_freq_het = min_allele_freq_het,
			additional_vcf_files = additional_vcf_files,
			reference_gff_annotations = reference_gff_annotations,
			reference_gff_exclude = reference_gff_exclude,
			snp_eff_results = snp_eff_results,
			run_mtnucratio = run_mtnucratio,
			mtnucratio_header = mtnucratio_header,
			run_sexdeterrmine = run_sexdeterrmine,
			sexdeterrmine_bedfile = sexdeterrmine_bedfile,
			run_nuclear_contamination = run_nuclear_contamination,
			contamination_chrom_name = contamination_chrom_name,
			metagenomic_complexity_filter = metagenomic_complexity_filter,
			metagenomic_complexity_entropy = metagenomic_complexity_entropy,
			run_metagenomic_screening = run_metagenomic_screening,
			metagenomic_tool = metagenomic_tool,
			database = database,
			metagenomic_min_support_reads = metagenomic_min_support_reads,
			percent_identity = percent_identity,
			malt_mode = malt_mode,
			malt_alignment_mode = malt_alignment_mode,
			malt_top_percent = malt_top_percent,
			malt_min_support_mode = malt_min_support_mode,
			malt_min_support_percent = malt_min_support_percent,
			malt_max_queries = malt_max_queries,
			malt_memory_mode = malt_memory_mode,
			malt_sam_output = malt_sam_output,
			run_maltextract = run_maltextract,
			maltextract_taxon_list = maltextract_taxon_list,
			maltextract_ncbifiles = maltextract_ncbifiles,
			maltextract_filter = maltextract_filter,
			maltextract_toppercent = maltextract_toppercent,
			maltextract_destackingoff = maltextract_destackingoff,
			maltextract_downsamplingoff = maltextract_downsamplingoff,
			maltextract_duplicateremovaloff = maltextract_duplicateremovaloff,
			maltextract_matches = maltextract_matches,
			maltextract_megansummary = maltextract_megansummary,
			maltextract_percentidentity = maltextract_percentidentity,
			maltextract_topalignment = maltextract_topalignment,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-eager/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		File samplesheet
		String udg_type = "none"
		Boolean? single_stranded
		Boolean? single_end
		Int colour_chemistry = 4
		Boolean? bam
		String? snpcapture_bed
		Boolean? run_convertinputbam
		File? fasta
		String? genome
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		String? bwa_index
		String? bt2_index
		String? fasta_index
		String? seq_dict
		Boolean? large_ref
		Boolean? save_reference
		String outdir = "./results"
		String publish_dir_mode = "copy"
		Boolean? help
		Boolean validate_params = true
		String? email
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean? show_hidden_params
		Boolean? enable_conda
		String schema_ignore_params = "genomes"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		String? awsqueue
		String awsregion = "eu-west-1"
		String? awscli
		Boolean? skip_fastqc
		Boolean? skip_adapterremoval
		Boolean? skip_preseq
		Boolean? skip_deduplication
		Boolean? skip_damage_calculation
		Boolean? skip_qualimap
		Boolean? complexity_filter_poly_g
		Int complexity_filter_poly_g_min = 10
		String clip_forward_adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
		String clip_reverse_adaptor = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
		String? clip_adapters_list
		Int clip_readlength = 30
		Int clip_min_read_quality = 20
		Int min_adap_overlap = 1
		Boolean? skip_collapse
		Boolean? skip_trim
		Boolean? preserve5p
		Boolean? mergedonly
		Int qualitymax = 41
		Boolean? run_post_ar_trimming
		Int post_ar_trim_front = 7
		Int post_ar_trim_tail = 7
		Int post_ar_trim_front2 = 7
		Int post_ar_trim_tail2 = 7
		String mapper = "bwaaln"
		Float bwaalnn = 0.01
		Int bwaalnk = 2
		Int bwaalnl = 1024
		Int bwaalno = 2
		Int circularextension = 500
		String circulartarget = "MT"
		Boolean? circularfilter
		String bt2_alignmode = "local"
		String bt2_sensitivity = "sensitive"
		Int? bt2n
		Int? bt2l
		Int? bt2_trim5
		Int? bt2_trim3
		Int bt2_maxins = 500
		Boolean? hostremoval_input_fastq
		String hostremoval_mode = "remove"
		Boolean? run_bam_filtering
		Int? bam_mapping_quality_threshold
		Int? bam_filter_minreadlength
		String bam_unmapped_type = "discard"
		String dedupper = "markduplicates"
		Boolean? dedup_all_merged
		String preseq_mode = "c_curve"
		Int preseq_step_size = 1000
		Int preseq_maxextrap = 10000000000
		Int preseq_terms = 100
		Int preseq_bootstrap = 100
		Float preseq_cval = 0.95
		Int damageprofiler_length = 100
		Int damageprofiler_threshold = 15
		Float damageprofiler_yaxis = 0.3
		Boolean? run_pmdtools
		Int pmdtools_range = 10
		Int pmdtools_threshold = 3
		String? pmdtools_reference_mask
		Int pmdtools_max_reads = 10000
		Boolean? pmdtools_platypus
		Boolean? run_mapdamage_rescaling
		Int rescale_length_5p = 12
		Int rescale_length_3p = 12
		Boolean? run_bedtools_coverage
		String? anno_file
		Boolean? run_trim_bam
		Int? bamutils_clip_double_stranded_half_udg_left
		Int? bamutils_clip_double_stranded_half_udg_right
		Int? bamutils_clip_double_stranded_none_udg_left
		Int? bamutils_clip_double_stranded_none_udg_right
		Int? bamutils_clip_single_stranded_half_udg_left
		Int? bamutils_clip_single_stranded_half_udg_right
		Int? bamutils_clip_single_stranded_none_udg_left
		Int? bamutils_clip_single_stranded_none_udg_right
		Boolean? bamutils_softclip
		Boolean? run_genotyping
		String? genotyping_tool
		String genotyping_source = "raw"
		Int gatk_call_conf = 30
		Int gatk_ploidy = 2
		Int gatk_downsample = 250
		String? gatk_dbsnp
		String gatk_hc_out_mode = "EMIT_VARIANTS_ONLY"
		String gatk_hc_emitrefconf = "GVCF"
		String gatk_ug_out_mode = "EMIT_VARIANTS_ONLY"
		String gatk_ug_genotype_model = "SNP"
		Boolean? gatk_ug_keep_realign_bam
		String? gatk_ug_defaultbasequalities
		Int freebayes_C = 1
		Int? freebayes_g
		Int freebayes_p = 2
		String? pileupcaller_bedfile
		String? pileupcaller_snpfile
		String pileupcaller_method = "randomHaploid"
		String pileupcaller_transitions_mode = "AllSites"
		Int pileupcaller_min_map_quality = 30
		Int pileupcaller_min_base_quality = 30
		String angsd_glmodel = "samtools"
		String angsd_glformat = "binary"
		Boolean? angsd_createfasta
		String angsd_fastamethod = "random"
		Boolean run_bcftools_stats = true
		Boolean? run_vcf2genome
		String? vcf2genome_outfile
		String? vcf2genome_header
		Int vcf2genome_minc = 5
		Int vcf2genome_minq = 30
		Float vcf2genome_minfreq = 0.8
		Boolean? run_multivcfanalyzer
		Boolean? write_allele_frequencies
		Int min_genotype_quality = 30
		Int min_base_coverage = 5
		Float min_allele_freq_hom = 0.9
		Float min_allele_freq_het = 0.9
		String? additional_vcf_files
		String reference_gff_annotations = "NA"
		String reference_gff_exclude = "NA"
		String snp_eff_results = "NA"
		Boolean? run_mtnucratio
		String mtnucratio_header = "MT"
		Boolean? run_sexdeterrmine
		String? sexdeterrmine_bedfile
		Boolean? run_nuclear_contamination
		String contamination_chrom_name = "X"
		Boolean? metagenomic_complexity_filter
		Float metagenomic_complexity_entropy = 0.3
		Boolean? run_metagenomic_screening
		String? metagenomic_tool
		String? database
		Int metagenomic_min_support_reads = 1
		Int percent_identity = 85
		String malt_mode = "BlastN"
		String malt_alignment_mode = "SemiGlobal"
		Int malt_top_percent = 1
		String malt_min_support_mode = "percent"
		Float malt_min_support_percent = 0.01
		Int malt_max_queries = 100
		String malt_memory_mode = "load"
		Boolean? malt_sam_output
		Boolean? run_maltextract
		String? maltextract_taxon_list
		String? maltextract_ncbifiles
		String maltextract_filter = "def_anc"
		Float maltextract_toppercent = 0.01
		Boolean? maltextract_destackingoff
		Boolean? maltextract_downsamplingoff
		Boolean? maltextract_duplicateremovaloff
		Boolean? maltextract_matches
		Boolean? maltextract_megansummary
		Float maltextract_percentidentity = 85.0
		Boolean? maltextract_topalignment

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /eager-2.4.2  -profile truwl,nfcore-eager  --input ~{samplesheet} 	~{"--samplesheet '" + samplesheet + "'"}	~{"--udg_type '" + udg_type + "'"}	~{true="--single_stranded  " false="" single_stranded}	~{true="--single_end  " false="" single_end}	~{"--colour_chemistry " + colour_chemistry}	~{true="--bam  " false="" bam}	~{"--snpcapture_bed '" + snpcapture_bed + "'"}	~{true="--run_convertinputbam  " false="" run_convertinputbam}	~{"--fasta '" + fasta + "'"}	~{"--genome '" + genome + "'"}	~{"--igenomes_base '" + igenomes_base + "'"}	~{true="--igenomes_ignore  " false="" igenomes_ignore}	~{"--bwa_index '" + bwa_index + "'"}	~{"--bt2_index '" + bt2_index + "'"}	~{"--fasta_index '" + fasta_index + "'"}	~{"--seq_dict '" + seq_dict + "'"}	~{true="--large_ref  " false="" large_ref}	~{true="--save_reference  " false="" save_reference}	~{"--outdir '" + outdir + "'"}	~{"--publish_dir_mode '" + publish_dir_mode + "'"}	~{true="--help  " false="" help}	~{true="--validate_params  " false="" validate_params}	~{"--email '" + email + "'"}	~{"--email_on_fail '" + email_on_fail + "'"}	~{true="--plaintext_email  " false="" plaintext_email}	~{"--max_multiqc_email_size '" + max_multiqc_email_size + "'"}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--multiqc_config '" + multiqc_config + "'"}	~{"--tracedir '" + tracedir + "'"}	~{true="--show_hidden_params  " false="" show_hidden_params}	~{true="--enable_conda  " false="" enable_conda}	~{"--schema_ignore_params '" + schema_ignore_params + "'"}	~{"--max_cpus " + max_cpus}	~{"--max_memory '" + max_memory + "'"}	~{"--max_time '" + max_time + "'"}	~{"--custom_config_version '" + custom_config_version + "'"}	~{"--custom_config_base '" + custom_config_base + "'"}	~{"--hostnames '" + hostnames + "'"}	~{"--config_profile_name '" + config_profile_name + "'"}	~{"--config_profile_description '" + config_profile_description + "'"}	~{"--config_profile_contact '" + config_profile_contact + "'"}	~{"--config_profile_url '" + config_profile_url + "'"}	~{"--awsqueue '" + awsqueue + "'"}	~{"--awsregion '" + awsregion + "'"}	~{"--awscli '" + awscli + "'"}	~{true="--skip_fastqc  " false="" skip_fastqc}	~{true="--skip_adapterremoval  " false="" skip_adapterremoval}	~{true="--skip_preseq  " false="" skip_preseq}	~{true="--skip_deduplication  " false="" skip_deduplication}	~{true="--skip_damage_calculation  " false="" skip_damage_calculation}	~{true="--skip_qualimap  " false="" skip_qualimap}	~{true="--complexity_filter_poly_g  " false="" complexity_filter_poly_g}	~{"--complexity_filter_poly_g_min " + complexity_filter_poly_g_min}	~{"--clip_forward_adaptor '" + clip_forward_adaptor + "'"}	~{"--clip_reverse_adaptor '" + clip_reverse_adaptor + "'"}	~{"--clip_adapters_list '" + clip_adapters_list + "'"}	~{"--clip_readlength " + clip_readlength}	~{"--clip_min_read_quality " + clip_min_read_quality}	~{"--min_adap_overlap " + min_adap_overlap}	~{true="--skip_collapse  " false="" skip_collapse}	~{true="--skip_trim  " false="" skip_trim}	~{true="--preserve5p  " false="" preserve5p}	~{true="--mergedonly  " false="" mergedonly}	~{"--qualitymax " + qualitymax}	~{true="--run_post_ar_trimming  " false="" run_post_ar_trimming}	~{"--post_ar_trim_front " + post_ar_trim_front}	~{"--post_ar_trim_tail " + post_ar_trim_tail}	~{"--post_ar_trim_front2 " + post_ar_trim_front2}	~{"--post_ar_trim_tail2 " + post_ar_trim_tail2}	~{"--mapper '" + mapper + "'"}	~{"--bwaalnn " + bwaalnn}	~{"--bwaalnk " + bwaalnk}	~{"--bwaalnl " + bwaalnl}	~{"--bwaalno " + bwaalno}	~{"--circularextension " + circularextension}	~{"--circulartarget '" + circulartarget + "'"}	~{true="--circularfilter  " false="" circularfilter}	~{"--bt2_alignmode '" + bt2_alignmode + "'"}	~{"--bt2_sensitivity '" + bt2_sensitivity + "'"}	~{"--bt2n " + bt2n}	~{"--bt2l " + bt2l}	~{"--bt2_trim5 " + bt2_trim5}	~{"--bt2_trim3 " + bt2_trim3}	~{"--bt2_maxins " + bt2_maxins}	~{true="--hostremoval_input_fastq  " false="" hostremoval_input_fastq}	~{"--hostremoval_mode '" + hostremoval_mode + "'"}	~{true="--run_bam_filtering  " false="" run_bam_filtering}	~{"--bam_mapping_quality_threshold " + bam_mapping_quality_threshold}	~{"--bam_filter_minreadlength " + bam_filter_minreadlength}	~{"--bam_unmapped_type '" + bam_unmapped_type + "'"}	~{"--dedupper '" + dedupper + "'"}	~{true="--dedup_all_merged  " false="" dedup_all_merged}	~{"--preseq_mode '" + preseq_mode + "'"}	~{"--preseq_step_size " + preseq_step_size}	~{"--preseq_maxextrap " + preseq_maxextrap}	~{"--preseq_terms " + preseq_terms}	~{"--preseq_bootstrap " + preseq_bootstrap}	~{"--preseq_cval " + preseq_cval}	~{"--damageprofiler_length " + damageprofiler_length}	~{"--damageprofiler_threshold " + damageprofiler_threshold}	~{"--damageprofiler_yaxis " + damageprofiler_yaxis}	~{true="--run_pmdtools  " false="" run_pmdtools}	~{"--pmdtools_range " + pmdtools_range}	~{"--pmdtools_threshold " + pmdtools_threshold}	~{"--pmdtools_reference_mask '" + pmdtools_reference_mask + "'"}	~{"--pmdtools_max_reads " + pmdtools_max_reads}	~{true="--pmdtools_platypus  " false="" pmdtools_platypus}	~{true="--run_mapdamage_rescaling  " false="" run_mapdamage_rescaling}	~{"--rescale_length_5p " + rescale_length_5p}	~{"--rescale_length_3p " + rescale_length_3p}	~{true="--run_bedtools_coverage  " false="" run_bedtools_coverage}	~{"--anno_file '" + anno_file + "'"}	~{true="--run_trim_bam  " false="" run_trim_bam}	~{"--bamutils_clip_double_stranded_half_udg_left " + bamutils_clip_double_stranded_half_udg_left}	~{"--bamutils_clip_double_stranded_half_udg_right " + bamutils_clip_double_stranded_half_udg_right}	~{"--bamutils_clip_double_stranded_none_udg_left " + bamutils_clip_double_stranded_none_udg_left}	~{"--bamutils_clip_double_stranded_none_udg_right " + bamutils_clip_double_stranded_none_udg_right}	~{"--bamutils_clip_single_stranded_half_udg_left " + bamutils_clip_single_stranded_half_udg_left}	~{"--bamutils_clip_single_stranded_half_udg_right " + bamutils_clip_single_stranded_half_udg_right}	~{"--bamutils_clip_single_stranded_none_udg_left " + bamutils_clip_single_stranded_none_udg_left}	~{"--bamutils_clip_single_stranded_none_udg_right " + bamutils_clip_single_stranded_none_udg_right}	~{true="--bamutils_softclip  " false="" bamutils_softclip}	~{true="--run_genotyping  " false="" run_genotyping}	~{"--genotyping_tool '" + genotyping_tool + "'"}	~{"--genotyping_source '" + genotyping_source + "'"}	~{"--gatk_call_conf " + gatk_call_conf}	~{"--gatk_ploidy " + gatk_ploidy}	~{"--gatk_downsample " + gatk_downsample}	~{"--gatk_dbsnp '" + gatk_dbsnp + "'"}	~{"--gatk_hc_out_mode '" + gatk_hc_out_mode + "'"}	~{"--gatk_hc_emitrefconf '" + gatk_hc_emitrefconf + "'"}	~{"--gatk_ug_out_mode '" + gatk_ug_out_mode + "'"}	~{"--gatk_ug_genotype_model '" + gatk_ug_genotype_model + "'"}	~{true="--gatk_ug_keep_realign_bam  " false="" gatk_ug_keep_realign_bam}	~{"--gatk_ug_defaultbasequalities '" + gatk_ug_defaultbasequalities + "'"}	~{"--freebayes_C " + freebayes_C}	~{"--freebayes_g " + freebayes_g}	~{"--freebayes_p " + freebayes_p}	~{"--pileupcaller_bedfile '" + pileupcaller_bedfile + "'"}	~{"--pileupcaller_snpfile '" + pileupcaller_snpfile + "'"}	~{"--pileupcaller_method '" + pileupcaller_method + "'"}	~{"--pileupcaller_transitions_mode '" + pileupcaller_transitions_mode + "'"}	~{"--pileupcaller_min_map_quality " + pileupcaller_min_map_quality}	~{"--pileupcaller_min_base_quality " + pileupcaller_min_base_quality}	~{"--angsd_glmodel '" + angsd_glmodel + "'"}	~{"--angsd_glformat '" + angsd_glformat + "'"}	~{true="--angsd_createfasta  " false="" angsd_createfasta}	~{"--angsd_fastamethod '" + angsd_fastamethod + "'"}	~{true="--run_bcftools_stats  " false="" run_bcftools_stats}	~{true="--run_vcf2genome  " false="" run_vcf2genome}	~{"--vcf2genome_outfile '" + vcf2genome_outfile + "'"}	~{"--vcf2genome_header '" + vcf2genome_header + "'"}	~{"--vcf2genome_minc " + vcf2genome_minc}	~{"--vcf2genome_minq " + vcf2genome_minq}	~{"--vcf2genome_minfreq " + vcf2genome_minfreq}	~{true="--run_multivcfanalyzer  " false="" run_multivcfanalyzer}	~{true="--write_allele_frequencies  " false="" write_allele_frequencies}	~{"--min_genotype_quality " + min_genotype_quality}	~{"--min_base_coverage " + min_base_coverage}	~{"--min_allele_freq_hom " + min_allele_freq_hom}	~{"--min_allele_freq_het " + min_allele_freq_het}	~{"--additional_vcf_files '" + additional_vcf_files + "'"}	~{"--reference_gff_annotations '" + reference_gff_annotations + "'"}	~{"--reference_gff_exclude '" + reference_gff_exclude + "'"}	~{"--snp_eff_results '" + snp_eff_results + "'"}	~{true="--run_mtnucratio  " false="" run_mtnucratio}	~{"--mtnucratio_header '" + mtnucratio_header + "'"}	~{true="--run_sexdeterrmine  " false="" run_sexdeterrmine}	~{"--sexdeterrmine_bedfile '" + sexdeterrmine_bedfile + "'"}	~{true="--run_nuclear_contamination  " false="" run_nuclear_contamination}	~{"--contamination_chrom_name '" + contamination_chrom_name + "'"}	~{true="--metagenomic_complexity_filter  " false="" metagenomic_complexity_filter}	~{"--metagenomic_complexity_entropy " + metagenomic_complexity_entropy}	~{true="--run_metagenomic_screening  " false="" run_metagenomic_screening}	~{"--metagenomic_tool '" + metagenomic_tool + "'"}	~{"--database '" + database + "'"}	~{"--metagenomic_min_support_reads " + metagenomic_min_support_reads}	~{"--percent_identity " + percent_identity}	~{"--malt_mode '" + malt_mode + "'"}	~{"--malt_alignment_mode '" + malt_alignment_mode + "'"}	~{"--malt_top_percent " + malt_top_percent}	~{"--malt_min_support_mode '" + malt_min_support_mode + "'"}	~{"--malt_min_support_percent " + malt_min_support_percent}	~{"--malt_max_queries " + malt_max_queries}	~{"--malt_memory_mode '" + malt_memory_mode + "'"}	~{true="--malt_sam_output  " false="" malt_sam_output}	~{true="--run_maltextract  " false="" run_maltextract}	~{"--maltextract_taxon_list '" + maltextract_taxon_list + "'"}	~{"--maltextract_ncbifiles '" + maltextract_ncbifiles + "'"}	~{"--maltextract_filter '" + maltextract_filter + "'"}	~{"--maltextract_toppercent " + maltextract_toppercent}	~{true="--maltextract_destackingoff  " false="" maltextract_destackingoff}	~{true="--maltextract_downsamplingoff  " false="" maltextract_downsamplingoff}	~{true="--maltextract_duplicateremovaloff  " false="" maltextract_duplicateremovaloff}	~{true="--maltextract_matches  " false="" maltextract_matches}	~{true="--maltextract_megansummary  " false="" maltextract_megansummary}	~{"--maltextract_percentidentity " + maltextract_percentidentity}	~{true="--maltextract_topalignment  " false="" maltextract_topalignment}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*html")
    }
    runtime {
        docker: "truwl/nfcore-eager:2.4.2_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    