// sample pipeline //
// Installs the tools using conda (conda enabled: true in config file)
// Accepts reads with prefix *_(1/2).(fastq/fq).gz
// Overwrite outdir name in cli using --output parameter
//params.reads = "*_R{1,2}.{fq,fastq}.gz"
params.output = "./testingResults5"
//params.pwd = "./"

in_ch = Channel.fromFilePairs( "./*_R{1,2}.{fq,fastq}.gz") // creating a channel for reads
new_ch = in_ch.toList()
new_in = Channel.fromFilePairs("./*_trimmed_{1,2}.fastq.gz")
new_in2 = new_in.toList()

process FASTQC {
    tag "FASTQC on ${sampleID} reads"	
	publishDir "${params.output}", mode: 'copy'
	// Use the correct conda directive
	conda "bioconda::fastqc=0.12.1"
	debug true

	input:
	tuple val(sampleID), path(reads)

	output:
	path "*_fastqc.{zip,html}"

	script:
	"""
	echo ${reads}
	fastqc -q ${reads}
	"""
}
process FASTP {
	tag "${sampleID} Read correction using FASTP"
    publishDir "${params.output}", mode: 'copy'
    cpus 4
	debug true
    conda "bioconda::fastp"

	input:
	tuple val(sampleID), path(reads)
	
	output:
	tuple val(sampleID), path("${sampleID}_trimmed_{1,2}.fastq.gz")

	script:
	"""
	fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_trimmed_1.fastq.gz -O ${sampleID}_trimmed_2.fastq.gz
	"""

}

process FLASH {
    tag "Overlapping/stitching PE reads of ${sampleID} with FLASH"
    publishDir "${params.output}", mode: 'copy'
    cpus 4
	debug true
    conda "bioconda::flash2"

    
	input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.{extendedFrags,notCombined_1,notCombined_2}.fastq.gz")

    script:
    """
    echo $reads
    flash2 -o ${sampleID} -z -t $task.cpus -M 150 ${reads[0]} ${reads[1]} 
    """
}


workflow {
	 in_ch.view()
	 fastqc_ch = FASTQC(in_ch)
	 fastp_ch = FASTP(in_ch)
	 fastp_ch.view()
	 flash_ch = FLASH(fastp_ch)
	 flash_ch.view()
	}
