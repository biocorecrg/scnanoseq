process SPLITPIPE_PRE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'biocorecrg/spipe:1.6.1'

    input:
    tuple val(meta),   path(reads)
    path(index)
    path(parfile_parse)
    
    output:
    tuple val(meta_updated), path("${meta.id}/process/barcode_head.fastq.gz")               , emit: out
    path  "versions.yml"                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def reads1        = []
    def reads2        = []
    meta.single_end ? [reads].flatten().each { reads1 << it } : reads.eachWithIndex { v, ix -> (ix & 1 ? reads2 : reads1) << v }
    def input_reads   = meta.single_end ? "-r ${reads1.join(" ")}" : "--fq1 ${reads1.join(" ")} --fq2 ${reads2.join(" ")}"
    def parfile      = parfile_parse ? "--parfile ${parfile_parse}" : ""
    meta_updated = meta + [single_end: true]

    """
    split-pipe \
       --mode pre --one_step ${args}  \
       ${input_reads} \
       --genome_dir ${index} \
       --output_dir ${prefix} \
       ${parfile} \
       --nthreads ${task.cpus}
  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split-pipe: \$(split-pipe --version | cut -d ' ' -f 2 | sed 's/v//g')
END_VERSIONS
    """


}