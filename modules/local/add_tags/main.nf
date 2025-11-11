process ADD_TAGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://biocorecrg/pysam:0.3'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path('*.tagged.bam')                       , emit: out
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def sample_name   = bam_file.name.replaceAll(/\.bam$/, '')

    """

    Tags2bam.py \
        ${args} \
        -i ${bam_file} \
        -o ${sample_name}.tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_tags: 1.0

    END_VERSIONS

    """


}
