process ARGENTAG_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    container 'docker://biocorecrg/taggy_demux:1.1'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('*.split.fq')                         , emit: out
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def sample_name   = input_file.name.replaceAll(/\.fq$/, '')

    """

    split.sh \
        ${args} \
        -i ${input_file} \
        -o ${sample_name}.split.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_tags: 1.0

    END_VERSIONS

    """


}
