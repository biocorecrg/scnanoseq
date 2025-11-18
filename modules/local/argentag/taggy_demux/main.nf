process ARGENTAG_TAGGY_DEMUX {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker://biocorecrg/taggy_demux:1.1'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('*.demux.fastq')                            , emit: out
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def sample_name   = input_file.name.replaceAll(/\.fastq$/, '')

    """
    cp -r /opt/scripts/res .

    taggy_demux \
        --in-fmt="fastq" \
        ${args} \
        --output-dir . \
        ${input_file} \

    cat *.fastq > ${sample_name}.demux.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_tags: 1.0

    END_VERSIONS

    """


}
