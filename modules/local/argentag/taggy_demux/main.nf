process ARGENTAG_TAGGY_DEMUX {
    tag "${meta.id}"
    label 'process_high'
    label 'process_long'
    label 'process_high_cpus'

    container 'docker://biocorecrg/taggy_demux:1.1'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path('*.tagged.fastq')                     , emit: out
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def sample_name   = input_file.name.replaceAll(/\.fq$/, '')

    """
    cp -r /opt/scripts/res .

    taggy_demux \
        --in-fmt="fastq" \
        ${args} \
        --output-dir . \
        ${input_file} \
        -T ${task.cpus}

    cat *.fastq > ${sample_name}.tagged.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_tags: 1.0

    END_VERSIONS

    """


}
