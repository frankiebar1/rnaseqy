process UNZIPPER {
    tag "${meta.id}"
    label 'process_low'
    publishDir "intermediate/unzipped", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/nf-core-base:2.2' :
    'community.wave.seqera.io/library/cutadapt_trim-galore_pigz:a98edd405b34582d' }"


    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.fq"), emit: reads
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        for f in ${reads}
        do
            base=\$(basename \$f .gz)
            gunzip -c \$f > \$base
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pigz: \$(pigz --version 2>&1 | head -n 1)
        END_VERSIONS
        """
}