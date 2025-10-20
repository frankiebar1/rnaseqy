process MERGY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_matplotlib_pandas_scikit-learn:6145ce6739a97347':
        'community.wave.seqera.io/library/pip_matplotlib_pandas_scikit-learn:17dbfac1cd918c7b' }"

    input:
    
    path(indir)
    path(outdir)
    path(placeholder)

    output:
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    count_plotting.py --indir ${indir} --outdir ${outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mergy: 1.0.0
    END_VERSIONS
    """
}
