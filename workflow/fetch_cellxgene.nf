process fetchCellxgene {
    input:
        val test_mode_flag

    output:
        path "collections_info.json"

    script:
    """
    python "${launchDir}/bin/fetch_cellxgene.py" ${test_mode_flag}
    """

    // ✅ Nextflow will copy collections_info.json to ${launchDir}/results/
    publishDir "${launchDir}/results/", mode: 'copy'
}


