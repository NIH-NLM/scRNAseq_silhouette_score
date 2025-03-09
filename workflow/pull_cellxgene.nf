process fetchDatasets {
    input:
        val test_mode_flag

    output:
        path "datasets_info.json"

    script:
    """
    python "${launchDir}/bin/fetch_cellxgene.py" "${test_mode_flag}"
    """
}
