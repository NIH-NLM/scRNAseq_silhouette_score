process fetchDatasets {
    tag 'fetch_datasets'

    publishDir "${launchDir}/${params.datadir}", mode: 'copy'

    input:
        val  collection_ready
        path collection_ch
	val  test_mode

    output:
        path "datasets_*.json", emit: datasets_jsons
	val  true, emit: datasets_ready_ch

    script:
    """
    echo "Fetching datasets from collection file and writing out a datasets.json file for each dataset_id (Test Mode: $test_mode)"
    
    python "${launchDir}/bin/fetch_datasets.py" "$collection_ch" "$test_mode"
    """
}
