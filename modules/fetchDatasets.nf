process fetchDatasets {
    tag 'fetch_datasets'

    publishDir "${launchDir}/${params.datadir}", mode: 'copy'

    input:
        path collection_json
	val test_mode

    output:
        path "dataset_*.json", emit: dataset_jsons

    script:
    """
    echo "Fetching datasets from collection file and writing out a datasets.json file for each dataset_id (Test Mode: $test_mode)"
    
    python "${launchDir}/bin/fetch_datasets.py" "$collection_json" "$test_mode"
    """
}
