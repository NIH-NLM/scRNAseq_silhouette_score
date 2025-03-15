process fetchDatasets {
    tag 'fetch_datasets'

    publishDir "${launchDir}/${params.datadir}", mode: 'copy'

    input:
        path collections_json
	val test_mode

    output:
        path "dataset*.json", emit: dataset_jsons

    script:
    """
    echo "Fetching datasets from collection file and writing out a datasets.json file for each dataset_id (Test Mode: $test_mode)"
    
    python "${launchDir}/bin/fetch_datasets.py" "$collections_json" "$test_mode"
    """
}
