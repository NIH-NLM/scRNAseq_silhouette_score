process splitDatasets {
    tag 'split_datasets'
    
    publishDir "${launchDir}/${params.datadir}", mode: 'copy'
    
    input:
        val  datasets_ready
        path dataset_json
        
    output:
        path "dataset_*.json", emit: dataset_jsons
        val  true, emit: versions_ready_ch
	
    script:
    """
    bash "${launchDir}/bin/splitDatasets.sh" "$dataset_json"
    """
}

