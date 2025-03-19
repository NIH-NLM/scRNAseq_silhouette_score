process splitDatasets {
    tag 'split_datasets'
    
    publishDir "${launchDir}/${params.datadir}", mode: 'copy'
    
    input:
        path dataset_json each
        
    output:
        path "dataset_*.json", emit: dataset_jsons

    script:
    """
    bash "${launchDir}/bin/splitDatasets.sh"
    """
}

