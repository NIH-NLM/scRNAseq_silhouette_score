process splitCollections {
    tag 'split_collections'

    publishDir "${launchDir}/${params.datadir}", mode: 'copy'

    input:
        path collection_json
        
    output:
        path "collection_*.json", emit: collections_json
        val true, emit: collections_ready_ch

    script:
    """
    bash "${launchDir}/bin/splitCollections.sh" "$collection_json"    
    """
}

