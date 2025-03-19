process computeSilhouette {
    tag 'computeSilhouette'

    publishDir "${launchDir}/${params.outdir}", mode: 'copy'

    input:
        val  dataset_versions_ready
        path dataset_json

    output:
        path "silhouette_scores_*.csv", emit: scores_csvs
	val  true, emit: scores_ready_ch

    script:
    """
    echo "Computing Silhouette Scores for datasets"

    python "${launchDir}/bin/compute_silhouette.py" "$dataset_json"
    """
}
