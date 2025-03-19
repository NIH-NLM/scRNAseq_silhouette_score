process generatePlots {
    tag 'generatePlots'

    publishDir "${launchDir}/${params.outdir}", mode: 'copy'

    input:
        val  scores_ready
        path scores_csv

    output:
        path "final_plots_*.png", emit: final_plots
        val  true, emit: plots_ready_ch

    script:
    """
    echo "Generating plots from silhouette scores"

    python "${launchDir}/bin/generate_plots.py" "$scores_csv"
    """
}
