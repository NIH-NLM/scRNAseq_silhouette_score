params.datadir              = "data"
params.outdir               = "results"
params.test_mode            = false
params.collection_info      = "collections_info.json"
params.final_report_pdf     = "final_report.pdf"
params.final_report_html    = "final_report.html"
params.collections_ready_ch = false
params.datasets_ready_ch    = false
params.versions_ready_ch    = false
params.scores_ready_ch      = false
params.plots_ready_ch       = false

// include all the modules
include { computeSilhouette } from './modules/computeSilhouette.nf'
include { fetchCollections }  from './modules/fetchCollections.nf'
include { fetchDatasets }     from './modules/fetchDatasets.nf'
include { generatePlots }     from './modules/generatePlots.nf'
include { mergeResults }      from './modules/mergeResults.nf'
include { splitCollections }  from './modules/splitCollections.nf'
include { splitDatasets }     from './modules/splitDatasets.nf'

// Define Workflow Execution Order

workflow {

    collection_json = fetchCollections(params.collection_info)
    collection_json.view()

    (collections_ch, collections_ready_ch) = splitCollections(collection_json)

    (datasets_ch, datasets_ready_ch) = fetchDatasets( collections_ready_ch.collect(),
                                                      collections_ch.flatten(),
						      params.test_mode )
						      
    ( datasets_versions_ch, versions_ready_ch ) = splitDatasets ( datasets_ready_ch.collect(),
                                                                  datasets_ch.flatten() )
    
    ( scores_csv_ch, scores_ready_ch ) = computeSilhouette ( versions_ready_ch.collect(),
                                                             datasets_versions_ch.flatten() )
							     
    ( plots_ch, plots_ready_ch ) = generatePlots ( scores_ready_ch,
                                                   scores_csv_ch.flatten() )

    mergeResults ( plots_ready_ch.collect(),
                   plots_ch,
		   params.final_report_pdf,
		   params.final_report_html )
    

}
