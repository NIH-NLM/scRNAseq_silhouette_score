params.datadir           = "data"
params.outdir            = "results"
params.test_mode         = false
params.collection_info   = "collections_info.json"
params.final_report_pdf  = "final_report.pdf"
params.final_report_html = "final_report.html"

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

    // **Step 1: Fetch collections (single file)**
    collection_json = fetchCollections(params.collection_info)

    // **Step 2: Split collections into multiple JSONs**
    collections_json = channel.of(splitCollections(collection_json))

    // **Step 3: Fetch dataset metadata for each collection**
    datasets_json = channel.of(fetchDatasets(collections_json,params.test_mode))
  
    // **Step 4: Split datasets into individual dataset versions**
    dataset_versions_json = channel.of(splitDatasets(datasets_json))

    // **Step 5: Compute silhouette scores for each dataset version**
    scores_csvs = channel.of(computeSilhouette(dataset_versions_json))

    // **Step 6: Generate plots**
    plots = channel.of(generatePlots(scores_csvs))

    // **Step 7: Merge results into a final report**
    mergeResults(plots, params.final_report_pdf, params.final_report_html)
}
