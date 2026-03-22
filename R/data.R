#' Sample data used in the FlexiDeconv vignette
#'
#' @format It contains a list having:
#'  1. Spatial Transcriptomics data generated from single-cell profiling of the
#'    hypothalamic preoptic region of mouse.
#'  2. Metadata for the Spatial Transcriptomics data (e.g., location).
#'  3. Relevant reference matrix.
#'  4. Ground truth gene expression profile for every cell type.
#'  5. Ground truth cell type proportion for every pixel.
#'
#' @source \url{https://datadryad.org/dataset/doi:10.5061/dryad.8t8s248}
#'
#' @usage data(mouseHypothalamus)
"mouseHypothalamus"



#' Sample data used in the inference vignette
#'
#' @format It contains a deconvolution result from applying imperfect reference
#' to the mouse hypothalamus data, in the standard FlexiDeconv output format
#'
#' @source \url{https://datadryad.org/dataset/doi:10.5061/dryad.8t8s248}
#'
#' @usage data(imperfectReferenceResult)
"imperfectReferenceResult"


#' Sample data used in the scRNAseq to reference conversion procedurem
#'
#' @format It contains scRNAseq data as well as the spatial data used to mimic
#' the procedure of generating reference based on scRNAseq data.
#'
#' @source \url{https://assets.nemoarchive.org/dat-61kfys3}
#'
#' @usage data(scRNAseq_to_reference)
"scRNAseq_to_reference"



#' Example for allowing user to input their own reference matrix
#'
#' @format It contains reference in the format of cell type x gene. The purpose
#' is to use this dataset to demonstrate that user can input their own reference.
#'
#' @source \url{https://assets.nemoarchive.org/dat-61kfys3}
#'
#' @usage data(reference_ready)
"reference_ready"
