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

