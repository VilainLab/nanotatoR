### Retrieves the genes for a given term/terms.
# Input: terms
# Output:genes and terms

library(nanotatoR)
test_extract_geneList <- function() {
  ### Checking a positive case
  result <- gene_list_generation(
    method = "Single", term = terms,
    returnMethod_GeneList = "dataFrame"
  )
  expect_equal(nrow(result), 305)
  ### Checking a negative case-wrong method
  result1 <- result <- gene_list_generation(
    method = "Triple", term = terms,
    returnMethod_GeneList = "dataFrame"
  )
  expect_equal(nrow(result1), 0)
  ### Checking a negative case-returnMethod
  result <- gene_list_generation(
    method = "Single", term = terms,
    returnMethod_GeneList = "Multiple"
  )
}
