# context("Manuka11")
# 
# test_that("Manuka dataset", {
# 
#   dset <- Manuka11()
#   
#   expect_equal(class(dset), "character")
#   expect_equal(length(dset), 1)
#   
#   # expect_equal(class(dset), "character")
#   # expect_output(str(dset), "List of 7")
#   # expect_equal(length(dset), 7)
#   # expect_equal(names(dset), c("genon", "depth_Ref", "depth_Alt", "chrom", "pos", "indID",  "SNP_Names"))
#   # 
#   # ## dataset byte size
#   # ## OS X
#   # if(length(grep("darwin", version$os))) {
#   #   expect_equal(as.numeric(object.size(dset)), 2522184)
#   # }
#   # 
#   # if(length(grep("mingw32", version$os))) {
#   #   expect_equal(as.numeric(object.size(dset)), 2501640)
#   # }
#   # 
#   # ## dimensions and length
#   # dset_dims <- lapply(dset, function(x)if(is.null(dim(x))){length(x)} else {dim(x)})
#   # 
#   # expect_equal(dset_dims, structure(list(genon = c(181L, 680L), depth_Ref = c(181L, 680L ),
#   #                                        depth_Alt = c(181L, 680L),
#   #                                        chrom = 680L, pos = 680L, indID = 181L, SNP_Names = 680L),
#   #                                   .Names = c("genon", "depth_Ref", "depth_Alt",  "chrom", "pos", "indID", "SNP_Names")))
# })
NULL
