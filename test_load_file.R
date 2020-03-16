#包testthat
#library("testthat")

#通过源文件检索
test_that(
  "parents", {
    expect_equal(parents[["GO:0000001"]], c("GO:0048308", "GO:0048311"))
    expect_equal(parents[["GO:0000039"]], NULL)
    expect_type(parents[["GO:0000003"]], "character")
  }
)

#与原作者数据 graph 比对
test_that(
  "children", {
    expect_equal(children[["GO:0000082"]], c("GO:0000083",
                                             "GO:0031575", "GO:2000045"))
    expect_equal(children[["GO:0001304"]], c("GO:0001307", "GO:0001308"))
    expect_type(children[["GO:0000001"]], "character")
  }
)

#与原作者数据 go_annotations['ancestors']
test_that(
  "ancestors", {
    tmp_1 <- ancestors[["GO:0000183"]]
    tmp_1 <- tmp_1[order(tmp_1)]
    expect_equal(tmp_1, 
                 c("GO:0000183", "GO:0006139", "GO:0006342", "GO:0006350",
                   "GO:0006351", "GO:0006355", "GO:0006807", "GO:0008150",
                   "GO:0008152", "GO:0009058", "GO:0009059", "GO:0009889",
                   "GO:0009890", "GO:0009892", "GO:0009987", "GO:0010467",
                   "GO:0010468", "GO:0010556", "GO:0010558", "GO:0010605",
                   "GO:0010629", "GO:0016070", "GO:0016458", "GO:0016481",
                   "GO:0019219", "GO:0019222", "GO:0031323", "GO:0031324",
                   "GO:0031326", "GO:0031327", "GO:0032774", "GO:0034641",
                   "GO:0034645", "GO:0040029", "GO:0043170", "GO:0044237",
                   "GO:0044238", "GO:0044249", "GO:0044260", "GO:0045449",
                   "GO:0045814", "GO:0045892", "GO:0045934", "GO:0048519",
                   "GO:0048523", "GO:0050789", "GO:0050794", "GO:0051171",
                   "GO:0051172", "GO:0051252", "GO:0051253", "GO:0060255",
                   "GO:0065007", "GO:0080090", "GO:0090304", "GO:2000112",
                   "GO:2000113")
                 )
    tmp_2 <- ancestors[["GO:0009142"]]
    tmp_2 <- tmp_2[order(tmp_2)]
    expect_equal(tmp_2,
                 c("GO:0006139", "GO:0006753", "GO:0006807", "GO:0008150",
                   "GO:0008152", "GO:0009058", "GO:0009117", "GO:0009141",
                   "GO:0009142", "GO:0009165", "GO:0009987", "GO:0034404",
                   "GO:0034641", "GO:0034654", "GO:0044237", "GO:0044238",
                   "GO:0044249", "GO:0044271", "GO:0044281", "GO:0044283",
                   "GO:0055086")
                 )
    expect_type(ancestors[["GO:0000001"]], "character")
  }
)

##与原作者数据 go_annotations['genes']
test_that(
  "final_annotations", {
    tmp_1 <- final_annotations[["GO:0000183"]]
    tmp_1 <- tmp_1[order(tmp_1)]
    expect_equal(tmp_1,
                 c("S000000017", "S000000890", "S000001161", "S000001718",
                   "S000001727", "S000002200", "S000002517", "S000002847",
                   "S000002856", "S000003005", "S000003612", "S000003663",
                   "S000004275", "S000004441", "S000004577", "S000004791",
                   "S000004876", "S000005274", "S000005364", "S000005366",
                   "S000005735", "S000005770", "S000005831", "S000005936")
                 )
    tmp_2 <- final_annotations[["GO:0031383"]]
    tmp_2 <- tmp_2[order(tmp_2)]
    expect_equal(tmp_2,
                 c("S000000532", "S000000951", "S000003944", "S000004219",
                   "S000004845", "S000005105")
                 )
    expect_type(final_annotations[["GO:0000001"]], "character")
  }
)