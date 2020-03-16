#包testthat

test_that(
  "parents", {
    expect_equal(parents[["GO:0000001"]], c("GO:0048308", "GO:0048311"))
    expect_equal(parents[["GO:0000039"]], NULL)
    expect_type(parents[["GO:0000003"]], "character")
  }
)

test_that(
  "children", {
    expect_equal(children[["GO:0000082"]], c("GO:0000083",
                                             "GO:0031575", "GO:2000045"))
    expect_equal(children[["GO:0001304"]], c("GO:0001307", "GO:0001308"))
    expect_type(children[["GO:0000001"]], "character")
  }
)
