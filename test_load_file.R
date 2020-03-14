#library("testthat")

test_that(
  "parents", {
    expect_equal(parents[["GO:0000001"]], c("GO:0048308", "GO:0048311"))
    expect_equal(parents[["GO:0000039"]], NULL)
    expect_type(parents[["GO:0000003"]], "character")
  }
)
