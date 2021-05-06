


context("testPSID")


test_that("check testPSID function", {
  
    n=20
    attr = 2
    td = testPSID(N=n,N.attr=attr)
    
    expect_true(all(names(td) == c("famvars1985","famvars1986", "IND2019ER")))
    expect_that(td, is_a("list"))
    expect_true(nrow(td$famvars1985) == n)
    expect_true(nrow(td$famvars1986) == 2*n - attr)
    expect_true(all(unlist(lapply(td,function(x) !is.null(x)))))

    expect_true(nrow( subset(td$IND2019ER,ER30498!=0) ) == 2*n-attr)
  
} )

