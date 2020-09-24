test_that("multiplication works", {

  p_1 = rnorm(100)
  p_1[1:10]=p_1[1:10]+2
  p_1 = 1-pnorm(p_1)
  p_2 = rnorm(100)
  p_2[1:10]=p_2[1:10]+2
  p_2 = 1-pnorm(p_2)
  out = SIM(p_1, p_2, alpha = 0.05, SIM_F_0_method=1, method_lambda=1, method_t=1)

  if(is.character(out)){
    expect_match(out, 'no rejection')
  }else{
    a_in_b = all(out %in% 1:100)
    expect_equal(a_in_b, TRUE)
  }

})
