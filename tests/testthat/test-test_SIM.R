test_that("main function works", {

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

test_that("lambda selection works", {

  lambda_set = NULL
  for(i in 1:100){
    p_1 = rnorm(100)
    p_1[1:10]=p_1[1:10]+2
    p_1 = 1-pnorm(p_1)

    option = list()
    option$method_selecting_lambda = 2
    option$search_grid_for_lambda = seq(0.1, 0.8, 0.05)
    out = SIM_lambda_selection(p_1, SIM_F_0_method=1, mu_theta=0, sigma_theta=1, option)
    lambda_set = c(lambda_set, out)
  }
  U = unique(lambda_set)
  a_b = (length(U)>1)
  expect_equal(a_b, TRUE)

})
