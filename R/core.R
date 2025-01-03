#' Parametric method for estimating the null distribution
#'
#' @param t vector of thresholds for SIM p-values p(theta)
#' @param p_theta vector of 1D SIM p-values; here is p(theta)
#' @param SIM_F_0_method 1: parametric method; 2: nonparametric method
#' @param mu_0_theta the mean of the single-index p-value
#' @param sigma_0_theta the sd of the single-index p-valu
#'
#' @return the null distribution, same size as t
#' @export
SIM_estimated_F_0 <- function(t, p_theta, SIM_F_0_method, mu_0_theta, sigma_0_theta){

  if(SIM_F_0_method==1){ #parametric method
    hat_F_0_t = stats::pnorm( (stats::qnorm(t) - mu_0_theta)/sigma_0_theta )
  }else if(SIM_F_0_method == 2){ #nonparametric method
    hat_F_0_t = nonparametric_estimated_F_0(t, p_theta)
  }else{
    stop('Incorrect method')
  }
  return(hat_F_0_t)
}



#' Nonparametric method for estimating the null distribution
#'
#' @param t vector of thresholds for p-values
#' @param p_vec vector of 1D p-values
#'
#' @return the null distribution, same size as t
#' @export
nonparametric_estimated_F_0 <- function(t, p_vec){

  denominator = 2*sum( p_vec > 0.5 & p_vec <= 1 ) + sum( p_vec == 0.5 )
  #
  func <- function(t_i){
    if(0 <= t_i & t_i <= 0.5){
      out = sum( p_vec >= (1-t_i) & p_vec <= 1 )/denominator
    }else{
      out = 1 - sum( p_vec >= t_i & p_vec <= 1 )/denominator
    }
    return(out)
  }
  #
  hat_F_0_t = sapply(t, func)
  hat_F_0_t = pmax(0, pmin(hat_F_0_t, 1))  # in [0, 1]
  return(hat_F_0_t)
}



#' Function for estimating the null proportion
#'
#' @param lambda scalar, auxiliary parameter
#' @param p_vec a vector of 1D p-values
#' @param hat_F_0_lambda \eqn{\hat{F}_0(\lambda)}
#'
#' @return \eqn{\hat\pi_0(\lambda)}
#' @export
estimated_pi_0 <- function(lambda, p_vec, hat_F_0_lambda){

  m= length(p_vec)
  R = sum( p_vec <= lambda )    # number of rejections using threshold lambda
  hat_pi_0 = (m - R) / ( (1 - hat_F_0_lambda)*m )
  hat_pi_0 = min(hat_pi_0, 1)  # bounded by 1

  return(hat_pi_0)
}


#' Function for estimating the FDR
#'
#' @param t vector of thresholds for p-values
#' @param p_vec vector of 1D p-values
#' @param hat_pi_0 scalar, \eqn{\hat{\pi}_0}
#' @param hat_F_0_t the null distribution, same size as t
#'
#' @return the estimated FDR, same size as t
#' @export
estimated_FDR <- function(t, p_vec, hat_pi_0, hat_F_0_t){

  m = length(p_vec)
  R_t <- sapply(t, function(x){ sum( p_vec <= x) })
  hat_V_t = m * hat_pi_0 * hat_F_0_t
  hat_V_t = pmin(hat_V_t, R_t)  #V(t) <= R(t)
  hat_FDR_t =  hat_V_t/pmax(R_t, 1)  # in [0, 1]

  return(hat_FDR_t)
}


#' Function to compute single index p-value
#'
#' @param theta projection direction in \eqn{[0, \pi/2]}
#' @param p_1 component p-value
#' @param p_2 component p-value
#'
#' @return SIM_p: single index modulated p-value
#' @export
SIM_p_value <- function(p_1, p_2, theta){

  if(theta == 0){
    SIM_p = p_1
  }else if(theta == pi/2){
    SIM_p = p_2
  }else{
    SIM_p = stats::pnorm( cos(theta)*stats::qnorm(p_1) + sin(theta)*stats::qnorm(p_2) )
  }

  return(SIM_p)
}


#' Function to estimate the mean of SIM p-value under true null
#'
#' @param SIM_p single-index SIM_p-value
#' @return \eqn{\hat\mu_0(\theta)}
#' @export
SIM_estimated_mu_0_theta <- function(SIM_p){

  z = stats::qnorm(SIM_p)
  hat_mu = mean(z[z>-2 & z<2])
  return(hat_mu)
}


#' Function to estimate the sd of the SIM p-value under true null
#'
#' @param SIM_p single-index SIM p-values
#' @return the sd of the SIM p-value
#' @export
SIM_estimated_sigma_0_theta <- function(SIM_p){

  SIM_p = pmin(SIM_p, 0.999999)
  z = stats::qnorm(SIM_p)
  positive = z[z>0]
  z_trancate = c(positive, positive*(-1))
  hat_sigma = stats::sd(z_trancate)

  return(hat_sigma)
}

#' Function: select lambda used in estimating the null proportion
#'
#' @param p_theta vector of SIM p-values
#' @param SIM_F_0_method 1: parametric method; 2: nonparametric method
#' @param mu_theta \eqn{\mu_0(\theta)}
#' @param sigma_theta \eqn{\sigma_0(\theta)}
#' @param option options for choosing lambda
#' @export
SIM_lambda_selection <-  function(p_theta, SIM_F_0_method, mu_theta, sigma_theta, option){

  if(option$method_selecting_lambda == 1){ #specified value
    lambda_select = option$specified_lambda
  }

  if(option$method_selecting_lambda == 2){ #adaptive method

    search_grid_for_lambda = option$search_grid_for_lambda
    search_grid_for_lambda = sort(search_grid_for_lambda)
    n_lambda_grid = length(search_grid_for_lambda)

    func <- function(lambda_i, pi_0_i_minus_1){
      F_0 = SIM_estimated_F_0(lambda_i, p_theta, SIM_F_0_method, mu_theta, sigma_theta)
      pi_0_i = estimated_pi_0(lambda_i, p_theta, F_0)
      return( list((pi_0_i >= pi_0_i_minus_1), pi_0_i) )
    }
    # initial value
    i = 0
    lambda_0 = 0
    F_0 = SIM_estimated_F_0(lambda_0, p_theta, SIM_F_0_method, mu_theta, sigma_theta)
    pi_0_i_minus_1 = estimated_pi_0(lambda_0, p_theta, F_0)

    while( i <= (n_lambda_grid - 1) ){
      i=i+1
      lambda_i = search_grid_for_lambda[i]
      out = func(lambda_i, pi_0_i_minus_1)
      if( out[[1]] ){
        lambda_select = search_grid_for_lambda[i]
        break
      }else{
        pi_0_i_minus_1 = out[[2]]
      }
    }

    if(i == n_lambda_grid) lambda_select = search_grid_for_lambda[n_lambda_grid]
  }
  return(lambda_select)
}


#' Function: threshold value for SIM p-valueS
#'
#' @param alpha  error rate
#' @param p_theta vector of SIM p-values
#' @param hat_pi_0_theta estimated \eqn{\pi_0}
#' @param SIM_F_0_method 1: parametric method; 2: nonparametric method.
#' @param mu_theta \eqn{\mu_0(\theta)}
#' @param sigma_theta \eqn{\sigma_0(\theta)}
#' @param option options for choosing the grid
#'
#' @return t_alpha_SIM: threshold value for SIM p-values
#' @export
SIM_threshold_for_p <- function(alpha, p_theta, hat_pi_0_theta, SIM_F_0_method, mu_theta, sigma_theta, option){

  if(option$method_t_alpha == 1){
    # step I: find the minimum value
    x_left = 0
    x_right = 1
    de = x_right-x_left
    while(de > 0.00001){
      length = x_right-x_left
      x_left1 = x_left+length*(1-0.618)
      x_right1 = x_right-length*(1-0.618)

      hat_F_0_x_left1  = SIM_estimated_F_0(x_left1, p_theta, SIM_F_0_method, mu_theta, sigma_theta)
      hat_F_0_x_right1 = SIM_estimated_F_0(x_right1, p_theta, SIM_F_0_method, mu_theta, sigma_theta)

      if(estimated_FDR(x_left1,  p_theta, 1, hat_F_0_x_left1) >
         estimated_FDR(x_right1, p_theta, 1, hat_F_0_x_right1)){
        x_left = x_left1
      }else{
        x_right = x_right1
      }

      de = x_right-x_left
    }

    # step II: find the cutoff for FDR control
    hat_F_0_x_left = SIM_estimated_F_0(x_left, p_theta, SIM_F_0_method, mu_theta, sigma_theta)
    if( estimated_FDR(x_left, p_theta, hat_pi_0_theta, hat_F_0_x_left) > alpha ){
      min = 0
    }else{
      min = x_left
      max_d = 1
      diff = 1
      mid = (min+max_d)/2
      while(diff > 0.000001){
        hat_F_0_mid = SIM_estimated_F_0(mid, p_theta, SIM_F_0_method, mu_theta, sigma_theta)
        est_FDR = estimated_FDR(mid, p_theta, hat_pi_0_theta, hat_F_0_mid)

        if(est_FDR > alpha){
          max_d = mid
        }else{
          min = mid
        }

        mid = (min+max_d)/2
        diff = abs(max_d-min)
      }
    }
    t_alpha_SIM = min

  }else if(option$method_t_alpha == 2){
    #------------use the sorted p_value as the t_grid--------------------
    t_grid=sort(p_theta)

    hat_F_0 = SIM_estimated_F_0(t_grid, p_theta, SIM_F_0_method, mu_theta, sigma_theta)

    hat_FDR = estimated_FDR(t_grid, p_theta, hat_pi_0_theta, hat_F_0)

    index_set = which( hat_FDR <= alpha )

    if(length(index_set)>0){ # not empty
      t_alpha_SIM = t_grid[ max(index_set) ]
    }else{
      t_alpha_SIM = 0
    }
  }

  return(t_alpha_SIM)
}



#' Function: optimal selection of theta
#'
#' @param p_1 the auxiliary p-values
#' @param p_2 the primary p-values
#' @param alpha the significance level
#' @param SIM_F_0_method the method for estimating the null distribution
#' @param option the grid for theta
#' @return the estimated projection direction
#' @export
SIM_theta_0_selection <- function(p_1, p_2, alpha, SIM_F_0_method, option){

  search_grid_for_theta_0 = option$search_grid_for_theta_0
  par = length(search_grid_for_theta_0)

  func <- function(theta){
    SIM_p = SIM_p_value(p_1, p_2, theta)

    if(SIM_F_0_method == 1){
      mu_theta = 0
      hat_sigma = SIM_estimated_sigma_0_theta(SIM_p)
    }else if(SIM_F_0_method == 2){
      mu_theta = NULL
      hat_sigma = NULL
    }

    out = sum(SIM_p <= SIM_threshold_for_p(alpha, SIM_p, 1, SIM_F_0_method, mu_theta, hat_sigma, option))
    return(out)
  }

  N_rej <- sapply(search_grid_for_theta_0, func)
  i_hat = which.max(N_rej)
  # refine the method
  if(i_hat == 1){
    x_left = search_grid_for_theta_0[1]
    x_right = search_grid_for_theta_0[2]
  }else if(i_hat == par){
    x_left = search_grid_for_theta_0[par-1]
    x_right = search_grid_for_theta_0[par]
  }else{
    x_left = search_grid_for_theta_0[i_hat-1]
    x_right = search_grid_for_theta_0[i_hat+1]
  }


  theta1 = seq(x_left, x_right, (x_right-x_left)/par)
  N_rej <- sapply(theta1, func)
  k_hat = which.max(N_rej)
  OUT = theta1[k_hat]

  return(OUT)
}


#' The main function in SIM paper
#'
#' @param p_1 the auxiliary  p-values
#' @param p_2 the primary p-values
#' @param alpha the significance level
#' @param SIM_F_0_method method for estimating the null distrbution (1: parametric; 2: nonparametic method)
#' @param method_lambda method for searching the lambda used in null proportion (1: specified; 2: adaptive)
#' @param method_t method for searching the p-value threshold (1: bisection; 2: grid search)
#' @return the indices of the hypotheses rejected
#' @examples
#' p_1 = rnorm(100)
#' p_1[1:10]=p_1[1:10]+2
#' p_1 = 1-pnorm(p_1)
#' p_2 = rnorm(100)
#' p_2[1:10]=p_2[1:10]+2
#' p_2 = 1-pnorm(p_2)
#' out = SIM(p_1, p_2, alpha = 0.1, SIM_F_0_method=1, method_lambda=1, method_t=2)
#' @export
SIM <- function(p_1, p_2, alpha, SIM_F_0_method=1, method_lambda=1, method_t=1){

  # theta-grid
  option = list()
  par1 = 30
  option$search_grid_for_theta_0 = seq(0, pi/2, (pi/2)/par1)


  # lambda-selection
  # input(' input option_method_selecting_lambda (1: specified; 2: adaptive) = ')
  option$method_selecting_lambda = method_lambda

  if(option$method_selecting_lambda == 1){
    option$specified_lambda = 0.5
  }else if(option$method_selecting_lambda == 2){
    search_grid_for_lambda = seq(0.1, 0.8, (0.8-0.1)/15)
    option$search_grid_for_lambda = c(0.02, 0.04, 0.06, 0.08, search_grid_for_lambda)
  }

  # t_grid search
  # input(' input option.method_t_alpha (1: bisection; 2: grid search)
  option$method_t_alpha = method_t

  theta_I = SIM_theta_0_selection(p_1, p_2, alpha, SIM_F_0_method, option)
  SIM_p = SIM_p_value(p_1, p_2, theta_I)

  if(SIM_F_0_method == 1){
    mu_0_theta = 0
    hat_sigma_0_theta = SIM_estimated_sigma_0_theta(SIM_p)
  }else{
    mu_0_theta = NULL
    hat_sigma_0_theta = NULL
  }


  lambda_select_I = SIM_lambda_selection(SIM_p, SIM_F_0_method, mu_0_theta, hat_sigma_0_theta, option)

  F_0 = SIM_estimated_F_0(lambda_select_I, SIM_p, SIM_F_0_method, mu_0_theta, hat_sigma_0_theta)
  hat_pi_0_lambda_I = estimated_pi_0(lambda_select_I, SIM_p, F_0)

  t_alpha_I = SIM_threshold_for_p(alpha, SIM_p, hat_pi_0_lambda_I,
                                  SIM_F_0_method, mu_0_theta, hat_sigma_0_theta, option)

  R_I = sum(SIM_p <= t_alpha_I)
  Index_I = which(SIM_p <= t_alpha_I)

  if(length(Index_I)>0){
    return(Index_I)
  }else{
    return('no rejection')
  }

}
