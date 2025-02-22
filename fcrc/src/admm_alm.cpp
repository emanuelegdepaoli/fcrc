// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;

// Soft thresholding operator
arma::vec soft_thre_int(arma::vec const x, double const par){
   double const x_norm = norm(x, 2);
   arma::vec out = zeros(x.n_rows);
   if((1-par/x_norm) > 0){
      out = x*(1-par/x_norm);
   }
   return(out);
}

// Alternative Direction Method of Multipliers (ADMM) for Group Lasso
// inspired by https://stanford.edu/~boyd/papers/admm/lasso/lasso.html
// and related paper Boyd, Stephen, Neal Parikh, and Eric Chu. 
// Distributed optimization and statistical learning via the alternating direction method of multipliers. 
// Now Publishers Inc, 2011.
// scaled augmented lagrangian
// rho can vary for max. 50 iterations and then at least 50 more iterations are done
// y, X centered
// loss function: 1/2*t(theta)*xtx*theta - t(theta)*xty + lambda*sum(||theta_j||_2)
// s.t theta_j = beta_j

// [[Rcpp::export]]
List admm_grplasso_path_int(arma::mat const xty, arma::mat const xtx, arma::uvec const index,
                            arma::vec const lambda_path, arma::vec beta_start, arma::vec u_start, 
                            bool const varying_rho = true, double rho = 1, 
                            double const abs_tol = 1e-4, double const rel_tol = 1e-2,
                            int const max_iter = 1000L){
  
  // abs_tol and rel_tol should be 1e-3 or 1e-4
  // inizialization results path
  int const n_lambdas = lambda_path.n_elem;
  int const n_coef = beta_start.n_elem;
  arma::mat res_beta(n_coef, n_lambdas);
  vec res_niter = zeros(n_lambdas);
  
  for(uword i = 0; i < n_lambdas; ++i){
  
    double const lambda = lambda_path(i);
    uword const p = xtx.n_cols;
    arma::uvec const grps = unique(index);
    uword const n_groups = grps.n_rows;
    arma::vec beta = beta_start;
    arma::vec theta = beta_start;
    arma::vec u = u_start;
    arma::mat xtx_V;
    arma::vec xtx_D;
    eig_sym(xtx_D, xtx_V, xtx);
    arma::mat xtxr_inv = xtx_V*diagmat(pow(xtx_D+rho, -1.0))*xtx_V.t();
    arma::vec objval = zeros(max_iter);
    arma::vec r_norm = zeros(max_iter);
    arma::vec s_norm = zeros(max_iter);
    arma::vec rho_list = zeros(max_iter);
    arma::vec eps_pri = zeros(max_iter);
    arma::vec eps_dual = zeros(max_iter);
    double rho_old = rho;
    arma::vec beta_old = beta;
    
    // iterations start
    uword iter = 1;
    while(iter <= max_iter){
      rho_list(iter-1) = rho;
      beta_old = beta;
      if(varying_rho && (rho_old != rho)){
        u = u*(rho_old/rho);
        xtxr_inv = xtx_V*diagmat(pow(xtx_D+rho, -1.0))*xtx_V.t();
      }
      
      // theta update
      theta = xtxr_inv*(xty+rho*(beta-u));
      
      // beta update
      for(uword j = 0; j < n_groups; ++j){
        arma::uvec const idx = find(index == grps(j));
        beta.elem(idx) = soft_thre_int(theta.elem(idx)+u.elem(idx), lambda/rho);
      }
      
      // u update
      u = u + theta-beta;
      
      // stopping criteria from Boyd et al. 2011
      r_norm(iter-1)= norm(theta-beta, 2); // primal residual
      s_norm(iter-1) = norm(-rho*(beta-beta_old), 2); // dual residual 
      // objval(iter-1) = obj_fun_grplasso_int(xty, xtx, theta, index, lambda);
      if(norm(theta, 2) > norm(-beta, 2)){
        eps_pri(iter-1) = sqrt(p)*abs_tol+rel_tol*norm(theta, 2);
      } else{
        eps_pri(iter-1) = sqrt(p)*abs_tol+rel_tol*norm(-beta, 2); 
      }
      eps_dual(iter-1) = sqrt(p)*abs_tol+rel_tol*rho*norm(u, 2); 
      
      // convergence check
      if(varying_rho){
        if((r_norm(iter-1) < eps_pri(iter-1)) && (s_norm(iter-1) < eps_dual(iter-1)) && iter >= 100) break;
      } else{
        if((r_norm(iter-1) < eps_pri(iter-1)) && (s_norm(iter-1) < eps_dual(iter-1))) break;
      }
      
      // varying penalty parameter
      rho_old = rho;
      if((iter <= 50) && varying_rho){
        if(r_norm(iter-1) > 10*s_norm(iter-1)){
          rho = 2*rho;
        } else if(s_norm(iter-1) > 10*r_norm(iter-1)){
          rho = rho/2;
        } else rho = rho;
      }
      
      iter += 1;
    }
    
    // if max_iter is reached, iter will be max_iter+1
    if(iter > max_iter) iter = max_iter;
    // store results for i-th lambdas
    res_beta.col(i) = beta;
    res_niter(i) = iter;
    beta_start = beta;
    u_start = u;
  }
  
  return(List::create(Named("res_beta") = res_beta, Named("res_niter") = res_niter));
}


// Augmented Lagrangian Method 
// solves Lρ(beta, u) = -t(beta)*J + 0.5*t(beta)*K*beta + lambda*sum(||theta_j||_2) 
// + t(u)*L*beta + rho/2*(||L*beta||_2)^2
// if varying_gamma = T, the penalty parameter associated to the internal ADMM varies each iteration 
// else it varies only at the first internal ADMM iteration

// [[Rcpp::export]]
List alm_cgl_path_int(arma::vec const J, arma::mat const K, arma::mat const L, arma::uvec const index, 
                      arma::vec u, arma::vec u_int, arma::vec beta, arma::vec const lambda_path, double rho, 
                      double const max_rho, double gamma, bool const varying_gamma = true, 
                      double const eps = 1e-5, int const max_iter = 100L, int const max_iter_int = 1000L, 
                      double const abs_tol_int = 1e-6, double const rel_tol_int = 1e-4) {
  
  // inizialization results path
  int const n_lambdas = lambda_path.n_elem;
  int const n_coef = beta.n_elem;
  arma::mat res_beta(n_coef, n_lambdas);
  arma::vec res_niter = zeros(n_lambdas);
  
  for(uword i = 0; i < n_lambdas; ++i){
    
    double const lambda = lambda_path(i);
    arma::vec err = zeros(max_iter);
    double rho_old = rho;
    
    arma::mat xtxr = K + rho*L.t()*L;
    uword const p = xtxr.n_cols;
    arma::uvec const grps = unique(index);
    uword const n_groups = grps.n_rows;
    arma::vec xty = zeros(size(J));
  
    // iterations start
    uword iter = 1;
    while (iter <= max_iter) {
      arma::vec beta_old = beta;
      if (rho_old != rho) {
        xtxr = K + rho*L.t()*L;
      } else {
        xty = J - L.t()*u; 
      }
      
      // beta update
      arma::vec theta = beta;
      arma::vec beta_int = beta;
      arma::mat xtxr_V;
      arma::vec xtxr_D;
      eig_sym(xtxr_D, xtxr_V, xtxr);
      arma::mat xtxr_inv = xtxr_V*diagmat(pow(xtxr_D+gamma, -1.0))*xtxr_V.t();
      arma::vec r_norm = zeros(max_iter_int);
      arma::vec s_norm = zeros(max_iter_int);
      arma::vec eps_pri = zeros(max_iter_int);
      arma::vec eps_dual = zeros(max_iter_int);
      double gamma_old = gamma;
      arma::vec beta_old_int = beta_int;
      
      // internal iterations start
      uword iter_int = 1;
      while(iter_int <= max_iter_int){
        beta_old_int = beta_int;
        if(gamma_old != gamma){
          u_int = u_int*(gamma_old/gamma);
          xtxr_inv = xtxr_V*diagmat(pow(xtxr_D+gamma, -1.0))*xtxr_V.t();
        }
        
        // theta update
        theta = xtxr_inv*(xty+gamma*(beta_int-u_int));
        
        // beta update
        for(uword j = 0; j < n_groups; ++j){
          arma::uvec const idx = find(index == grps(j));
          beta_int.elem(idx) = soft_thre_int(theta.elem(idx)+u_int.elem(idx), lambda/gamma);
        }
        
        // u update
        u_int = u_int + theta-beta_int;
        
        // stopping criteria from Boyd et al. 2011
        r_norm(iter_int-1)= norm(theta-beta_int, 2); // primal residual
        s_norm(iter_int-1) = norm(-gamma*(beta_int-beta_old_int), 2); // dual residual 
        if(norm(theta, 2) > norm(-beta_int, 2)){
          eps_pri(iter_int-1) = sqrt(p)*abs_tol_int+rel_tol_int*norm(theta, 2);
        } else{
          eps_pri(iter_int-1) = sqrt(p)*abs_tol_int+rel_tol_int*norm(-beta_int, 2); 
        }
        eps_dual(iter_int-1) = sqrt(p)*abs_tol_int+rel_tol_int*gamma*norm(u_int, 2); 
        
        // convergence check
        if(varying_gamma | iter == 1){
          if((r_norm(iter_int-1) < eps_pri(iter_int-1)) && (s_norm(iter_int-1) < eps_dual(iter_int-1)) && iter_int >= 100) break;
        } else{
          if((r_norm(iter_int-1) < eps_pri(iter_int-1)) && (s_norm(iter_int-1) < eps_dual(iter_int-1))) break;
        }
        
        // varying penalty parameter
        gamma_old = gamma;
        if(varying_gamma | iter == 1) {
          if(iter_int <= 50){
            if(r_norm(iter_int-1) > 10*s_norm(iter_int-1)){
              gamma = 2*gamma;
            } else if(s_norm(iter_int-1) > 10*r_norm(iter_int-1)){
              gamma = gamma/2;
            } else gamma = gamma;
          }
        }
        
        iter_int += 1;
      }
      // if max_iter is reached, iter will be max_iter+1
      if(iter_int > max_iter_int) iter_int = max_iter_int;
      beta = beta_int;
      
      // u update
      // varying rho according to Bertsekas, 1981 - pag 123
      rho_old = rho;
      err(iter-1) = max(abs(L*beta));
      if (err(iter-1) < eps) {
        break;
      } else if (err(iter-1) > 0.25*max(abs(L*beta_old)) & rho < max_rho) {
        rho = 10*rho_old;
      } else {
        u = u+rho*L*beta;
      }
      
      iter += 1;
    }
    
    if(iter > max_iter) iter = max_iter; // if max_iter is reached, iter will be max_iter+1
    
    // store results for i-th lambdas
    res_beta.col(i) = beta;
    res_niter(i) = iter;
  }
  
  return(List::create(Named("res_beta") = res_beta, Named("res_niter") = res_niter));
}


