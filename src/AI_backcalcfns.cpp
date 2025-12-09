// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat normalize_mat(arma::mat m)
{
  const double eps = 1e-10;
  arma::mat out = m;
  arma::vec row_sums = arma::sum(m, 1);
  for (unsigned int i = 0; i < m.n_rows; ++i)
  {
    if (row_sums[i] > eps)
    {
      out.row(i) = (m.row(i) / row_sums[i]);
    }
    else
    {
      out.row(i).zeros();
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat normalize_mat_1(arma::mat m)
{
  const double eps = 1e-10;
  arma::mat out = m;
  double s = arma::accu(m);
  if (s > eps)
  {
    for (unsigned int i = 0; i < m.n_rows; ++i)
    {
      out.row(i) = (m.row(i) / s);
    }
  }
  else
  {
    out.zeros();
  }
  return out;
}

// [[Rcpp::export]]
arma::cube trans_mat(const arma::vec &d1, const arma::vec &d2, const arma::vec &d3, const arma::vec &d4, const arma::vec &q)
{
  int t_max = d1.n_elem; // epidemic times
  int ncd4_st = q.n_elem;
  int nst = 2 * ncd4_st + 1;

  arma::cube trans_mats(nst, nst, t_max);

  arma::mat empty_mat(nst, nst, arma::fill::eye);

  for (int t = 0; t < t_max; ++t)
  {
    trans_mats.slice(t) = empty_mat; // Initializing...

    // Creating vector of dx probabilities at time t
    arma::vec d(ncd4_st);
    d(0) = d1(t);
    d(1) = d2(t);
    d(2) = d3(t);
    d(3) = d4(t);

    for (int j = 0; j < ncd4_st; ++j)
    {
      trans_mats(j, j, t) = (1 - q(j)) * (1 - d(j));
      trans_mats(j, j + 1, t) = q(j) * (1 - d(j));
      trans_mats(j, j + ncd4_st + 1, t) = d(j);
    }
  }
  return trans_mats;
}

// [[Rcpp::export]]
arma::mat prev_t(int t, const arma::vec &CD4_prev, const arma::vec &h, const arma::cube &trans_mats)
{
  int t_max = h.n_elem;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2 * ncd4_st + 1;

  // pre-allocate with zeros
  arma::mat prev(t_max, nst, arma::fill::zeros);

  // initialise vector of infs
  arma::vec infs(nst, arma::fill::zeros);
  infs(0) = h(t);

  if (t == 0)
  {
    // initial prevalence case
    arma::vec prev0 = join_cols(CD4_prev, arma::zeros(ncd4_st + 1));
    prev.row(t) = (prev0.t() * trans_mats.slice(0)) + infs.t();
  }
  else
  {
    prev.row(t) = infs.t(); // Can change to row/col vector to avoid .t()
  }

  for (int i = t + 1; i < t_max; ++i)
  {
    prev.row(i) = prev.row(i - 1) * trans_mats.slice(i);
  }
  return prev;
}

// [[Rcpp::export]]
arma::mat prev_t_m(int t, const arma::vec &CD4_prev, const arma::mat &om, const arma::cube &trans_mats)
{
  int t_max = om.n_cols;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2 * ncd4_st + 1;

  // pre-allocate with zeros
  arma::mat prev(t_max, nst, arma::fill::zeros);

  // initialise vector of migrations
  arma::vec infs = join_cols(om.col(t), arma::zeros(ncd4_st + 1));

  if (t == 0)
  {
    // initial prevalence case
    arma::vec prev0 = join_cols(CD4_prev, arma::zeros(ncd4_st + 1));
    prev.row(t) = (prev0.t() * trans_mats.slice(0)) + infs.t();
  }
  else
  {
    prev.row(t) = infs.t(); // Can change to row/col vector to avoid .t()
  }

  for (int i = t + 1; i < t_max; ++i)
  {
    prev.row(i) = prev.row(i - 1) * trans_mats.slice(i);
  }
  return prev;
}

// [[Rcpp::export]]
List forward_calculate(const arma::vec &CD4_prev, const arma::vec &h, const arma::vec &d1, const arma::vec &d2, const arma::vec &d3, const arma::vec &d4, const arma::vec &q)
{
  int t_max = h.n_elem;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2 * ncd4_st + 1;

  arma::cube out(t_max, nst, t_max, arma::fill::zeros);
  arma::cube tmats = trans_mat(d1, d2, d3, d4, q);
  arma::rowvec empty_vec(nst - ncd4_st, arma::fill::zeros);

  arma::mat HIV_diags(t_max, t_max);
  arma::mat AIDS_diags(t_max, t_max);
  arma::cube exp_p(nst - ncd4_st - 1, t_max, t_max, arma::fill::zeros);

  for (int t = 0; t < t_max; ++t)
  {
    out.slice(t) = prev_t(t, CD4_prev, h, tmats);
    // Obtaining diagnoses
    arma::mat X = out.slice(t);
    arma::mat X1 = X(arma::span::all, arma::span(ncd4_st, nst - 1));
    arma::mat diags = arma::join_cols(empty_vec, arma::diff(X1));
    arma::mat HIV_d = diags.cols(1, ncd4_st);
    HIV_diags.row(t) = arma::sum(HIV_d, 1).t();
    AIDS_diags.row(t) = arma::sum(diags.col(0), 1).t();
    exp_p.slice(t) = normalize_mat(HIV_d).t();
  }

  // Now I need to return all quantities of interest
  return List::create(
      _["prevalence"] = out,
      _["HIV.diagnoses"] = HIV_diags,
      _["AIDS.diagnoses"] = AIDS_diags,
      _["expected.p"] = exp_p);
}

// [[Rcpp::export]]
List forward_calculate_m(const arma::vec &CD4_prev, const arma::mat &om, const arma::vec &d1, const arma::vec &d2, const arma::vec &d3, const arma::vec &d4, const arma::vec &q)
{
  int t_max = om.n_cols;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2 * ncd4_st + 1;

  arma::cube out(t_max, nst, t_max, arma::fill::zeros);
  arma::cube tmats = trans_mat(d1, d2, d3, d4, q);
  arma::rowvec empty_vec(nst - ncd4_st, arma::fill::zeros);

  arma::mat HIV_diags(t_max, t_max);
  arma::mat AIDS_diags(t_max, t_max);
  arma::cube exp_p(nst - ncd4_st - 1, t_max, t_max, arma::fill::zeros);
  for (int t = 0; t < t_max; ++t)
  {
    out.slice(t) = prev_t_m(t, CD4_prev, om, tmats);
    // Obtaining diagnoses
    arma::mat X = out.slice(t);
    arma::mat X1 = X(arma::span::all, arma::span(ncd4_st, nst - 1));
    arma::mat diags = arma::join_cols(empty_vec, arma::diff(X1));
    arma::mat HIV_d = diags.cols(1, ncd4_st);
    HIV_diags.row(t) = arma::sum(HIV_d, 1).t();
    AIDS_diags.row(t) = arma::sum(diags.col(0), 1).t();
    exp_p.slice(t) = normalize_mat(HIV_d).t();
  }

  // Now I need to return all quantities of interest
  return List::create(
      _["prevalence"] = out,
      _["HIV.diagnoses"] = HIV_diags,
      _["AIDS.diagnoses"] = AIDS_diags,
      _["expected.p"] = exp_p);
}

// [[Rcpp::export]]
arma::cube temp_diags_fun(List diagnoses)
{
  arma::mat diag = diagnoses(1);
  arma::cube prev = diagnoses(3);
  int nst = prev.n_rows;
  int nquar = diag.n_rows;
  arma::cube out(nquar, nst, nquar, arma::fill::zeros);
  for (int i = 0; i < nquar; ++i)
  {
    for (int j = 0; j < nquar; ++j)
    {
      for (int k = 0; k < nst; ++k)
      {
        out(j, k, i) = diag(i, j) * prev(k, j, i);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat diag_to_infec_agg(int diag_time, const arma::cube &temp_HIV)
{
  int t_max = temp_HIV.n_rows;
  arma::mat diag_mat(temp_HIV.n_cols, diag_time - 1);

  for (int i = 0; i < (diag_time - 1); ++i)
  {
    arma::mat tmp = temp_HIV.slice(i);
    diag_mat.col(i) = tmp.row(diag_time - 1).t(); // -1 R/C++ different indices
  }

  // Adding 0 from diag_time onwards
  // So that all objects have the same size
  arma::mat zero_mat(temp_HIV.n_cols, t_max - diag_time + 1, arma::fill::zeros);

  // cbind matrices
  arma::mat out = arma::join_rows(diag_mat, zero_mat);

  return out;
}

// [[Rcpp::export]]
arma::cube diag_to_infec_wrap(const arma::cube &temp_HIV)
{
  int t_max = temp_HIV.n_rows;
  arma::cube out(temp_HIV.n_cols, t_max, t_max, arma::fill::zeros);

  for (int i = 1; i < t_max; ++i)
  {
    out.slice(i) = normalize_mat_1(diag_to_infec_agg(i + 1, temp_HIV));
  }
  return out;
}

// [[Rcpp::export]]
arma::vec diag_to_infec_agg1(int diag_time, const arma::cube &prev, int state)
{
  int t_max = prev.n_slices;
  arma::vec out(t_max, arma::fill::zeros);

  for (int i = 0; i < (diag_time - 1); ++i)
  {
    out[i] = prev(diag_time - 2, state - 1, i); //-1 for different R/C++ conventions  //t-2 emulates what Paul did with adding zeros (prevalence 1 quarter before dx)
  }

  return out;
}

// [[Rcpp::export]]
arma::mat diag_to_infec_wrap1(const arma::cube &prev, int state)
{
  int t_max = prev.n_slices;
  arma::mat out(t_max, t_max, arma::fill::zeros);

  for (int i = (state - 1); i < (t_max - 1); ++i)
  {                                                                // -1 for different R/C++ conventions
    arma::rowvec tmp = diag_to_infec_agg1(i + 2, prev, state).t(); //-1 for different R/C++ conventions  //t-2 emulates what Paul did with adding zeros (prevalence 1 quarter before dx)
    out.row(i - state + 1) = tmp / arma::accu(tmp);                // Normalization step
  }

  return out;
}

// The idea is the following
// m1[yr.inf,yr.dx] gives imputed yr.inf of AIDS dx
// so consider all yr.inf < t (1:y) and remove all dx that have happened
// before t (-1:-y)... these are the actually undiagnosed prevalent

// [[Rcpp::export]]
arma::vec prev_fun(arma::mat m1)
{
  int n_c = m1.n_cols - 1;
  int n_r = m1.n_rows;
  arma::vec out(n_r, arma::fill::zeros);
  for (int i = 0; i < n_r - 1; ++i)
  {
    out(i) = arma::accu(m1.submat(0, i + 1, i, n_c));
  }
  return out;
}

// The idea is the same as for prev_fun
// Just that these have cube
// because HIV diagnoses are further stratified by state

// [[Rcpp::export]]
arma::vec prev_fun1(arma::cube c1)
{
  int n_c = c1.n_cols - 1;
  int n_r = c1.n_rows - 1;
  int n_s = c1.n_slices - 1;
  arma::vec out(n_c + 1, arma::fill::zeros);
  for (int i = 0; i < n_c; ++i)
  {
    out(i) = arma::accu(c1.subcube(0, 0, i + 1, n_r, i, n_s));
  }
  return out;
}

// [[Rcpp::export]]
arma::cube trans_mat_rita(const arma::vec &d1, const arma::vec &d2, const arma::vec &d3, const arma::vec &d4, const arma::vec &d5, const arma::vec &q)
{
  int t_max = d1.n_elem; // epidemic times
  int ncd4_st = 7;
  int nst = 13;

  arma::cube trans_mats(nst, nst, t_max);

  arma::mat empty_mat(nst, nst, arma::fill::eye);

  for (int t = 0; t < t_max; ++t)
  {
    trans_mats.slice(t) = empty_mat; // Initializing...

    // Creating vector of dx probabilities at time t
    arma::vec d(5);
    d(0) = d1(t);
    d(1) = d2(t);
    d(2) = d3(t);
    d(3) = d4(t);
    d(4) = d5(t);

    trans_mats(0, 0, t) = 0;
    trans_mats(0, 1, t) = (1 - q(0)) * (1 - d(0));
    trans_mats(0, 2, t) = q(0) * (1 - d(0));
    trans_mats(0, 8, t) = d(0);

    trans_mats(1, 1, t) = 0;
    trans_mats(1, 3, t) = (1 - q(0)) * (1 - d(0));
    trans_mats(1, 4, t) = q(0) * (1 - d(0));
    trans_mats(1, 8, t) = d(0);

    trans_mats(2, 2, t) = 0;
    trans_mats(2, 4, t) = (1 - q(1)) * (1 - d(0));
    trans_mats(2, 5, t) = q(1) * (1 - d(0));
    trans_mats(2, 8, t) = d(0);

    for (int j = 3; j < ncd4_st; ++j)
    {
      trans_mats(j, j, t) = (1 - q(j - 3)) * (1 - d(j - 2));
      trans_mats(j, j + 1, t) = q(j - 3) * (1 - d(j - 2));
      trans_mats(j, j + 6, t) = d(j - 2);
    }
  }
  return trans_mats;
}

// [[Rcpp::export]]
arma::mat prev_t_rita(int t, const arma::vec &CD4_prev, const arma::vec &h, const arma::cube &trans_mats)
{
  int t_max = h.n_elem;
  int ncd4_st = 7;
  int nst = 13;

  arma::vec empty_vec(ncd4_st - 1, arma::fill::zeros);
  arma::vec prev0 = arma::join_cols(CD4_prev, empty_vec);
  arma::mat prev(t_max, nst, arma::fill::zeros);

  arma::vec infs(nst, arma::fill::zeros);
  infs(0) = h(t);

  prev.row(t) = infs.t(); // Can change to row/col vector to avoid .t()
  if (t == 0)
  { // In case of initial prevalence
    prev.row(t) = (infs + prev0).t();
  }

  for (int i = t + 1; i < t_max; ++i)
  {
    prev.row(i) = prev.row(i - 1) * trans_mats.slice(i);
  }

  // Rcpp::Rcout << prev << std::endl;

  return prev;
}

// [[Rcpp::export]]
List forward_calculate_rita(const arma::vec &CD4_prev, const arma::vec &h, const arma::vec &d1, const arma::vec &d2, const arma::vec &d3, const arma::vec &d4, const arma::vec &d5, const arma::vec &q)
{
  int t_max = h.n_elem;
  int ncd4_st = 7;
  int nst = 13;

  arma::cube out(t_max, nst, t_max, arma::fill::zeros);
  arma::cube tmats = trans_mat_rita(d1, d2, d3, d4, d5, q);
  arma::rowvec empty_vec(nst - ncd4_st, arma::fill::zeros);

  arma::mat HIV_diags(t_max, t_max);
  arma::mat AIDS_diags(t_max, t_max);
  arma::cube exp_p(nst - ncd4_st - 1, t_max, t_max, arma::fill::zeros);

  for (int t = 0; t < t_max; ++t)
  {
    out.slice(t) = prev_t_rita(t, CD4_prev, h, tmats);
    // Obtaining diagnoses
    arma::mat X = out.slice(t);
    arma::mat X1 = X(arma::span::all, arma::span(ncd4_st, nst - 1));
    // Rcpp::Rcout << X << std::endl;
    arma::mat diags = arma::join_cols(empty_vec, arma::diff(X1));
    arma::mat HIV_d = diags.cols(1, 5);
    HIV_diags.row(t) = arma::sum(HIV_d, 1).t();
    AIDS_diags.row(t) = arma::sum(diags.col(0), 1).t();
    exp_p.slice(t) = normalize_mat(HIV_d).t();
  }

  // Now I need to return all quantities of interest
  return List::create(
      _["prevalence"] = out,
      _["HIV.diagnoses"] = HIV_diags,
      _["AIDS.diagnoses"] = AIDS_diags,
      _["expected.p"] = exp_p);
}