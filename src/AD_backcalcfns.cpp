// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
arma::mat Q1_fct(const arma::mat &d, const arma::mat &q, int t, int a)
{
  arma::mat Q(4, 4);
  Q(0, 0) = (1 - d(0, t)) * (1 - q(0, a));
  Q(0, 1) = (1 - d(0, t)) * q(0, a);
  Q(0, 2) = 0;
  Q(0, 3) = 0;
  Q(1, 0) = 0;
  Q(1, 1) = (1 - d(1, t)) * (1 - q(1, a));
  Q(1, 2) = (1 - d(1, t)) * q(1, a);
  Q(1, 3) = 0;
  Q(2, 0) = 0;
  Q(2, 1) = 0;
  Q(2, 2) = (1 - d(2, t)) * (1 - q(2, a));
  Q(2, 3) = (1 - d(2, t)) * q(2, a);
  Q(3, 0) = 0;
  Q(3, 1) = 0;
  Q(3, 2) = 0;
  Q(3, 3) = (1 - d(3, t)) * (1 - q(3, a));
  return Q;
}

// [[Rcpp::export]]
arma::mat Q1_fct_rita(const arma::mat &d, const arma::mat &q, int t, int a)
{
  arma::mat Q(7, 7, arma::fill::zeros);
  Q(0, 1) = (1 - d(0, t)) * (1 - q(0, a));
  Q(0, 2) = (1 - d(0, t)) * q(0, a);
  Q(1, 3) = (1 - d(0, t)) * (1 - q(0, a));
  Q(1, 4) = (1 - d(0, t)) * q(0, a);
  Q(2, 4) = (1 - d(0, t)) * (1 - q(1, a));
  Q(2, 5) = (1 - d(0, t)) * q(1, a);
  Q(3, 3) = (1 - d(1, t)) * (1 - q(0, a));
  Q(3, 4) = (1 - d(1, t)) * q(0, a);
  Q(4, 4) = (1 - d(2, t)) * (1 - q(1, a));
  Q(4, 5) = (1 - d(2, t)) * q(1, a);
  Q(5, 5) = (1 - d(3, t)) * (1 - q(2, a));
  Q(5, 6) = (1 - d(3, t)) * q(2, a);
  Q(6, 6) = (1 - d(4, t)) * (1 - q(3, a));
  return Q;
}

// [[Rcpp::export]]
int age_fct(int a0, int t, int t0, int nquar, int nage)
{
  // Individuals age after quarters 3,7,11 etc
  // So at quarters 4,8,12 etc
  // C++ notation starts at 0 (unlike STANs)

  int a = a0;
  for (int t1 = t0 + 1; t1 < (t + 1); ++t1)
  {
    int rem = (t1 + 1) % 4;
    int addyr = 0;
    if (rem == 1)
    {
      addyr = 1;
    }
    else
    {
      addyr = 0;
    }
    a += addyr;
  }
  return std::min(a, nage - 1);
}

int age_to_agroup(int a)
{
  int ag;
  if (a <= 9)
  {
    ag = 0;
  }
  else if (a > 9 && a <= 19)
  {
    ag = 1;
  }
  else if (a > 19 && a <= 29)
  {
    ag = 2;
  }
  else
  {
    ag = 3;
  }
  return ag;
}

// Function to make any [t][a0][t0] array
// into a [t][a][t0] array
std::vector<std::vector<std::vector<arma::vec>>> a0toa_fct(const std::vector<std::vector<std::vector<arma::vec>>> &array, int nquar, int nage)
{
  std::vector<std::vector<std::vector<arma::vec>>> out(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  arma::vec emptyvec(array[0][0][0].n_elem, arma::fill::zeros);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      for (int t0 = 0; t0 < nquar; ++t0)
      {
        out[t][a][t0] = emptyvec;
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      for (int t = t0; t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        out[t][a][t0] += array[t][a0][t0];
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::cube gof_iter_fct(const arma::vec &H, const arma::cube &d, const arma::mat &q, const arma::mat &prev, int nquar, int nage)
{
  arma::cube exp_diag(nquar, nage, 5);
  const arma::vec empty_vec(5, arma::fill::zeros);
  const arma::vec empty_vec1(4, arma::fill::zeros);
  int nagegrp = d.n_slices;

  std::vector<std::vector<std::vector<arma::mat>>> PA(nquar, std::vector<std::vector<arma::mat>>(nage, std::vector<arma::mat>(nagegrp)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr1(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      for (int ag = 0; ag < nagegrp; ++ag)
      {
        arma::mat d1 = d.slice(ag);
        PA[t][a][ag] = Q1_fct(d1, q, t, a);
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      // Adding initial prevalence at t==0
      arma::vec diag_p(5, arma::fill::zeros);
      arma::vec lat_p(4, arma::fill::zeros);
      arma::vec tmp(4);
      int ag0 = age_to_agroup(a0);

      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;
      diag_arr[t0][a0][t0] = empty_vec;

      if (t0 == 0)
      {
        tmp(0) = prev(0, a0);
        tmp(1) = prev(1, a0);
        tmp(2) = prev(2, a0);
        tmp(3) = prev(3, a0);

        // Progressed initial prevalence
        lat_p = lat_p + (PA[0][a0][ag0].t() * tmp);

        // Diagnosed initial prevalence
        diag_p(0) = tmp(0) * d(0, 0, ag0);
        diag_p(1) = tmp(1) * d(1, 0, ag0);
        diag_p(2) = tmp(2) * d(2, 0, ag0);
        diag_p(3) = tmp(3) * d(3, 0, ag0);
        diag_p(4) = tmp(3) * (1 - d(3, 0, ag0)) * q(3, a0);
      }
      else
      { // Probably dont need this
        diag_p = empty_vec;
        lat_p = empty_vec1;
      }

      lat_arr[t0][a0][t0][0] = H((a0 * nquar) + t0) + lat_p(0);
      lat_arr[t0][a0][t0][1] = lat_p(1);
      lat_arr[t0][a0][t0][2] = lat_p(2);
      lat_arr[t0][a0][t0][3] = lat_p(3);

      diag_arr[t0][a0][t0] = diag_p; // use row_vec instead?

      for (int t = (t0 + 1); t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);

        lat_arr[t][a0][t0] = PA[t][a0][ag].t() * lat_arr[t - 1][a0][t0];

        arma::vec tmpvec(5);

        tmpvec[0] = lat_arr[t - 1][a0][t0][0] * d(0, t, ag);
        tmpvec[1] = lat_arr[t - 1][a0][t0][1] * d(1, t, ag);
        tmpvec[2] = lat_arr[t - 1][a0][t0][2] * d(2, t, ag);
        tmpvec[3] = lat_arr[t - 1][a0][t0][3] * d(3, t, ag);
        tmpvec[4] = lat_arr[t - 1][a0][t0][3] * (1 - d(3, t, ag)) * q(3, a0);

        diag_arr[t][a0][t0] = tmpvec;
      }
    }
  }

  diag_arr1 = a0toa_fct(diag_arr, nquar, nage);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      exp_diag(t, a, 0) = 0;
      exp_diag(t, a, 1) = 0;
      exp_diag(t, a, 2) = 0;
      exp_diag(t, a, 3) = 0;
      exp_diag(t, a, 4) = 0;

      for (int t0 = 0; t0 <= t; ++t0)
      {
        exp_diag(t, a, 0) += diag_arr1[t][a][t0][0];
        exp_diag(t, a, 1) += diag_arr1[t][a][t0][1];
        exp_diag(t, a, 2) += diag_arr1[t][a][t0][2];
        exp_diag(t, a, 3) += diag_arr1[t][a][t0][3];
        exp_diag(t, a, 4) += diag_arr1[t][a][t0][4];
      }
    }
  }

  return exp_diag;
}

// [[Rcpp::export]]
arma::cube gof_iter_fct_rita(const arma::vec &H, const arma::cube &d, const arma::mat &q, const arma::mat &prev, int nquar, int nage)
{
  arma::cube exp_diag(nquar, nage, 6);
  const arma::vec empty_vec(6, arma::fill::zeros);
  const arma::vec empty_vec1(7, arma::fill::zeros);
  int nagegrp = d.n_slices;

  std::vector<std::vector<std::vector<arma::mat>>> PA(nquar, std::vector<std::vector<arma::mat>>(nage, std::vector<arma::mat>(nagegrp)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr1(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      for (int ag = 0; ag < nagegrp; ++ag)
      {
        arma::mat d1 = d.slice(ag);
        PA[t][a][ag] = Q1_fct_rita(d1, q, t, a);
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      // Adding initial prevalence at t==0
      arma::vec diag_p(6, arma::fill::zeros);
      arma::vec lat_p(7, arma::fill::zeros);
      arma::vec tmp(7);
      int ag0 = age_to_agroup(a0);

      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;
      diag_arr[t0][a0][t0] = empty_vec;

      if (t0 == 0)
      {
        tmp(0) = prev(0, a0);
        tmp(1) = prev(1, a0);
        tmp(2) = prev(2, a0);
        tmp(3) = prev(3, a0);
        tmp(4) = prev(4, a0);
        tmp(5) = prev(5, a0);
        tmp(6) = prev(6, a0);

        // Progressed initial prevalence
        lat_p = lat_p + (PA[0][a0][ag0].t() * tmp);

        // Diagnosed initial prevalence
        diag_p(0) = (tmp(0) * d(0, 0, ag0)) + (tmp(1) * d(0, 0, ag0)) + (tmp(2) * d(0, 0, ag0));
        diag_p(1) = tmp(3) * d(1, 0, ag0);
        diag_p(2) = tmp(4) * d(2, 0, ag0);
        diag_p(3) = tmp(5) * d(3, 0, ag0);
        diag_p(4) = tmp(6) * d(4, 0, ag0);
        diag_p(5) = tmp(6) * (1 - d(4, 0, ag0)) * q(3, a0);
      }
      else
      { // Probably dont need this
        diag_p = empty_vec;
        lat_p = empty_vec1;
      }

      lat_arr[t0][a0][t0][0] = H((a0 * nquar) + t0) + lat_p(0);
      lat_arr[t0][a0][t0][1] = lat_p(1);
      lat_arr[t0][a0][t0][2] = lat_p(2);
      lat_arr[t0][a0][t0][3] = lat_p(3);
      lat_arr[t0][a0][t0][4] = lat_p(4);
      lat_arr[t0][a0][t0][5] = lat_p(5);
      lat_arr[t0][a0][t0][6] = lat_p(6);

      diag_arr[t0][a0][t0] = diag_p; // use row_vec instead?

      for (int t = (t0 + 1); t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);
        lat_arr[t][a0][t0] = PA[t][a0][ag].t() * lat_arr[t - 1][a0][t0];

        arma::vec tmpvec(6);

        tmpvec[0] = (lat_arr[t - 1][a0][t0][0] * d(0, t, ag)) + (lat_arr[t - 1][a0][t0][1] * d(0, t, ag)) + (lat_arr[t - 1][a0][t0][2] * d(0, t, ag));
        tmpvec[1] = lat_arr[t - 1][a0][t0][3] * d(1, t, ag);
        tmpvec[2] = lat_arr[t - 1][a0][t0][4] * d(2, t, ag);
        tmpvec[3] = lat_arr[t - 1][a0][t0][5] * d(3, t, ag);
        tmpvec[4] = lat_arr[t - 1][a0][t0][6] * d(4, t, ag);
        tmpvec[5] = lat_arr[t - 1][a0][t0][6] * (1 - d(4, t, ag)) * q(3, a0);

        diag_arr[t][a0][t0] = tmpvec;
      }
    }
  }

  diag_arr1 = a0toa_fct(diag_arr, nquar, nage);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      exp_diag(t, a, 0) = 0;
      exp_diag(t, a, 1) = 0;
      exp_diag(t, a, 2) = 0;
      exp_diag(t, a, 3) = 0;
      exp_diag(t, a, 4) = 0;
      exp_diag(t, a, 5) = 0;

      for (int t0 = 0; t0 <= t; ++t0)
      {
        exp_diag(t, a, 0) += diag_arr1[t][a][t0][0];
        exp_diag(t, a, 1) += diag_arr1[t][a][t0][1];
        exp_diag(t, a, 2) += diag_arr1[t][a][t0][2];
        exp_diag(t, a, 3) += diag_arr1[t][a][t0][3];
        exp_diag(t, a, 4) += diag_arr1[t][a][t0][4];
        exp_diag(t, a, 5) += diag_arr1[t][a][t0][5];
      }
    }
  }

  return exp_diag;
}

// [[Rcpp::export]]
arma::cube prev_iter_fct(const arma::vec &H, const arma::cube &d, const arma::mat &q, const arma::mat &prev, int nquar, int nage)
{

  arma::cube exp_diag(nquar, nage, 4);
  const arma::vec empty_vec(5, arma::fill::zeros);
  const arma::vec empty_vec1(4, arma::fill::zeros);
  int nagegrp = d.n_slices;

  std::vector<std::vector<std::vector<arma::mat>>> PA(nquar, std::vector<std::vector<arma::mat>>(nage, std::vector<arma::mat>(nagegrp)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr1(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      for (int ag = 0; ag < nagegrp; ++ag)
      {
        arma::mat d1 = d.slice(ag);
        PA[t][a][ag] = Q1_fct(d1, q, t, a);
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      // Adding initial prevalence at t==0
      arma::vec diag_p(5, arma::fill::zeros);
      arma::vec lat_p(4, arma::fill::zeros);
      arma::vec tmp(4);
      int ag0 = age_to_agroup(a0);

      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;
      diag_arr[t0][a0][t0] = empty_vec;

      if (t0 == 0)
      {
        tmp(0) = prev(0, a0);
        tmp(1) = prev(1, a0);
        tmp(2) = prev(2, a0);
        tmp(3) = prev(3, a0);

        // Progressed initial prevalence
        lat_p = lat_p + (PA[0][a0][ag0].t() * tmp);

        // Diagnosed initial prevalence
        diag_p(0) = tmp(0) * d(0, 0, ag0);
        diag_p(1) = tmp(1) * d(1, 0, ag0);
        diag_p(2) = tmp(2) * d(2, 0, ag0);
        diag_p(3) = tmp(3) * d(3, 0, ag0);
        diag_p(4) = tmp(3) * (1 - d(3, 0, ag0)) * q(3, a0);
      }
      else
      { // Probably dont need this
        diag_p = empty_vec;
        lat_p = empty_vec1;
      }

      lat_arr[t0][a0][t0][0] = H((a0 * nquar) + t0) + lat_p(0);
      lat_arr[t0][a0][t0][1] = lat_p(1);
      lat_arr[t0][a0][t0][2] = lat_p(2);
      lat_arr[t0][a0][t0][3] = lat_p(3);

      diag_arr[t0][a0][t0] = diag_p; // use row_vec instead?

      for (int t = (t0 + 1); t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);

        lat_arr[t][a0][t0] = PA[t][a0][ag].t() * lat_arr[t - 1][a0][t0];

        arma::vec tmpvec(5);

        tmpvec[0] = lat_arr[t - 1][a0][t0][0] * d(0, t, ag);
        tmpvec[1] = lat_arr[t - 1][a0][t0][1] * d(1, t, ag);
        tmpvec[2] = lat_arr[t - 1][a0][t0][2] * d(2, t, ag);
        tmpvec[3] = lat_arr[t - 1][a0][t0][3] * d(3, t, ag);
        tmpvec[4] = lat_arr[t - 1][a0][t0][3] * (1 - d(3, t, ag)) * q(3, a0);

        diag_arr[t][a0][t0] = tmpvec;
      }
    }
  }

  lat_arr1 = a0toa_fct(lat_arr, nquar, nage);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      exp_diag(t, a, 0) = 0;
      exp_diag(t, a, 1) = 0;
      exp_diag(t, a, 2) = 0;
      exp_diag(t, a, 3) = 0;

      for (int t0 = 0; t0 <= t; ++t0)
      {
        exp_diag(t, a, 0) += lat_arr1[t][a][t0][0];
        exp_diag(t, a, 1) += lat_arr1[t][a][t0][1];
        exp_diag(t, a, 2) += lat_arr1[t][a][t0][2];
        exp_diag(t, a, 3) += lat_arr1[t][a][t0][3];
      }
    }
  }

  return exp_diag;
}

// [[Rcpp::export]]
arma::cube prev_iter_fct_rita(const arma::vec &H, const arma::cube &d, const arma::mat &q, const arma::mat &prev, int nquar, int nage)
{
  arma::cube exp_diag(nquar, nage, 7);
  const arma::vec empty_vec(6, arma::fill::zeros);
  const arma::vec empty_vec1(7, arma::fill::zeros);
  int nagegrp = d.n_slices;

  std::vector<std::vector<std::vector<arma::mat>>> PA(nquar, std::vector<std::vector<arma::mat>>(nage, std::vector<arma::mat>(nagegrp)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr1(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));
  std::vector<std::vector<std::vector<arma::vec>>> diag_arr(nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      for (int ag = 0; ag < nagegrp; ++ag)
      {
        arma::mat d1 = d.slice(ag);
        PA[t][a][ag] = Q1_fct_rita(d1, q, t, a);
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      // Adding initial prevalence at t==0
      arma::vec diag_p(6, arma::fill::zeros);
      arma::vec lat_p(7, arma::fill::zeros);
      arma::vec tmp(7);
      int ag0 = age_to_agroup(a0);

      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;
      diag_arr[t0][a0][t0] = empty_vec;

      if (t0 == 0)
      {
        tmp(0) = prev(0, a0);
        tmp(1) = prev(1, a0);
        tmp(2) = prev(2, a0);
        tmp(3) = prev(3, a0);
        tmp(4) = prev(4, a0);
        tmp(5) = prev(5, a0);
        tmp(6) = prev(6, a0);

        // Progressed initial prevalence
        lat_p = lat_p + (PA[0][a0][ag0].t() * tmp);

        // Diagnosed initial prevalence
        diag_p(0) = (tmp(0) * d(0, 0, ag0)) + (tmp(1) * d(0, 0, ag0)) + (tmp(2) * d(0, 0, ag0));
        diag_p(1) = tmp(3) * d(1, 0, ag0);
        diag_p(2) = tmp(4) * d(2, 0, ag0);
        diag_p(3) = tmp(5) * d(3, 0, ag0);
        diag_p(4) = tmp(6) * d(4, 0, ag0);
        diag_p(5) = tmp(6) * (1 - d(4, 0, ag0)) * q(3, a0);
      }
      else
      { // Probably dont need this
        diag_p = empty_vec;
        lat_p = empty_vec1;
      }

      lat_arr[t0][a0][t0][0] = H((a0 * nquar) + t0) + lat_p(0);
      lat_arr[t0][a0][t0][1] = lat_p(1);
      lat_arr[t0][a0][t0][2] = lat_p(2);
      lat_arr[t0][a0][t0][3] = lat_p(3);
      lat_arr[t0][a0][t0][4] = lat_p(4);
      lat_arr[t0][a0][t0][5] = lat_p(5);
      lat_arr[t0][a0][t0][6] = lat_p(6);

      diag_arr[t0][a0][t0] = diag_p; // use row_vec instead?

      for (int t = (t0 + 1); t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);
        lat_arr[t][a0][t0] = PA[t][a0][ag].t() * lat_arr[t - 1][a0][t0];

        arma::vec tmpvec(6);

        tmpvec[0] = (lat_arr[t - 1][a0][t0][0] * d(0, t, ag)) + (lat_arr[t - 1][a0][t0][1] * d(0, t, ag)) + (lat_arr[t - 1][a0][t0][2] * d(0, t, ag));
        tmpvec[1] = lat_arr[t - 1][a0][t0][3] * d(1, t, ag);
        tmpvec[2] = lat_arr[t - 1][a0][t0][4] * d(2, t, ag);
        tmpvec[3] = lat_arr[t - 1][a0][t0][5] * d(3, t, ag);
        tmpvec[4] = lat_arr[t - 1][a0][t0][6] * d(4, t, ag);
        tmpvec[5] = lat_arr[t - 1][a0][t0][6] * (1 - d(4, t, ag)) * q(3, a0);

        diag_arr[t][a0][t0] = tmpvec;
      }
    }
  }

  lat_arr1 = a0toa_fct(lat_arr, nquar, nage);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      exp_diag(t, a, 0) = 0;
      exp_diag(t, a, 1) = 0;
      exp_diag(t, a, 2) = 0;
      exp_diag(t, a, 3) = 0;
      exp_diag(t, a, 4) = 0;
      exp_diag(t, a, 5) = 0;
      exp_diag(t, a, 6) = 0;

      for (int t0 = 0; t0 <= t; ++t0)
      {
        exp_diag(t, a, 0) += lat_arr1[t][a][t0][0];
        exp_diag(t, a, 1) += lat_arr1[t][a][t0][1];
        exp_diag(t, a, 2) += lat_arr1[t][a][t0][2];
        exp_diag(t, a, 3) += lat_arr1[t][a][t0][3];
        exp_diag(t, a, 4) += lat_arr1[t][a][t0][4];
        exp_diag(t, a, 5) += lat_arr1[t][a][t0][5];
        exp_diag(t, a, 6) += lat_arr1[t][a][t0][6];
      }
    }
  }

  return exp_diag;
}

// [[Rcpp::export]]
arma::cube prev_iter_fct_migration(
    const arma::vec &H,               // length nquar*nage (age-blocked: t varies fastest)
    const arma::cube &d,              // dims: 4 (state) x nquar (time) x 4 (agegrp)
    const arma::mat &q,               // 4 x nage
    const arma::mat &prev_total,      // 4 x nage (init_prev_age + init_prev_abroad_age)
    const arma::vec &omega,           // length nquar
    const arma::vec &p,               // length nquar (UK-acquired prop among migrants)
    const arma::mat &m,               // 5 x nquar (migration probs; rows 1..4 are abroad CD4 states)
    const arma::mat &abroad_age_prop, // nquar x nage (rows sum to 1; fallback uniform handled in R)
    int nquar,
    int nage)
{
  arma::cube exp_prev(nquar, nage, 4, arma::fill::zeros);
  const arma::vec empty_vec1(4, arma::fill::zeros);

  int nagegrp = d.n_slices; // should be 4

  // Precompute transition matrices PA[t][a0][ag] where a0 is age-at-entry index
  std::vector<std::vector<std::vector<arma::mat>>> PA(
      nquar, std::vector<std::vector<arma::mat>>(nage, std::vector<arma::mat>(nagegrp)));

  for (int t = 0; t < nquar; ++t)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      for (int ag = 0; ag < nagegrp; ++ag)
      {
        arma::mat d1 = d.slice(ag);           // 4 x nquar
        PA[t][a0][ag] = Q1_fct(d1, q, t, a0); // uses q(.,a0) = progression by age at entry
      }
    }
  }

  // Cohort latent arrays (t, a0, t0) with 4-state vectors
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr(
      nquar, std::vector<std::vector<arma::vec>>(nage, std::vector<arma::vec>(nquar)));

  // Init with zeros
  for (int t = 0; t < nquar; ++t)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      for (int t0 = 0; t0 < nquar; ++t0)
      {
        lat_arr[t][a0][t0] = empty_vec1;
      }
    }
  }

  for (int t0 = 0; t0 < nquar; ++t0)
  {
    for (int a0 = 0; a0 < nage; ++a0)
    {
      arma::vec lat_p(4, arma::fill::zeros);
      int ag0 = age_to_agroup(a0);

      // Add (combined) initial prevalence only at t0==0 (quarter 1)
      if (t0 == 0)
      {
        arma::vec tmp(4);
        tmp(0) = prev_total(0, a0);
        tmp(1) = prev_total(1, a0);
        tmp(2) = prev_total(2, a0);
        tmp(3) = prev_total(3, a0);

        // survivors after diagnosis+progression in quarter 1
        lat_p += (PA[0][a0][ag0].t() * tmp);
      }

      // Inflows at entry time t0 and entry age a0
      double inflow_native = H((a0 * nquar) + t0);
      double inflow_mig_uk = inflow_native * p(t0);

      double inflow_ab1 = omega(t0) * abroad_age_prop(t0, a0) * m(0, t0);
      double inflow_ab2 = omega(t0) * abroad_age_prop(t0, a0) * m(1, t0);
      double inflow_ab3 = omega(t0) * abroad_age_prop(t0, a0) * m(2, t0);
      double inflow_ab4 = omega(t0) * abroad_age_prop(t0, a0) * m(3, t0);

      // Combined (native + migrant) undiagnosed latent state at t0
      lat_arr[t0][a0][t0](0) = inflow_native + inflow_mig_uk + inflow_ab1 + lat_p(0);
      lat_arr[t0][a0][t0](1) = inflow_ab2 + lat_p(1);
      lat_arr[t0][a0][t0](2) = inflow_ab3 + lat_p(2);
      lat_arr[t0][a0][t0](3) = inflow_ab4 + lat_p(3);

      for (int t = (t0 + 1); t < nquar; ++t)
      {
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);

        // propagate undiagnosed (diagnosis probs depend on current age group via ag, progression on a0)
        lat_arr[t][a0][t0] = PA[t][a0][ag].t() * lat_arr[t - 1][a0][t0];
      }
    }
  }

  // Convert cohort-by-entry-age arrays to cohort-by-current-age and aggregate t0
  std::vector<std::vector<std::vector<arma::vec>>> lat_arr1 = a0toa_fct(lat_arr, nquar, nage);

  for (int t = 0; t < nquar; ++t)
  {
    for (int a = 0; a < nage; ++a)
    {
      arma::vec tot(4, arma::fill::zeros);
      for (int t0 = 0; t0 <= t; ++t0)
      {
        tot += lat_arr1[t][a][t0];
      }
      exp_prev(t, a, 0) = tot(0);
      exp_prev(t, a, 1) = tot(1);
      exp_prev(t, a, 2) = tot(2);
      exp_prev(t, a, 3) = tot(3);
    }
  }

  return exp_prev;
}