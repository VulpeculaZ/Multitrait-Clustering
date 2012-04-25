#include <Rcpp.h>

Rcpp::List seq(data);
Rcpp::List cvSeq(cvData);
// copy the data to armadillo structures
arma::mat P_a = Rcpp::as<arma::mat>(P);
arma::mat contrast_a = Rcpp::as<arma::mat>(contrast);
arma::colvec bf = Rcpp::as<arma::colvec>(bfs);
Rcpp::NumericMatrix ecp(ecps);
int ncluster = ecp.ncol();
int nr = Rcpp::as<int>(nrs);
int n = seq.size();
int cvn = cvSeq.size();
int nc = Rcpp::as<int>(ncs);
int nco  = Rcpp::as<int>(ncos);
Rcpp::NumericMatrix ecn(cvn, ncluster);
arma::mat lookup = arma::log(contrast_a * P_a);
arma::mat post(nr, nc);
arma::mat tmpPost = arma::zeros<arma::mat>(nr,nc);

for(int l=0; l<ncluster; l++){
    post.zeros();
    for(int i=0; i<n; i++){
        Rcpp::IntegerVector tmpSeq = Rcpp::as<Rcpp::IntegerVector>(seq[i]);
        for(int k=0; k<nc; k++){
            for(int j=0; j<nr; j++){
                tmpPost[k * nr + j] = lookup[k*nco + tmpSeq[j] -1L];
            }
        }
        post = post + (tmpPost * ecp[i+l*n]);
    }
    post = arma::exp(post);
    arma::colvec post_rowsum = arma::sum(post,1);
    for(int k=0; k<nc; k++){
        for(int j=0; j<nr; j++){
            post[k * nr + j] = post[k * nr + j] / post_rowsum[j];
        }
    }
    for(int i=0; i<cvn; i++){
        Rcpp::IntegerVector tmpSeq = Rcpp::as<Rcpp::IntegerVector>(cvSeq[i]);
        for(int k=0; k<nc; k++){
            for(int j=0; j<nr; j++){
                tmpPost[k * nr + j] = lookup[k*nco + tmpSeq[j] -1L];
            }
        }
        ecn[i+l*cvn] = arma::accu(arma::log((post % arma::exp(tmpPost)) * bf));
    }
}

return Rcpp::wrap(ecn);

