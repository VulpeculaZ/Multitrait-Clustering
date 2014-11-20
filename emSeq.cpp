#include <Rcpp.h>
#include <algorithm>

Rcpp::List seq(data);
// copy the input data to Rcpp or armadillo structures
arma::mat contrast_a = Rcpp::as<arma::mat>(contrast);
Rcpp::NumericMatrix ecp(ecps);
Rcpp::IntegerVector weight(weightS);
int ncluster = ecp.ncol();
int nr = Rcpp::as<int>(nrs); // number of different site patterns.
int n = seq.size(); //number of sequences
int nc = Rcpp::as<int>(ncs); //number of distinct codes
int nco  = Rcpp::as<int>(ncos); // number of rows of contrast matrix
int nt = Rcpp::as<int>(nts); // total number of characters in the a sequence
double MachineEpsC = Rcpp::as<double>(MachineEps);
// Private varibles for the function
arma::mat nbi(nr, nc);
arma::mat tmpCo(nr, nc);
// Output variables
Rcpp::NumericMatrix logp(Rcpp::clone(ecp));
Rcpp::NumericVector theta(ncluster);
Rcpp::IntegerMatrix y(nr, ncluster);
// Estimate parameter y and theta
for(int l=0; l<ncluster; l++){
    nbi.zeros();
    for(int i=0; i<n; i++){
        Rcpp::IntegerVector tmpSeq = Rcpp::as<Rcpp::IntegerVector>(seq[i]);
        for(int k=0; k<nc; k++){
            for(int j=0; j<nr; j++){
                tmpCo[k * nr + j] = contrast_a[k*nco + tmpSeq[j] -1L];
            }
        }
        nbi = nbi + (tmpCo * ecp[i+l*n]);
    }
    double Ns=0.0;
    for(int j=0; j<nr; j++){
        double nbiRowj[4];
        for(int k=0; k<4; k++){
            nbiRowj[k] =nbi(j,k);
        }
        int Nrow = sizeof(nbiRowj) / sizeof(double);
        y(j,l) = std::distance(nbiRowj, std::max_element(nbiRowj, nbiRowj+Nrow)) + 1L;
        Ns = Ns + nbi(j,y(j,l)-1L) * weight[j];
    }
    Rcpp::NumericVector zl = ecp(_, l);
    double sumZ  = std::accumulate(zl.begin(), zl.end(), 0.0);
    theta[l] = 1.0/3.0 - (Ns/sumZ) / 3.0 / nt;
    if(theta[l] <= 0){
        theta[l] = MachineEpsC;
    }
 }
// Estimate the log of conditional expectation Ez
for(int l=0; l<ncluster; l++){
    for(int i=0; i<n; i++){
        Rcpp::IntegerVector tmpSeq = Rcpp::as<Rcpp::IntegerVector>(seq[i]);
        int nSame = 0;
        for(int j=0; j<nr; j++){
            if(y(j,l)==tmpSeq[j])
                nSame = nSame + weight[j];
        }
        logp(i,l) = nSame * log(1.0 - 3.0 * theta[l]) +  (nt - nSame) * log(theta[l]);
    }
 }

return Rcpp::List::create(Rcpp::Named("y", y),
                          Rcpp::Named("theta", theta),
                          Rcpp::Named("ecps", logp));
