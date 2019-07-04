#include<cmath>
#include<algorithm>
#include<Rcpp.h>

using namespace std;
using namespace Rcpp;

double factorial(double n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double approxLogFact(double n)
{
  if(n < 150)
    return log(factorial(n));
  return (n+0.5)*n-n+0.5*log(2*M_PI);
}

int compute_instance_index(int n,int rowIndex, IntegerMatrix data, IntegerVector r, int pLen, IntegerVector p_i)
{
  int j = 0;
  int c = r[p_i[pLen-1]];
  for(int l = pLen-1; l >=0; l--)
  {
    if(l==pLen-1)
    {
      j += data[p_i[l]+n*rowIndex];
    }
    else
    {
      j += c*data[p_i[l]+n*rowIndex];
      c *= r[p_i[l]];
    }
  }
  return j;
}

IntegerMatrix compute_alpha(int i, int n, IntegerVector r, int pLen, IntegerVector p_i, int m, IntegerMatrix data, int& nrows)
{
  nrows = 1;
  for(int l = 0; l < pLen; l++)
  {
    nrows *= r[p_i[l]];
  }
  int r_i = r[i];
  IntegerMatrix alpha(nrows, r_i);
  for(int a = 0; a < m; a++)
  {
    int k = data[i+a*n];
    int j = compute_instance_index(n,a,data, r, pLen, p_i);
    alpha[k+j*r_i]++;
  }
  return alpha;
}

IntegerMatrix compute_alpha_nop(int i, int n, IntegerVector r, int m, IntegerMatrix data)
{
  int r_i = r[i];
  IntegerMatrix alpha(r_i,1);
  for(int a = 0; a < m; a++)
  {
    alpha[data[i+n*a]]++;
  }
  return alpha;
}

double compute_f_nop(int i, int n, IntegerVector r,  int m, IntegerMatrix data)
{
  IntegerMatrix alpha = compute_alpha_nop(i, n, r, m, data);
  int N = 0;
  double prod = 1;
  for(int k = 0; k < r[i]; k++)
  {
    N += alpha[k];
    prod *= factorial(alpha[k]);
  }
  return factorial(r[i]-1)*prod/factorial(N+r[i]-1);
}
const double exp_scale = 100000;
double compute_f(int i, int n, IntegerVector r, int pLen, IntegerVector p_i, int m, IntegerMatrix data)
{
  int nrows = 0;
  IntegerMatrix alpha = compute_alpha(i,n,r,pLen, p_i, m, data, nrows);
  IntegerVector N(nrows);
  int r_i = r[i];
  for(int j = 0; j < nrows; j++)
  {
    for(int k = 0; k < r_i; k++)
    {
      N[j] += alpha[k+j*r_i];
    }
  }
  double f = 1;
  //double lf = nrows*approxLogFact(r_i-1);
  for(int j = 0; j < nrows; j++)
  {
    f *= factorial(r_i-1);
    
    f /=  factorial(N[j]+r_i-1);
    //lf -= approxLogFact(N[j]+r_i-1);
    for(int k = 0; k < r_i; k++)
    {
      f *= factorial(alpha[k+j*r_i]);
      //lf += approxLogFact(alpha[k+j*r_i]);
    }
  }
  //Rcout << pow(exp(lf/100000),100000) << endl;
  //return pow(exp(lf/100000),100000);
  return f;
}

IntegerVector k2alg(IntegerVector& cp, int u,int i, int n, IntegerVector r, int m, IntegerMatrix data)
{
  IntegerVector p;
  int pLen = 0;
  double pOld = compute_f_nop(i, n, r, m, data);
  bool flag = true;
  while(flag && pLen < u && cp.size() > 0)
  {
    int cMax = -1;
    double pMax = 0;
    for(int j = 0; j<cp.size(); j++)
    {
      int candidate = cp[j];
      p.push_back(candidate);
      double pNew = compute_f(i, n, r, pLen+1, p, m, data);
      if(pNew > pMax)
      {
        cMax = candidate;
        pMax = pNew;
      }
      p.erase(pLen);
    }
    if(cMax==-1)
    {
      flag = false;
    }
    else
    {
      if(pMax > pOld)
      {
        p.push_back(cMax);
        pLen++;
        std::remove(cp.begin(), cp.end(), cMax);
        pOld = pMax;
      }
      else
      {
        flag=false;
      }
      
    }        
  }
  return p;
}


bool checkDataset(IntegerMatrix data, IntegerVector r)
{

  int m = data.rows();
  int n = data.cols();
  if(r.size()!=n)
  {
    return false;
  }
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m; j++)
    {
      if(data(j,i)>=r[i])
      {
        return false;
      }
    }
  }
  return true;
}


// [[Rcpp::export]]
List k2procedure(SEXP x,SEXP dims, SEXP varOrder, int u = -1)
{
  IntegerMatrix data(x);
  IntegerVector order(varOrder);
  IntegerVector r(dims);
  bool valid = checkDataset(data, r);
  if(!valid)
  {
    Rcout << "Invalid dataset" << endl;
    return NULL;
  }
  List result;
  int n = data.cols();
  int m = data.rows();
  if(u==-1)
  {
    u=n-1;
  }
  for(int i = 0; i < n; i++)
  {
    IntegerVector cp;
    for(int j = 0; j < n; j++)
    {
      if(order[j]!=i)
      {
        cp.push_back(order[j]);
      }
      else
      {
        break;
      }
    }
    if(cp.size()!=0)
    {
      IntegerVector p_i = k2alg(cp,u,i,n,r,m,data);
      result.push_back(p_i);
    }
    else
    {
      result.push_back(-1);
    }
  }
  return result;
}


/*** R
k2procedureSplit <- function(x, dims, initVarOder, u = -1, subSize = 500)
{
  nRows <- dim(x)[0]
  nSplits <- nRows / subSize
  vOrder <- initVarOder
  for(i in 1:nSplits)
  {
    lastRow <- i + subSize - 1
    if(i==nSplits)
    {
      lastRow <- nRows
    }
    subData <- x[i:lastRow,]
    res <- k2procedure(subData, dims, vOrder, u)
    
  }
}
*/



/*
 //test using these commands
 order <- c(0,1,2)  
 r <- c(2,2,2)
   data <- matrix(c(1,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0), 10,3)
 k2procedure(data, r, order, 2)
 */