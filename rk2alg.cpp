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
      j += data(rowIndex, p_i[l]);
    }
    else
    {
      j += c*data(rowIndex, p_i[l]);
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
    int k = data(a,i);
    int j = compute_instance_index(n,a,data, r, pLen, p_i);
    alpha(j,k)++;
  }
  return alpha;
}

IntegerVector compute_alpha_nop(int i, int n, IntegerVector r, int m, IntegerMatrix data)
{
  int r_i = r[i];
  IntegerVector alpha(r_i);
  for(int a = 0; a < m; a++)
  {
    alpha[data(a,i)]++;
  }
  return alpha;
}

double compute_f_nop(int i, int n, IntegerVector r,  int m, IntegerMatrix data)
{
  IntegerVector alpha = compute_alpha_nop(i, n, r, m, data);
  int N = 0;
  double prod = 1;
  for(int k = 0; k < r[i]; k++)
  {
    N += alpha[k];
    prod *= factorial(alpha[k]);
  }
  return factorial(r[i]-1)*prod/factorial(N+r[i]-1);
}

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
      N[j] += alpha(j,k);
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
      f *= factorial(alpha(j,k));
      //lf += approxLogFact(alpha[k+j*r_i]);
    }
  }
  //Rcout << pow(exp(lf/100000),100000) << endl;
  //return pow(exp(lf/100000),100000);
  return f;
}

IntegerVector k2alg(IntegerVector& cp, double& score,int u,int i, int n, IntegerVector r, int m, IntegerMatrix data)
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
  score = pOld;
  if(pLen==0)
  {
    p.push_back(i);
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


List k2procedureInternal(SEXP x,SEXP dims, SEXP varOrder, NumericVector& scores, int u = -1)
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
  scores = NumericVector(n);
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
      double score = 0;
      IntegerVector p_i = k2alg(cp,score,u,i,n,r,m,data);
      scores.push_back(score);
      std::sort(p_i.begin(), p_i.end());
      result.push_back(p_i);
    }
    else
    {
      result.push_back(rep(i,1));
    }
  }
  return result;
}

// [[Rcpp::export]]
SEXP k2procedure(SEXP x,SEXP dims, SEXP varOrder, int u =-1,int returnType = 0, int verbose = 0, int splitSize=100)
{
  IntegerMatrix data(x);
  IntegerVector order(varOrder);
  IntegerVector r(dims);
  int nRows = data.nrow();
  int nCols = data.ncol();
  int nSplits = nRows / splitSize;
  if(nSplits == 0)
  {
    nSplits++;
    splitSize=nRows;
  }
  
  NumericMatrix wadj(nCols, nCols);
  double normScore = 0;
  double bestScore = 0;
  int bestRichness = 0;
  List bestList;
  int from = 0;
  for(int s = 0; s < nSplits; s++)
  {
    int sz = splitSize;
    if(s==nSplits-1 && nRows % splitSize != 0)
    {
      sz = nRows % splitSize;
    }
    //now copy the data matrix
    IntegerMatrix sData(sz, nCols);
    for(int i = 0; i < sz*nCols; i++)
    {
      sData[i] = data[from+i];
    }
    from += sz*nCols;
    
    NumericVector scs(nCols);
    List res =  k2procedureInternal(sData, r, order, scs, u);
    
    double totalSc = mean(scs);
    int richness = 0;
    for(int i = 0; i < nCols; i++)
    {
      if(res[i]!=nullptr)
      {
        IntegerVector temp(res[i]);
        richness += temp.size();
      }
    }
    if(verbose > 0 && s % verbose == 0)
    {
      Rcout << "Split n° " << (s+1) << " with size " << sz << " with score " << totalSc << endl;
      //Rcpp:print(res);
    }
    if(totalSc+richness>bestScore+bestRichness)
    {
      bestScore = totalSc;
      bestRichness = richness;
      bestList = res;
    }
    if(returnType == 2)
    {
      normScore += totalSc;
      for(int i = 0; i < nCols; i++)
      {
        IntegerVector p_i(res.at(i));
        for(int j = 0; j < p_i.size(); j++)
        {
          if(p_i[j]!=i)
            wadj(p_i[j],i)+=normScore;
        }
      }
    }
  }
  
  CharacterVector names(nCols);
  for(int i = 0; i < nCols; i++)
  {
    names[i]=("x"+to_string(i));
  }
  if(returnType==1)
  {
    IntegerMatrix adj(nCols,nCols);
    
    for(int i = 0; i < nCols; i++)
    {
      IntegerVector p_i(bestList.at(i));
      for(int j = 0; j < p_i.size(); j++)
      {
        if(p_i[j]!=i)
          adj(p_i[j],i)=1;
      }
    }
    rownames(adj)=names;
    colnames(adj) =names;
    return adj;
  }
  else if(returnType==2)
  {
    for(int i = 0; i < nCols; i++)
    {
      for(int j = 0; j < nCols; j++)
      {
        wadj[i*nCols+j] /= (normScore);
      }
    }
    rownames(wadj)=names;
    colnames(wadj) =names;
    return wadj;
  }
  else
  {
    return bestList;
  }
}


// [[Rcpp::export]]
NumericMatrix condProb(int i, IntegerMatrix x,IntegerVector dims, IntegerVector parents, int method = 0)
{
  int nRows = 1;
  int m = x.nrow();
  int n= x.ncol();
  int pLen = parents.size();
  //Rcout << "n° samples: "<< m << endl;
  //Rcout << "n° variables " << n << endl;
  if(i < 0 || i >= n)
  {
    return NumericMatrix(0,0);
  }
  for(int l = 0; l < pLen; l++)
  {
    if(parents[l]<0 || parents[l]>= n || parents[l]==i)
    {
      return NumericMatrix(0,0);
    }
  }
  
  if(pLen >= n)
  {
    return NumericMatrix(0,0);
  }
  int r_i = dims[i];
  if(pLen==0)
  {
    return NumericMatrix(0,0);
  }
  for(int l = 0; l < pLen; l++)
  {
    nRows *= dims[parents[l]];
  }
  NumericVector P(nRows*r_i);
  
  if(method == 0)
  {
    
  
  
  
    IntegerVector counts(nRows);
    for(int a = 0; a < m; a++)
    {
      int j = compute_instance_index(n, a, x, dims, pLen, parents); //the actual index in P is different
      P[j*r_i+x(a,i)]++;
      counts[j]++;
    }
    for(int j = 0; j < nRows; j++)
    {
      for(int l = 0; l < r_i; l++)
      {
        P[j*r_i+l] /= counts[j];
      }
    }
  }
  else
  {
    IntegerMatrix Nijk = compute_alpha(i,n,dims,pLen,parents,m,x,nRows);
    IntegerVector Nij(nRows);
    for(int j = 0; j < nRows; j++)
    {
      for(int k = 0; k < r_i; k++)
      {
        Nij(j) += Nijk[j*r_i+k];
        P[j*r_i+k] = Nijk[j*r_i+k]+1;
      }
    }
    for(int j = 0; j < nRows; j++)
    {
      for(int k = 0; k < r_i; k++)
      {
        P[j*r_i+k] /= Nij[j]+r_i;
      }
    }
  }
  NumericMatrix T(nRows*r_i, 2+pLen);
  CharacterVector names(2+pLen);
  names[pLen] = "x" + to_string(i);
  names[pLen+1] = "P";
  for(int l = 0; l < pLen; l++)
  {
    names[l] = "cond. x" + to_string(parents[l]);
  }
  colnames(T) = names;
  for(int j = 0; j < nRows; j++)
  {
    for(int l = 0; l < r_i; l++)
    {
      T ((j*r_i+l),pLen+1) =P[j*r_i+l];
      T (j*r_i+l,pLen)=l;
    }
    int c = nRows;
    for(int l = 0; l < pLen; l++)
    {
      c /= dims[parents[l]];
      int val;

      if(l==pLen-1)
      {
        val = (j % dims[parents[l]]);
      }
      else
      {
        val = (j / c)%dims[parents[l]]; 
      }
      
      for(int k = 0; k < r_i; k++)
      {
        T(j*r_i+k,l)=val;
      }
    }
  }
  return T;
}




/*
 //test using these commands
 order <- c(0,1,2)  
 r <- c(2,2,2)
   data <- matrix(c(1,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0), 10,3)
 k2procedure(data, r, order, 2)
 */