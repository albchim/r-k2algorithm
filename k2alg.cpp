#include<cstdlib>
#include<iostream>
#include<vector>
#include<cmath>
#include<set>
#include<algorithm>
using namespace std;
//linear time in N=\prod_l r[p_i[l]]+pLen
//we dot it N*n times

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

/* 
int* compute_instances(int i, int* r, int pLen, int* p_i, int& nrows)
{
    int ncols = pLen; //the number of possible instances of x_i
    nrows = 1;
    for(int l = 0; l < pLen; l++)
    {
        nrows *= r[p_i[l]];
    }
    int* instances = (int*)malloc(nrows*pLen*sizeof(int));
    int c = 1;
    for(int j = 0; j <nrows; j++)
    {
        
        //we have to make a correspondance between
        //the instantiation and the actual values taken
        //we have to resolve the enumeration
        int c = r[p_i[pLen-1]] ;
        for(int l = pLen-1; l >=0; l--)
        {
            int val;
           if(l==pLen-1)
           {
               val = (j%c);
           }
           else
           {
               val = j/c;
               c*= r[p_i[l]];
           }
            instances[l+ncols*j] = val%r[p_i[l]];
        }
    }
    return instances;
}
*/
int compute_instance_index(int* row, int* r, int pLen, vector<int>& p_i)
{
    int j = 0;
    int c = r[p_i[pLen-1]];
    for(int l = pLen-1; l >=0; l--)
    {
        if(l==pLen-1)
        {
            j += row[p_i[l]];
        }
        else
        {
            j += c*row[p_i[l]];
            c *= r[p_i[l]];
        }
    }
    return j;
}

int* compute_alpha(int i, int n, int* r, int pLen, vector<int>& p_i, int m, int* data, int& nrows)
{

    nrows = 1;
    for(int l = 0; l < pLen; l++)
    {
        nrows *= r[p_i[l]];
    }
    int r_i = r[i];
    int* alpha = (int*)malloc(nrows*r_i*sizeof(int));
    //instantiate alpha to zero
    for(int a = 0; a < nrows*r_i; a++)
    {
        alpha[a] = 0;
    }
  
    for(int a = 0; a < m; a++)
    {
        int k = data[i+a*n];
        int j = compute_instance_index(data+a*n, r, pLen, p_i);
        alpha[k+j*r_i]++;
    }
    
    return alpha;
}

int* compute_alpha_nop(int i, int n, int* r, int m, int* data)
{
    int r_i = r[i];
    int* alpha = (int*) malloc(r_i*sizeof(int));
    //instantiate alpha to zero
    for(int a = 0; a < r_i; a++)
    {
        alpha[a] = 0;
    }
    for(int k = 0; k < r_i; k++)
    {
        for(int a = 0; a < m; a++)
        {
            if(data[i+n*a]==k)
            {
                alpha[k]++;
            }
        }
    }
    return alpha;
}

double compute_f_nop(int i, int n, int* r,  int m, int* data)
{
    int* alpha = compute_alpha_nop(i, n, r, m, data);
    int N = 0;
    double prod = 1;
    for(int k = 0; k < r[i]; k++)
    {
        N += alpha[k];
        prod *= factorial(alpha[k]);
    }
    free(alpha);
    return factorial(r[i]-1)*prod/factorial(N+r[i]-1);
}

double compute_f(int i, int n, int* r, int pLen, vector<int> p_i, int m, int* data)
{
    int nrows = 0;
    int* alpha = compute_alpha(i,n,r,pLen, p_i, m, data, nrows);
    int* N = (int*)malloc(nrows*sizeof(int));
    int r_i = r[i];
    for(int j = 0; j < nrows; j++)
    {
        N[j]=0;
        for(int k = 0; k < r_i; k++)
        {
            N[j] += alpha[k+j*r_i];
        }
    }
    double f = 1;
    for(int j = 0; j < nrows; j++)
    {
        f *= factorial(r_i-1);
        f /=  factorial(N[j]+r_i-1);
        for(int k = 0; k < r_i; k++)
        {
            f *= factorial(alpha[k+j*r_i]);
        }
    }
    free(alpha);
    return f;
}

void k2alg(vector<int>& cp, int u,int i, int n, int* r, int m, int* data, int& pLen, vector<int>& p)
{
    p = vector<int>();
    p.reserve(u);
    pLen = 0;
    double pOld = compute_f_nop(i, n, r, m, data);
    bool flag = true;
    while(flag && pLen < u && cp.size() > 0)
    {
        cout << "Old: " << pOld << endl;
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
            p.pop_back();
        }
        if(cMax==-1)
        {
            flag = false;
        }
        else
        {
            cout << "Candidate: " << cMax << " with prob. " << pMax << endl;
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
}


vector<vector<int>> k2procedure(int u, int n, int* r, int*order,int m, int* data)
{
    vector<vector<int>> result;
    result.reserve(n);
    for(int i = 0; i < n; i++)
    {
        vector<int> cp;
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
            int pLen;
            vector<int> p_i;
            k2alg(cp,u,i,n,r,m,data, pLen, p_i);
            result.push_back(p_i);
        }
    }
    return result;
}



int main()
{
    int u = 2;
    int r[3] = {2,2,2};
    int n = 3;
    int data[30] = {1,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0};
    int m = 10;
    int order[3] = {0,1,2};
    auto result = k2procedure(u, n, r, order, m, data);
}

/*

order <- c(0,1,2)
r <- c(2,2,2)
data <- matrix(c(1,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0), 10,3)
k2procedure(2,3,r,order,10,data)
 */