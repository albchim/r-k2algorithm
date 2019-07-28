#include<memory>
template<typename T>
struct anode
{
    T x;
    int k;
    anode** sib;
};

template<typename T>
class adjnet
{

private:
    anode* v;
    int n;
    bool fixed;

    //node initialization; assuming maximum size of sibilings
    void init_nodes()
    {
        v = (anode*)malloc(sizeof(anode)*n);
        for(int i = 0; i < n; i++)
        {
            v[i].k = 0;
            v[i].sib = (anode**)malloc(sizeof(anode*)*n);
        }
        fixed = false;
    }

    //node initialization; fixed degree
    void init_nodes(int* degrees, int** sibilings)
    {
        v = (anode*)malloc(sizeof(anode)*n);
        for(int i = 0; i < n; i++)
        {
            v[i] = anode();
            v[i].k = degrees[i];
            if(v[i].k>0)
            {
                v[i].sib = (anode**)malloc(sizeof(anode*)*v[i].k);
            }
        }
        fixed = true;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < v[i].k; j++)
            {
                v[i].sib[j] = v + sibilings[i][j];
            }
        }
    }

public:
    adjnet(int nNodes)
    {
        n = nNodes;
        init_nodes();
    }
    adjnet(int nNodes, int* degs, int** sibs)
    {
        n = nNodes;
        init_nodes(degs, sibs);
    }

    ~adjnet()
    {
        for(int i = 0; i < n; i++)
        {
            free(v[i].sib);
        }
        free(v);
    }


    bool add_edge(int from, int to, bool bothWays = false)
    {
        if(fixed)
        {
            return false;
        }
        v[from].sib[to] = (v+to);
        if(bothWays)
        {
            v[to].sib[from] = (v+from);
        }
        return true;
    }
    bool remove_edge(int from, int to, bool bothWays = false)
    {
        if(fixed)
        {
            return false;
        }
        v[from].sib[to] = NULL;
        if(bothWays)
        {
            v[to].sib[from] = NULL;
        }
        return true;
    }

    T& operator[](const int index)
    {
        return v[i].x;
    }
    const T& operator[](const int index) const
    {
        return v[i].x;
    } 
};
