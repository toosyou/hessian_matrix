#include <iostream>
#include <vector>
#include <time.h>
#include <cstdio>
#include <cstdlib>

#include "hessian_matrix.h"

using namespace std;

int main(void){

    srand(time(NULL));

    // init testing data
    vector<vector<vector<float> > > v;
    v.resize(20);
    for(int i=0;i<20;++i){
        v[i].resize(20);
        for(int j=0;j<20;++j){
            v[i][j].resize(20);
            for(int k=0;k<20;++k){
                v[i][j][k] = ((float)rand() / (float)RAND_MAX); // 0.0 - 1.0
            }
        }
    }

    volume lung(v);
    lung.calculate_hessian(3); // window size = 3
    // Each (x, y, z) has a vector<float>(e1, e2, e3) for sorted eigen value,
    // that is, e1 < e2 < e3.
    // ( Corresponding eigen vectors v1, v2, v3 are all normalized with eigen values,
    //   so |v1| = |v2| = |v3| = 1.0 )
    vector<vector<vector<vector<float> > > > eigen_values = lung.get_eigen_value();

    for(int i=0;i<20;++i){
        for(int j=0;j<20;++j){
            for(int k=0;k<20;++k){
                for(int m=0;m<3;++m)
                    cout << eigen_values[i][j][k][m] << '\t';
                cout << endl;
            }
        }
    }

    return 0;
}
