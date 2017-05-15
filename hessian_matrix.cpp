#include "hessian_matrix.h"

vector<float>& slice::operator [](int index_y){
    return this->gray_scale_[index_y];
}

void volume::make_gaussian_window_(const int size, const float standard_deviation){

    //init
    this->gaussian_window_.resize(size);
    for(int i=0;i<size;++i){
        this->gaussian_window_[i].resize(size);
        for(int j=0;j<size;++j){
            this->gaussian_window_[i][j].resize(size,0.0);
        }
    }

    //sd : standard_deviation
    //g(x,y,z) = N * exp[ -(x^2 + y^2 + z^2)/sd^2 ];
    //where N = 1 / ( sd^3 * (2pi)^(3/2) )

    float maximum = 0.0;
    float summation = 0.0;

    //make N first
    float N = 1.0 / ( pow(standard_deviation,3) * pow( (2.0 * M_PI), 1.5 ) );

    //do g(x,y,z)
    for(int i=0;i<size;++i){
        for(int j=0;j<size;++j){
            for(int k=0;k<size;++k){

                float fi = (float)i - (float)(size-1) / 2.0;
                float fj = (float)j - (float)(size-1) / 2.0;
                float fk = (float)k - (float)(size-1) / 2.0;

                float exp_part = -( fi*fi + fj*fj + fk*fk ) / (standard_deviation*standard_deviation);
                gaussian_window_[i][j][k] = N * exp(exp_part);
                //maximum = max( maximum, gaussian_window[i][j][k] );
                summation += gaussian_window_[i][j][k];
            }
        }
    }

    //normalize the summation of gaussian_window to 1.0
    float ratio = 1.0/summation;
    for(int i=0;i<size;++i){
        for(int j=0;j<size;++j){
            for(int k=0;k<size;++k){
                gaussian_window_[i][j][k] *= ratio;
                maximum = max( maximum, gaussian_window_[i][j][k] );
            }
        }
    }

    return;
}

volume::volume(vector<vector<vector<float> > > v){
    // init
    this->intensity_.resize(v.size());
    for(unsigned int i=0;i<v.size();++i){
        this->intensity_[i] = slice( v[i] );
    }
}

float volume::Ix_(int x, int y, int z){
    if( x+1 >= this->intensity_[z][y].size() )
        return this->intensity_[z][y][x] - this->intensity_[z][y][x-1];
    else if( x-1 < 0)
        return this->intensity_[z][y][x+1] - this->intensity_[z][y][x];
    else
        return ( this->intensity_[z][y][x+1] - this->intensity_[z][y][x-1] ) / 2.0;
}

float volume::Iy_(int x, int y, int z){
    if( y+1 >= this->intensity_[z].size() )
        return this->intensity_[z][y][x] - this->intensity_[z][y-1][x];
    else if( y-1 < 0)
        return this->intensity_[z][y+1][x] - this->intensity_[z][y][x];
    else
        return ( this->intensity_[z][y+1][x] - this->intensity_[z][y-1][x] ) / 2.0;
}

float volume::Iz_(int x, int y, int z){
    if( z+1 >= this->intensity_.size() )
        return this->intensity_[z][y][x] - this->intensity_[z-1][y][x];
    else if( z-1 < 0)
        return this->intensity_[z+1][y][x] - this->intensity_[z][y][x];
    else
        return ( this->intensity_[z+1][y][x] - this->intensity_[z-1][y][x] ) / 2.0;
}

float volume::summation_within_window_gaussianed_(int x, int y, int z, int size){

    float summation = 0.0;

    for(int i=0;i<size;++i){ // x
        for(int j=0;j<size;++j){ // y
            for(int k=0;k<size;++k){ // z

                int sx = x+i;
                int sy = y+j;
                int sz = z+k;

                sz = sz < 0 ? 0 : sz;
                sy = sy < 0 ? 0 : sy;
                sx = sx < 0 ? 0 : sx;

                sz = sz >= this->intensity_.size() ? this->intensity_.size()-1 : sz;
                sy = sy >= this->intensity_[sz].size() ? this->intensity_[sz].size()-1 : sy;
                sx = sx >= this->intensity_[sz][sy].size() ? this->intensity_[sz][sy].size()-1 : sx;

                summation += this->intensity_[sz][sy][sx] * this->gaussian_window_[k][j][i];
            }
        }
    }

    return summation;
}

void volume::make_differential_matrix_(){

    //init
    progressbar *progress = progressbar_new("Initializing",this->intensity_.size());
    differential_matrix_.resize(this->intensity_.size());
    #pragma omp parallel for
    for(int i=0;i<differential_matrix_.size();++i){
        this->differential_matrix_[i].resize(this->intensity_[i].size());
        for(int j=0;j<differential_matrix_[i].size();++j){
            differential_matrix_[i][j].resize(this->intensity_[i][j].size(),matrix(3,0));
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    /* differential_matrix      j->
     * __                           __
     * |    IxIx    IxIy    IxIz     |  i
     * |                             |  |
     * |    IxIy    IyIy    IyIz     |  v
     * |                             |
     * |    IxIz    IyIz    IzIz     |
     * L_                           _|
     */

    progress = progressbar_new("Calculating",this->intensity_.size());
    #pragma omp parallel for
    for(int z=0;z<this->intensity_.size();++z){
        for(int y=0;y<this->intensity_[z].size();++y){
            for(int x=0;x<this->intensity_[z][y].size();++x){

                matrix &this_matrix = differential_matrix_[z][y][x];
                float Ix = this->Ix_(x,y,z);
                float Iy = this->Iy_(x,y,z);
                float Iz = this->Iz_(x,y,z);

                this_matrix[0][0] = Ix*Ix;
                this_matrix[1][1] = Iy*Iy;
                this_matrix[2][2] = Iz*Iz;

                this_matrix[0][1] = this_matrix[1][0] = Ix*Iy;
                this_matrix[0][2] = this_matrix[2][0] = Ix*Iz;
                this_matrix[1][2] = this_matrix[2][1] = Iy*Iz;
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

void volume::make_tensor_(const int window_size){

    //init
    this->tensor_.resize(this->differential_matrix_.size());
    progressbar *progress = progressbar_new("Initializing",this->differential_matrix_.size());
    #pragma omp parallel for
    for(int i=0;i<this->tensor_.size();++i){
        this->tensor_[i].resize(this->differential_matrix_[i].size());
        for(int j=0;j<this->tensor_[i].size();++j){
            this->tensor_[i][j].resize(this->differential_matrix_[i][j].size());
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    // struct tensor A = sum_u_v_w( gaussian(u,v,w) * differential(u,v,w) )

    //for every points
    progress = progressbar_new("Calculating",this->tensor_.size());
    #pragma omp parallel for
    for(int z=0;z<this->tensor_.size();++z){
        for(int y=0;y<this->tensor_[z].size();++y){
            for(int x=0;x<this->tensor_[z][y].size();++x){

                matrix temp(3,0);
                //inside the window
                for(int k=z-window_size/2;k<z+(window_size+1)/2;++k){
                    for(int j=y-window_size/2;j<y+(window_size+1)/2;++j){
                        for(int i=x-window_size/2;i<x+(window_size+1)/2;++i){
                            //check boundary
                            if( k < 0 || k >= this->tensor_.size() ||
                                    j < 0 || j >= this->tensor_[k].size() ||
                                    i < 0 || i >= this->tensor_[k][j].size())
                                continue;
                            //sum it up with gaussian ratio
                            int k_g = k - (z-window_size/2);
                            int j_g = j - (y-window_size/2);
                            int i_g = i - (x-window_size/2);
                            temp += differential_matrix_[k][j][i] * gaussian_window_[k_g][j_g][i_g];
                        }
                    }
                }
                this->tensor_[z][y][x] = temp;
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

void volume::make_eigen_values_(){

    cout << "making eigen values..." <<endl;

    //init
    this->eigen_values_.resize(this->tensor_.size());

    progressbar *progress = progressbar_new("Initializing",this->eigen_values_.size());
    #pragma omp parallel for
    for(int i=0;i<this->eigen_values_.size();++i){
        this->eigen_values_[i].resize( this->tensor_[i].size() );
        for(int j=0;j<this->eigen_values_[i].size();++j){
            this->eigen_values_[i][j].resize( this->tensor_[i][j].size() );
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                this->eigen_values_[i][j][k].resize(3,0.0);
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //using gsl for eigenvalue
    progress = progressbar_new("Calculating",this->tensor_.size());

    #pragma omp parallel for
    for(int i=0;i<this->tensor_.size();++i){
        for(int j=0;j<this->tensor_[i].size();++j){
            for(int k=0;k<this->tensor_[i][j].size();++k){

                matrix &this_matrix = tensor_[i][j][k];

                //allocate needed
                gsl_matrix *tensor_matrix = gsl_matrix_alloc(3,3);
                gsl_vector *eigen_value = gsl_vector_alloc(3);
                gsl_matrix *eigen_vector = gsl_matrix_alloc(3,3);
                gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

                //convert tensor_[i][j][k] to gsl_matrix
                for(int x=0;x<3;++x){
                    for(int y=0;y<3;++y){
                        gsl_matrix_set(tensor_matrix,x,y,this_matrix[x][y]);
                    }
                }

                //do the eigenvalue thing
                gsl_eigen_symmv(tensor_matrix,eigen_value,eigen_vector,w);

                gsl_eigen_symmv_sort(eigen_value,eigen_vector,GSL_EIGEN_SORT_ABS_ASC);

                //save absolute of it to eigen_values_
                for(int x=0;x<3;++x){
                    float ev = gsl_vector_get(eigen_value,x);
                    this->eigen_values_[i][j][k][x] = ev > 0.0 ? ev : -ev;
                }

                //free everything
                gsl_matrix_free(tensor_matrix);
                gsl_matrix_free(eigen_vector);
                gsl_vector_free(eigen_value);
                gsl_eigen_symmv_free(w);
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return ;
}

void volume::make_differential_matrix_(int start_z, int number_z){

    //init
    this->differential_matrix_.resize(this->intensity_.size());

    /* differential_matrix      j->
     * __                           __
     * |    IxIx    IxIy    IxIz     |  i
     * |                             |  |
     * |    IxIy    IyIy    IyIz     |  v
     * |                             |
     * |    IxIz    IyIz    IzIz     |
     * L_                           _|
     */

    for(int z=0;z<this->intensity_.size();++z){

        if( z < start_z || z >= start_z+number_z ){ // clear it because it's not needed
            this->differential_matrix_[z].clear();

        }
        else if(this->differential_matrix_[z].size() == 0){ // only calculate one which not calculated before

            //init
            this->differential_matrix_[z].resize(this->intensity_[z].size());
            #pragma omp parallel for
            for(int y=0;y<this->differential_matrix_[z].size();++y){
                this->differential_matrix_[z][y].resize(this->intensity_[z][y].size(),matrix(3,0.0));
            }

            //calculating
            #pragma omp parallel for
            for(int y=0;y<this->intensity_[z].size();++y){
                for(int x=0;x<this->intensity_[z][y].size();++x){

                    matrix &this_matrix = differential_matrix_[z][y][x];
                    float Ix = this->Ix_(x,y,z);
                    float Iy = this->Iy_(x,y,z);
                    float Iz = this->Iz_(x,y,z);

                    this_matrix[0][0] = Ix*Ix;
                    this_matrix[1][1] = Iy*Iy;
                    this_matrix[2][2] = Iz*Iz;

                    this_matrix[0][1] = this_matrix[1][0] = Ix*Iy;
                    this_matrix[0][2] = this_matrix[2][0] = Ix*Iz;
                    this_matrix[1][2] = this_matrix[2][1] = Iy*Iz;
                }
            }
        }
    }

    return;
}

void volume::make_tensor_(const int window_size, int index_z){
    //init
    this->tensor_.clear();
    this->tensor_.resize(this->intensity_.size());

    // struct tensor A = sum_u_v_w( gaussian(u,v,w) * differential(u,v,w) )

    //init
    this->tensor_[index_z].resize(this->intensity_[index_z].size());
    #pragma omp parallel for
    for(int y=0;y<this->tensor_[index_z].size();++y){
        this->tensor_[index_z][y].resize(this->intensity_[index_z][y].size());
    }

    //calculating
    #pragma omp parallel for
    for(int y=0;y<this->intensity_[index_z].size();++y){
        for(int x=0;x<this->intensity_[index_z][y].size();++x){

            matrix temp(3,0);
            //inside the window
            for(int k=index_z-window_size/2;k<index_z+(window_size+1)/2;++k){
                for(int j=y-window_size/2;j<y+(window_size+1)/2;++j){
                    for(int i=x-window_size/2;i<x+(window_size+1)/2;++i){
                        //check boundary
                        if( k < 0 || k >= this->tensor_.size() ||
                                j < 0 || j >= this->tensor_[k].size() ||
                                i < 0 || i >= this->tensor_[k][j].size())
                            continue;
                        //sum it up with gaussian ratio
                        int k_g = k - (index_z-window_size/2);
                        int j_g = j - (y-window_size/2);
                        int i_g = i - (x-window_size/2);
                        temp += differential_matrix_[k][j][i] * gaussian_window_[k_g][j_g][i_g];
                    }
                }
            }
            this->tensor_[index_z][y][x] = temp;

        }
    }

    return;
}

void volume::eigen_values_initialize_(){

    //init
    this->eigen_values_.resize(this->intensity_.size());

    progressbar *progress = progressbar_new("EigenValueInit",this->eigen_values_.size());
    #pragma omp parallel for
    for(int i=0;i<this->eigen_values_.size();++i){
        this->eigen_values_[i].resize( this->intensity_[i].size() );

        for(int j=0;j<this->eigen_values_[i].size();++j){
            this->eigen_values_[i][j].resize( this->intensity_[i][j].size() );

            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                this->eigen_values_[i][j][k].resize(3,0.0);
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

void volume::make_eigen_values_(int index_z){

    if(this->intensity_.size() >= TIFF_IMAGE_LARGE_SIZE){ // for the super large data
        //allocate the needed and free others
        this->eigen_values_.clear();
        this->eigen_values_.resize( this->intensity_.size() );
        this->eigen_values_[index_z].resize( this->intensity_[index_z].size() );

        #pragma omp for
        for(int j=0;j<this->eigen_values_[index_z].size();++j){
            this->eigen_values_[index_z][j].resize( this->intensity_[index_z][j].size() );

            for(int k=0;k<this->eigen_values_[index_z][j].size();++k){
                this->eigen_values_[index_z][j][k].resize(3,0.0);
            }
        }
    }

    //using gsl for eigenvalue
    #pragma omp parallel for
    for(int j=0;j<this->tensor_[index_z].size();++j){
        for(int k=0;k<this->tensor_[index_z][j].size();++k){

            matrix &this_matrix = tensor_[index_z][j][k];

            //allocate needed
            gsl_matrix *tensor_matrix = gsl_matrix_alloc(3,3);
            gsl_vector *eigen_value = gsl_vector_alloc(3);
            gsl_matrix *eigen_vector = gsl_matrix_alloc(3,3);
            gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

            //convert tensor_[i][j][k] to gsl_matrix
            for(int x=0;x<3;++x){
                for(int y=0;y<3;++y){
                    gsl_matrix_set(tensor_matrix,x,y,this_matrix[x][y]);
                }
            }

            //do the eigenvalue thing
            gsl_eigen_symmv(tensor_matrix,eigen_value,eigen_vector,w);

            gsl_eigen_symmv_sort(eigen_value,eigen_vector,GSL_EIGEN_SORT_ABS_ASC);

            //save absolute of it to eigen_values_
            for(int x=0;x<3;++x){
                float ev = gsl_vector_get(eigen_value,x);
                this->eigen_values_[index_z][j][k][x] = ev > 0.0 ? ev : -ev;
            }

            //free everything
            gsl_matrix_free(tensor_matrix);
            gsl_matrix_free(eigen_vector);
            gsl_vector_free(eigen_value);
            gsl_eigen_symmv_free(w);
        }
    }

    return;
}

void volume::calculate_hessian(const int window_size, float threshold, const float standard_deviation){

    cout << "making gaussian window with window_size : " << window_size;
    (cout << "\tstandard_deviation : " << standard_deviation ).flush();

    this->make_gaussian_window_(window_size,standard_deviation*(float)window_size/2.0);
    cout << "\tdone!"<<endl;

    //init
    this->eigen_values_initialize_();

    //load data when needed, free it otherwise
    progressbar *progress = progressbar_new("Calculating",this->intensity_.size());
    for(int i=0;i<this->intensity_.size();++i){

        int number_z = window_size;
        int start_z = (i - window_size/2) >= 0 ? (i - window_size/2) : 0 ;
        start_z = (start_z+number_z) <= this->intensity_.size() ? start_z  : this->intensity_.size() - number_z;

        this->make_differential_matrix_(start_z, number_z);
        this->make_tensor_(window_size, i);
        this->make_eigen_values_(i);

        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

vector<float> operator -(vector<float> &a, vector<float> &b){
    if(a.size() != b.size()){
        cout << "ERROR : vector size not compatible" <<endl;
        exit(-1);
    }

    vector<float> result(a.size(),0.0);

    for(int i=0;i<a.size();++i){
        result[i] = a[i] - b[i];
    }

    return result;

}

vector<float> operator +(vector<float> &a, vector<float> &b){
    if(a.size() != b.size()){
        cout << "ERROR : vector size not compatible" <<endl;
        exit(-1);
    }

    vector<float> result(a.size(),0.0);

    for(int i=0;i<a.size();++i){
        result[i] = a[i] + b[i];
    }

    return result;
}

float vector_dot(vector<float> &a, vector<float> &b){
    if(a.size() != b.size()){
        cout << "ERROR : vector size not compatible" <<endl;
        exit(-1);
    }

    float result = 0.0;

    for(int i=0;i<a.size();++i){
        result += a[i] * b[i];
    }

    return result;

}

float vector_length(vector<float> &a){

    float result = 0.0;

    for(int i=0;i<a.size();++i){
        result += a[i] * a[i];
    }

    result = pow(result, 0.5);

    return result;
}
