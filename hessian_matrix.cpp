#include "tomo_tiff.h"

tomo_tiff::tomo_tiff(const char* address){
    TIFF *tif = TIFFOpen( address, "r" );
    if(tif == NULL){
        cerr << "ERROR : cannot open " << address << endl;
        return;
    }

    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &this->height_);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,  &this->width_);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &this->bits_per_sample_);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &this->samples_per_pixel_);

    //init
    this->address_ = string(address);
    this->gray_scale_.resize(this->height_);
    for(unsigned int i=0;i<this->height_;++i){
        this->gray_scale_[i].resize(this->width_,0.0);
    }

    //read raw-data
    unsigned int line_size = TIFFScanlineSize(tif);
    char* buf = new char[ line_size * this->height_];
    for(unsigned int i=0;i<this->height_;++i){
        TIFFReadScanline( tif, &buf[i*line_size], i );
    }

    if(this->bits_per_sample_ == 16 && this->samples_per_pixel_ == 1){
        for(unsigned int i=0;i<this->height_;++i){
            for(unsigned int j=0;j<this->width_;++j){
                this->gray_scale_[i][j] = (float)((uint16_t*)buf)[ i*this->width_ + j ] / 65535.0;
            }
        }
    }
    else{
        cerr << "ERROR : " << address << " not handled!" <<endl;
        cerr << "bits_per_sample : " << this->bits_per_sample_ << " ";
        cerr << "samples_per_pixel : " << this->samples_per_pixel_ <<endl;
    }

    delete [] buf;
    TIFFClose(tif);

    return;
}

void tomo_tiff::save(const char* address, int max_gray_scale ){
    TIFF *tif = TIFFOpen(address, "w");
    if(tif == NULL){
        cerr << "ERROR : cannot create file " << address <<endl;
        return;
    }

    this->height_ = this->gray_scale_.size();
    this->width_ = this->gray_scale_[0].size();

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, this->width_);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, this->height_);

    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, this->bits_per_sample_);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, this->samples_per_pixel_);

    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    TIFFSetField(tif, TIFFTAG_XRESOLUTION, 0);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, 0);
    TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);

    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    if(this->bits_per_sample_ == 16 && this->samples_per_pixel_ == 1){
        vector<uint16_t> data;
        for(unsigned int i=0;i<this->height_;++i){
            for(unsigned int j=0;j<this->width_;++j){
                data.push_back( gray_scale_[i][j] * (float)max_gray_scale );
            }
        }
        TIFFWriteEncodedStrip(tif, 0, &data[0], this->height_*this->width_*2);
    }
    else{
        cerr << "ERROR : " << address << " not handled!" <<endl;
        cerr << "bits_per_sample : " << this->bits_per_sample_ << " ";
        cerr << "samples_per_pixel : " << this->samples_per_pixel_ <<endl;
    }

    TIFFClose(tif);
    return;
}

vector<float>& tomo_tiff::operator [](int index_y){
    return this->gray_scale_[index_y];
}

void tomo_super_tiff::make_gaussian_window_(const int size, const float standard_deviation){

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

    //normalize the maximum to 1 for output
    float normalize_ratio = 1.0 / maximum;
    vector< vector< vector<float> > >output_test( gaussian_window_ );
    for(int i=0;i<size;++i){
        for(int j=0;j<size;++j){
            for(int k=0;k<size;++k){
                output_test[i][j][k] *= normalize_ratio;
            }
        }
    }
    mkdir("gaussian",0755);
    string prefix_test("gaussian/");
    for(int i=0;i<size;++i){
        tomo_tiff test_tiff( output_test[i] );
        char number_string[50];
        sprintf(number_string,"%d",i);
        string address_test = prefix_test+string(number_string)+string(".tiff");
        test_tiff.save( address_test.c_str() );
    }

    return;
}

tomo_super_tiff::tomo_super_tiff(const char *address_filelist){

    fstream in_filelist(address_filelist,fstream::in);

    int size_tiffs = -1;
    char prefix[100]={0};
    char original_dir[100]={0};
    getcwd(original_dir,100);

    in_filelist >> size_tiffs;
    in_filelist >> prefix;

    //read filelist first for parallel
    this->tiffs_.resize(size_tiffs);
    this->address_tiffs_.resize(size_tiffs);
    for(int i=0;i<size_tiffs;++i){
        in_filelist >> this->address_tiffs_[i];
    }
    this->prefix_ = string(prefix);

    if(size_tiffs < TIFF_IMAGE_LARGE_SIZE){
        cout << "change working directory to " << prefix <<endl;
        chdir(prefix);

        progressbar *progress = progressbar_new("Reading .tifs",size_tiffs);
        #pragma omp parallel for
        for(int i=0;i<size_tiffs;++i){
            tiffs_[i] = tomo_tiff(this->address_tiffs_[i].c_str());
            #pragma omp critical
            {
                progressbar_inc(progress);
            }
        }
        progressbar_finish(progress);

        cout << "change working directory back to " << original_dir <<endl;
        chdir(original_dir);
        this->tiffs_[0].save("favicon.tif");

    }else{
        cout << "size_tiffs = " << size_tiffs <<endl;
        cout << "reading address only due to the lack of memory." <<endl;
    }
    return;
}

float tomo_super_tiff::Ix_(int x, int y, int z){
    if( x+1 >= this->tiffs_[z][y].size() )
        return this->tiffs_[z][y][x] - this->tiffs_[z][y][x-1];
    else if( x-1 < 0)
        return this->tiffs_[z][y][x+1] - this->tiffs_[z][y][x];
    else
        return ( this->tiffs_[z][y][x+1] - this->tiffs_[z][y][x-1] ) / 2.0;
}

float tomo_super_tiff::Iy_(int x, int y, int z){
    if( y+1 >= this->tiffs_[z].size() )
        return this->tiffs_[z][y][x] - this->tiffs_[z][y-1][x];
    else if( y-1 < 0)
        return this->tiffs_[z][y+1][x] - this->tiffs_[z][y][x];
    else
        return ( this->tiffs_[z][y+1][x] - this->tiffs_[z][y-1][x] ) / 2.0;
}

float tomo_super_tiff::Iz_(int x, int y, int z){
    if( z+1 >= this->tiffs_.size() )
        return this->tiffs_[z][y][x] - this->tiffs_[z-1][y][x];
    else if( z-1 < 0)
        return this->tiffs_[z+1][y][x] - this->tiffs_[z][y][x];
    else
        return ( this->tiffs_[z+1][y][x] - this->tiffs_[z-1][y][x] ) / 2.0;
}

float tomo_super_tiff::summation_within_window_gaussianed_(int x, int y, int z, int size){

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

                sz = sz >= this->tiffs_.size() ? this->tiffs_.size()-1 : sz;
                sy = sy >= this->tiffs_[sz].size() ? this->tiffs_[sz].size()-1 : sy;
                sx = sx >= this->tiffs_[sz][sy].size() ? this->tiffs_[sz][sy].size()-1 : sx;

                summation += this->tiffs_[sz][sy][sx] * this->gaussian_window_[k][j][i];
            }
        }
    }

    return summation;
}

void tomo_super_tiff::down_size(int magnification, const char *save_prefix, float sample_sd){

    vector< tomo_tiff > result;
    int process = 0;

    //init
    cout << "allocting result of down_size..." <<endl;
    result.resize( this->tiffs_.size()/magnification );
    #pragma omp parallel for
    for(int i=0;i<result.size();++i){
        result[i].resize( this->tiffs_[i*magnification].size()/magnification );
        for(int j=0;j<result[i].size();++j){
            result[i][j].resize( this->tiffs_[i*magnification][j*magnification].size()/magnification, 0.0 );
        }
        cout << process << " / " << result.size() <<endl;
        #pragma omp critical
        {
            process++;
        }
    }
    for(int i=0;i<result.size();++i){
        result[i].bits_per_sample_ = 16;
        result[i].samples_per_pixel_ = 1;
    }
    cout << "\t\tdone!" <<endl;

    //gaussian & sampling
    cout << "gaussian and sampling..." <<endl;

    this->make_gaussian_window_(magnification, sample_sd*(float)magnification/2.0);

    process = 0;
    #pragma omp parallel for
    for(int z=0;z<result.size();++z){
        for(int y=0;y<result[z].size();++y){
            for(int x=0;x<result[z][y].size();++x){

                int sx = x*magnification;
                int sy = y*magnification;
                int sz = z*magnification;

                result[z][y][x] = this->summation_within_window_gaussianed_( sx-magnification/2,
                                                                             sy-magnification/2,
                                                                             sz-magnification/2,
                                                                             magnification);
            }
        }
        cout << process << " / " << result.size() <<endl;
        #pragma omp critical
        {
            process++;
        }
    }
    cout << "\t\tdone!" <<endl;

    //save it to save_prefix
    cout << "saving files..." <<endl;
    mkdir(save_prefix, 0755);

    process = 0;
    #pragma omp parallel for
    for(int i=0;i<result.size();++i){
        //make address
        char number_string[50]={0};
        string address(save_prefix);
        sprintf(number_string,"%d",i);
        address += string(number_string) + string(".tiff");
        //save
        result[i].save( address.c_str() );
        cout << address << " saved\t\t" << process << " / " << result.size() <<endl;
        #pragma omp critical
        {
            process++;
        }
    }
    cout << "\t\tdone!" <<endl;

    return;
}

void tomo_super_tiff::make_differential_matrix_(){

    //init
    progressbar *progress = progressbar_new("Initializing",this->tiffs_.size());
    differential_matrix_.resize(this->tiffs_.size());
    #pragma omp parallel for
    for(int i=0;i<differential_matrix_.size();++i){
        this->differential_matrix_[i].resize(this->tiffs_[i].size());
        for(int j=0;j<differential_matrix_[i].size();++j){
            differential_matrix_[i][j].resize(this->tiffs_[i][j].size(),matrix(3,0));
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

    progress = progressbar_new("Calculating",this->tiffs_.size());
    #pragma omp parallel for
    for(int z=0;z<this->tiffs_.size();++z){
        for(int y=0;y<this->tiffs_[z].size();++y){
            for(int x=0;x<this->tiffs_[z][y].size();++x){

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

void tomo_super_tiff::make_tensor_(const int window_size){

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

void tomo_super_tiff::make_nobles_measure_(float measure_constant){

    cout << "making nobles measures..."<<endl;

    //init

    this->measure_.resize( this->tensor_.size() );
    progressbar *progress = progressbar_new("Initializing", this->tensor_.size());
    #pragma omp parallel for
    for(int i=0;i<measure_.size();++i){
        this->measure_[i].resize(this->tensor_[i].size());
        for(int j=0;j<measure_[i].size();++j){
            this->measure_[i][j].resize(this->tensor_[i][j].size(),0.0);
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //calculate measure
    // Noble's cornor measure :
    //      Mc = 2* det(tensor) / ( trace(tensor) + c )

    progress = progressbar_new("Calculating",this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<measure_.size();++i){
        for(int j=0;j<measure_[i].size();++j){
            for(int k=0;k<measure_[i][j].size();++k){
                float trace = this->tensor_[i][j][k].trace();
                this->measure_[i][j][k] = 2 * this->tensor_[i][j][k].det();
                this->measure_[i][j][k] /= trace*trace + measure_constant;
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

void tomo_super_tiff::make_eigen_values_(){

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

void tomo_super_tiff::experimental_measurement(float threshold){

    cout << "making measurement..." <<endl;

    //resize & init
    this->measure_.resize(this->eigen_values_.size());
    progressbar *progress = progressbar_new("Initializing", this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        this->measure_[i].resize(this->eigen_values_[i].size());
        for(int j=0;j<this->measure_[i].size();++j){
            this->measure_[i][j].resize(this->eigen_values_[i][j].size(),0.0);
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //measurement
    progress = progressbar_new("Calculating",this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        for(int j=0;j<this->measure_[i].size();++j){
            for(int k=0;k<this->measure_[i][j].size();++k){
                vector<float> &ev = this->eigen_values_[i][j][k];
                measure_[i][j][k] = 0.3 * ( ev[0] + ev[1] + ev[2]) * ( ev[0] + ev[1] + ev[2]) - ev[0] * ev[1] * ev[2];
                if( threshold > 0 ){
                    if(measure_[i][j][k] >= threshold )
                        measure_[i][j][k] = 1.0;
                    else
                        measure_[i][j][k] = 0.0;
                }
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //normalize
    this->experimental_measurement_normalize_();

    return ;

}

void tomo_super_tiff::make_differential_matrix_(int start_z, int number_z){

    //init
    this->differential_matrix_.resize(this->tiffs_.size());

    /* differential_matrix      j->
     * __                           __
     * |    IxIx    IxIy    IxIz     |  i
     * |                             |  |
     * |    IxIy    IyIy    IyIz     |  v
     * |                             |
     * |    IxIz    IyIz    IzIz     |
     * L_                           _|
     */

    for(int z=0;z<this->tiffs_.size();++z){

        if( z < start_z || z >= start_z+number_z ){ // clear it because it's not needed
            this->differential_matrix_[z].clear();

        }
        else if(this->differential_matrix_[z].size() == 0){ // only calculate one which not calculated before

            //init
            this->differential_matrix_[z].resize(this->tiffs_[z].size());
            #pragma omp parallel for
            for(int y=0;y<this->differential_matrix_[z].size();++y){
                this->differential_matrix_[z][y].resize(this->tiffs_[z][y].size(),matrix(3,0.0));
            }

            //calculating
            #pragma omp parallel for
            for(int y=0;y<this->tiffs_[z].size();++y){
                for(int x=0;x<this->tiffs_[z][y].size();++x){

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

void tomo_super_tiff::make_tensor_(const int window_size, int index_z){
    //init
    this->tensor_.clear();
    this->tensor_.resize(this->tiffs_.size());

    // struct tensor A = sum_u_v_w( gaussian(u,v,w) * differential(u,v,w) )

    //init
    this->tensor_[index_z].resize(this->tiffs_[index_z].size());
    #pragma omp parallel for
    for(int y=0;y<this->tensor_[index_z].size();++y){
        this->tensor_[index_z][y].resize(this->tiffs_[index_z][y].size());
    }

    //calculating
    #pragma omp parallel for
    for(int y=0;y<this->tiffs_[index_z].size();++y){
        for(int x=0;x<this->tiffs_[index_z][y].size();++x){

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

void tomo_super_tiff::eigen_values_initialize_(){

    //init
    this->eigen_values_.resize(this->tiffs_.size());

    progressbar *progress = progressbar_new("EigenValueInit",this->eigen_values_.size());
    #pragma omp parallel for
    for(int i=0;i<this->eigen_values_.size();++i){
        this->eigen_values_[i].resize( this->tiffs_[i].size() );

        for(int j=0;j<this->eigen_values_[i].size();++j){
            this->eigen_values_[i][j].resize( this->tiffs_[i][j].size() );

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

void tomo_super_tiff::make_eigen_values_(int index_z){

    if(this->tiffs_.size() >= TIFF_IMAGE_LARGE_SIZE){ // for the super large data
        //allocate the needed and free others
        this->eigen_values_.clear();
        this->eigen_values_.resize( this->tiffs_.size() );
        this->eigen_values_[index_z].resize( this->tiffs_[index_z].size() );

        #pragma omp for
        for(int j=0;j<this->eigen_values_[index_z].size();++j){
            this->eigen_values_[index_z][j].resize( this->tiffs_[index_z][j].size() );

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

void tomo_super_tiff::experimental_measurement_initialize_(){

    //resize & init
    this->measure_.resize(this->eigen_values_.size());
    progressbar *progress = progressbar_new("MeasurementInit", this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        this->measure_[i].resize(this->eigen_values_[i].size());
        for(int j=0;j<this->measure_[i].size();++j){
            this->measure_[i][j].resize(this->eigen_values_[i][j].size(),0.0);
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    return;
}

void tomo_super_tiff::experimental_measurement_normalize_(){

    //normalize
    float maximum = 0.0;
    for(int i=0;i<this->measure_.size();++i){
        for(int j=0;j<this->measure_[i].size();++j){
            for(int k=0;k<this->measure_[i][j].size();++k){
                maximum = measure_[i][j][k] > maximum ? measure_[i][j][k] : maximum;
            }
        }
    }
    cout << "normalized by " << maximum <<endl;
    this->normalized_measure_ = maximum;

    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        for(int j=0;j<this->measure_[i].size();++j){
            for(int k=0;k<this->measure_[i][j].size();++k){
                measure_[i][j][k] /= maximum;
            }
        }
    }

    return;
}

void tomo_super_tiff::experimental_measurement_(int index_z, float threshold){

    if( this->tiffs_.size() >= TIFF_IMAGE_LARGE_SIZE ){ // for the super large data
        //allocate the needed and free others
        this->measure_.clear();
        this->measure_.resize( this->tiffs_.size() );
        this->measure_[index_z].resize( this->tiffs_[index_z].size() );

        #pragma omp for
        for(int j=0;j<this->measure_[index_z].size();++j){
            this->measure_[index_z][j].resize(this->eigen_values_[index_z][j].size(),0.0);
        }
    }

    #pragma omp parallel for
    for(int j=0;j<this->measure_[index_z].size();++j){
        for(int k=0;k<this->measure_[index_z][j].size();++k){
            vector<float> &ev = this->eigen_values_[index_z][j][k];
            measure_[index_z][j][k] = 0.3 * ( ev[0] + ev[1] + ev[2]) * ( ev[0] + ev[1] + ev[2]) - ev[0] * ev[1] * ev[2];
            if( threshold > 0 ){
                if( measure_[index_z][j][k] >= threshold )
                    measure_[index_z][j][k] = 1.0;
                else
                    measure_[index_z][j][k] = 0.0;
            }
        }
    }


    return;
}

void tomo_super_tiff::neuron_detection(const int window_size, float threshold, const float standard_deviation){

    cout << "making gaussian window with window_size : " << window_size;
    (cout << "\tstandard_deviation : " << standard_deviation ).flush();

    this->make_gaussian_window_(window_size,standard_deviation*(float)window_size/2.0);
    cout << "\tdone!"<<endl;

    if(this->tiffs_.size() < TIFF_IMAGE_MEDIUM_SIZE){ // prevent starvation
        cout << "making differential matrix..." <<endl;
        this->make_differential_matrix_();

        cout << "making struct tensor..." <<endl;
        this->make_tensor_(window_size);

        this->make_eigen_values_();

        this->experimental_measurement( threshold );

    }else if(this->tiffs_.size() < TIFF_IMAGE_LARGE_SIZE){//too large to process normally, using half serial processing

        //init
        this->eigen_values_initialize_();
        this->experimental_measurement_initialize_();

        //load data when needed, free it otherwise
        progressbar *progress = progressbar_new("Calculating",this->tiffs_.size());
        for(int i=0;i<this->tiffs_.size();++i){

            int number_z = window_size;
            int start_z = (i - window_size/2) >= 0 ? (i - window_size/2) : 0 ;
            start_z = (start_z+number_z) <= this->tiffs_.size() ? start_z  : this->tiffs_.size() - number_z;

            this->make_differential_matrix_(start_z, number_z);
            this->make_tensor_(window_size, i);
            this->make_eigen_values_(i);
            this->experimental_measurement_(i, threshold);

            progressbar_inc(progress);
        }
        progressbar_finish(progress);

        //normalize
        this->experimental_measurement_normalize_();

    }else{ //super large, using full serial processing

        //vector<float> maximums_eigen_values(this->tiffs_.size(),0.0);
        vector<float> maximums_measurements(this->tiffs_.size(),0.0);

        //resize eigen value and measurement
        this->measure_.resize(this->tiffs_.size());

        //load data when needed, free it otherwise
        progressbar *progress = progressbar_new("Calculating",this->tiffs_.size());
        for(int i=0;i<this->tiffs_.size();++i){

            int number_z = window_size;
            int start_z = (i - window_size/2) >= 0 ? (i - window_size/2) : 0 ;
            start_z = (start_z+number_z) <= this->tiffs_.size() ? start_z  : this->tiffs_.size() - number_z;

            //load original data needed and free it otherwise
            FILE* err_redir = freopen("tiff_reading_err.txt", "w", stderr);// redirect stderr to err_file

            char original_directory[100] = {0};
            getcwd(original_directory,100);
            chdir(this->prefix_.c_str()); // change to the directory of original data
            #pragma omp for
            for(int j=0;j<this->tiffs_.size();++j){
                if( j < start_z-2 || j >= start_z+number_z+2 ){ // free it
                    this->tiffs_[j].clear();

                }else if(this->tiffs_[j].size() == 0){ // load it
                    this->tiffs_[j] = tomo_tiff( this->address_tiffs_[j].c_str() );
                }
            }
            chdir(original_directory);//change it back
            fclose(err_redir);
            freopen("/dev/tty", "a", stderr); // redirect stderr back to screen

            this->make_differential_matrix_(start_z, number_z);
            this->make_tensor_(window_size, i);
            this->make_eigen_values_(i);
            // todo : save eigen_values[i] for tmp. and find maximum
            this->experimental_measurement_(i, threshold);
            // save measurements[i] for tmp. and find maximum for the first normalization
            for(int j=0;j<this->measure_[i].size();++j){
                for(int k=0;k<this->measure_[i][j].size();++k){
                    maximums_measurements[i] = maximums_measurements[i] > this->measure_[i][j][k] ?
                                maximums_measurements[i] : this->measure_[i][j][k];
                }
            }
            #pragma omp for
            for(int j=0;j<this->measure_[i].size();++j){
                for(int k=0;k<this->measure_[i][j].size();++k){
                    this->measure_[i][j][k] /= maximums_measurements[i];
                }
            }
            char address_tiff[100] = {0};
            mkdir("measurement",0755);
            sprintf(address_tiff, "measurement/%d.tif", i);
            tomo_tiff tiff_mearsure(this->measure_[i]);
            tiff_mearsure.save( address_tiff );

            progressbar_inc(progress);
        }
        progressbar_finish(progress);

        // renormalize the tmp. eigen_values and tmp. measurements
        // find maximum of maximums
        float final_maximum_measurements = 0.0;
        for(int i=0;i<maximums_measurements.size();++i){
            final_maximum_measurements = final_maximum_measurements > maximums_measurements[i] ?
                        final_maximum_measurements : maximums_measurements[i];
        }
        //save info.txt
        fstream out_info("info.txt", fstream::out);
        if(out_info.is_open() == false){
            cerr << "ERROR : cannot open info.txt" <<endl;
            exit(-1);
        }
        out_info << "xyz-size " << this->measure_.back().back().size() << " " << this->measure_.back().size() << " " << this->measure_.size() <<endl;
        out_info << "normalized " << fixed << setprecision(8) << final_maximum_measurements <<endl;
        out_info << "order xyz"<<endl;
        out_info.close();

        #pragma omp for
        for(int i=0;i<this->measure_.size();++i){
            char address_tiff[100] = {0};
            sprintf(address_tiff, "measurement/%d.tif", i);
            tomo_tiff tiff_measure(address_tiff);
            for(int j=0;j<tiff_measure.size();++j){
                for(int k=0;k<tiff_measure[j].size();++k){
                    tiff_measure[j][k] *= maximums_measurements[i];
                    tiff_measure[j][k] /= final_maximum_measurements;
                }
            }
            tiff_measure.save(address_tiff);
        }
    }

    return;
}

void tomo_super_tiff::save_measure(const char *prefix){

    char original_directory[100];
    getcwd(original_directory,100);
    mkdir(prefix, 0755);
    chdir(prefix);
    cout << "changing working directory to " << prefix <<endl;

    //save info.txt
    fstream out_info("info.txt", fstream::out);
    if(out_info.is_open() == false){
        cerr << "ERROR : cannot open info.txt" <<endl;
        exit(-1);
    }
    out_info << "xyz-size " << this->measure_[0][0].size() << " " << this->measure_[0].size() << " " << this->measure_.size() <<endl;
    out_info << "normalized " << fixed << setprecision(8) << this->normalized_measure_ <<endl;
    out_info << "order xyz"<<endl;
    out_info.close();

    //normalize, merge & save
    progressbar *progress = progressbar_new("Saving",this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        //init, normalize & merge
        vector< vector<float> > output_image(this->measure_[i].size());
        for(int j=0;j<output_image.size();++j){
            output_image[j].resize(this->measure_[i][j].size());
            for(int k=0;k<output_image[j].size();++k){
                output_image[j][k] = this->measure_[i][j][k];
            }
        }
        //make address
        char number_string[50]={0};
        sprintf(number_string, "%d", i);
        string address = string(number_string) + string(".tif");
        //save
        tomo_tiff output_tiff(output_image);
        output_tiff.save(address.c_str());

        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    chdir(original_directory);
    cout << "changing working directory back to " << original_directory <<endl;

    return;

}

void tomo_super_tiff::save_measure_merge(const char *prefix){

    char original_directory[100];
    getcwd(original_directory,100);
    mkdir(prefix, 0755);
    chdir(prefix);
    cout << "changing working directory to " << prefix <<endl;

    //normalize, merge & save
    progressbar *progress = progressbar_new("Saving",this->measure_.size());
    #pragma omp parallel for
    for(int i=0;i<this->measure_.size();++i){
        //init, normalize & merge
        vector< vector<float> > output_image(this->measure_[i].size());
        for(int j=0;j<output_image.size();++j){
            output_image[j].resize(this->measure_[i][j].size() + this->tiffs_[i][j].size());
            for(int k=0;k<output_image[j].size();++k){
                output_image[j][k] = k < this->tiffs_[i][j].size() ?
                            this->tiffs_[i][j][k] : this->measure_[i][j][k - this->tiffs_[i][j].size()];
            }
        }
        //making address
        char number_string[50]={0};
        sprintf(number_string, "%d", i);
        string address = string(number_string) + string(".tif");

        //save
        tomo_tiff output_tiff(output_image);
        output_tiff.save(address.c_str());

        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    chdir(original_directory);
    cout << "changing working directory back to " << original_directory <<endl;

    return;
}

void tomo_super_tiff::save_eigen_values_rgb(const char *prefix){

    //find maximum of eigen_values_
    float maximum = -1.0;
    for(int i=0;i<this->eigen_values_.size();++i){
        for(int j=0;j<this->eigen_values_[i].size();++j){
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    maximum = maximum > this->eigen_values_[i][j][k][m] ? maximum : this->eigen_values_[i][j][k][m];
                }
            }
        }
    }

    //save them
    char original_dir[100] = {0};
    mkdir(prefix,0755);
    getcwd(original_dir,100);
    chdir(prefix);

    #pragma omp parallel for
    for(int i=0;i<this->eigen_values_.size();++i){
        //make file name
        char number_string[50] = {0};
        sprintf(number_string,"%d",i);
        string address = string(number_string) + string(".tiff");

        //open file for saving
        TIFF *tif = TIFFOpen(address.c_str(),"w");
        if(tif == NULL){
            cerr << "ERROR : cannot open to save " << prefix << "/" << address <<endl;
        }

        int height = this->eigen_values_[i].size();
        int width = this->eigen_values_[i][0].size();

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

        vector<uint16_t> tmp_data(width * height * 3);
        int index_tmp = 0;
        for(int j=0;j<height;++j){
            for(int k=0;k<width;++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    tmp_data[index_tmp++] = (uint16_t)(this->eigen_values_[i][j][k][m] / maximum * 65535.0);
                }
            }
        }

        TIFFWriteEncodedStrip(tif, 0, &tmp_data[0], width * height * 6);

        TIFFClose(tif);
    }
    chdir(original_dir);

    return;
}

void tomo_super_tiff::save_eigen_values_rgb_merge(const char *prefix){
    //find maximum of eigen_values_
    float maximum = -1.0;
    for(int i=0;i<this->eigen_values_.size();++i){
        for(int j=0;j<this->eigen_values_[i].size();++j){
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    maximum = maximum > this->eigen_values_[i][j][k][m] ? maximum : this->eigen_values_[i][j][k][m];
                }
            }
        }
    }

    //save them
    char original_dir[100] = {0};
    mkdir(prefix,0755);
    getcwd(original_dir,100);
    chdir(prefix);

    #pragma omp parallel for
    for(int i=0;i<this->eigen_values_.size();++i){
        //make file name
        char number_string[50] = {0};
        sprintf(number_string,"%d",i);
        string address = string(number_string) + string(".tiff");

        //open file for saving
        TIFF *tif = TIFFOpen(address.c_str(),"w");
        if(tif == NULL){
            cerr << "ERROR : cannot open to save " << prefix << "/" << address <<endl;
        }

        int width = this->eigen_values_[i].size() + this->tiffs_[i].size();
        int height = this->eigen_values_[i][0].size();

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

        vector<uint16_t> tmp_data(width * height * 3);
        int index_tmp = 0;
        for(int j=0;j<height;++j){
            for(int k=0;k<width;++k){
                for(int m=0;m<3;++m){
                    if(k < this->tiffs_[i].size())
                        tmp_data[index_tmp++] = (uint16_t)(this->tiffs_[i][j][k] * 65535.0);
                    else
                        tmp_data[index_tmp++] = (uint16_t)(this->eigen_values_[i][j][k-tiffs_[i].size()][m] / maximum * 65535.0);
                }
            }
        }

        TIFFWriteEncodedStrip(tif, 0, &tmp_data[0], width * height * 6);

        TIFFClose(tif);
    }
    chdir(original_dir);

    return;
}

void tomo_super_tiff::save_eigen_values_separated(const char *prefix){


    //find maximum of eigen_values_
    float maximum = -1.0;
    for(int i=0;i<this->eigen_values_.size();++i){
        for(int j=0;j<this->eigen_values_[i].size();++j){
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    maximum = maximum > this->eigen_values_[i][j][k][m] ? maximum : this->eigen_values_[i][j][k][m];
                }
            }
        }
    }

    //save them
    char original_dir[100] = {0};
    mkdir(prefix,0755);
    getcwd(original_dir,100);
    chdir(prefix);

    //save info.txt
    fstream out_info("info.txt", fstream::out);
    if(out_info.is_open() == false){
        cerr << "ERROR : cannot open info.txt" <<endl;
        exit(-1);
    }
    out_info << "exyz-size " << this->eigen_values_[0][0][0].size() << " " << this->eigen_values_[0][0].size() << " " << this->eigen_values_[0].size() << " " << this->eigen_values_.size() <<endl;
    out_info << "normalized " << fixed << setprecision(8) << maximum <<endl;
    out_info << "order xyz"<<endl;
    out_info.close();

    //ev0
    for(int t=0;t<3;++t){
        char original_dir_t[100] = {0};
        char number_string[50] = {0};
        sprintf(number_string,"%d",t);
        mkdir(number_string,0755);
        getcwd(original_dir_t,100);
        chdir(number_string);

        #pragma omp parallel for
        for(int i=0;i<this->eigen_values_.size();++i){
            //make file name
            char address[100] = {0};
            sprintf(address,"%d.tiff",i);

            vector< vector<float> > data(this->eigen_values_[i].size());
            for(int j=0;j<data.size();++j){
                data[j].resize(this->eigen_values_[i][j].size(),0.0);
                for(int k=0;k<data[j].size();++k){
                    data[j][k] = this->eigen_values_[i][j][k][t] / maximum;
                }
            }

            tomo_tiff tmp(data);
            tmp.save(address);
        }
        chdir(original_dir_t);
    }
    chdir(original_dir);

    return;
}

void tomo_super_tiff::save_eigen_values_ev(const char *address){

    cout << "saving " << address << "..." <<endl;

    fstream out_ev(address, fstream::out);
    if(out_ev.is_open() == false){
        cerr << "ERROR : cannot open " <<address <<endl;
        exit(-1);
    }

    //find maximum
    progressbar *progress = progressbar_new("Maximum",this->eigen_values_.size());
    float maximum = 0.0;
    for(int i=0;i<this->eigen_values_.size();++i){
        for(int j=0;j<this->eigen_values_[i].size();++j){
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    maximum = maximum > this->eigen_values_[i][j][k][m] ? maximum : this->eigen_values_[i][j][k][m];
                }
            }
        }
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    out_ev << "exyz-size " << this->eigen_values_[0][0][0].size() << " " << this->eigen_values_[0][0].size() << " " << this->eigen_values_[0].size() << " " << this->eigen_values_.size() <<endl;
    out_ev << "normalized " << fixed << setprecision(8) <<  maximum <<endl;
    out_ev << "order xyz"<<endl;

    progress = progressbar_new("Saving",this->eigen_values_.size());
    for(int i=0;i<this->eigen_values_.size();++i){
        for(int j=0;j<this->eigen_values_[i].size();++j){
            for(int k=0;k<this->eigen_values_[i][j].size();++k){
                for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                    out_ev << fixed << setprecision(8) << (float)(this->eigen_values_[i][j][k][m]/maximum) << " ";
                }
            }
        }
        out_ev << endl;
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    out_ev.close();
    return;
}

void tomo_super_tiff::load_eigen_values_ev(const char *address){

    cout << "reading " << address << "..." <<endl;

    fstream in_ev(address, fstream::in);
    if(in_ev.is_open() == false){
        cerr << "ERROR : cannot open " <<address <<endl;
        exit(-1);
    }

    string buffer;
    string order;
    float normalized = 0.0;
    int size_e = 0;
    int size_x = 0;
    int size_y = 0;
    int size_z = 0;

    //exyz-size
    in_ev >> buffer >> size_e >> size_x >> size_y >> size_z;
    //normalized
    in_ev >> buffer >> normalized;
    //order
    in_ev >> buffer >> order;

    //init
    this->eigen_values_.resize(size_z);
    progressbar *progress = progressbar_new("Initialize",this->eigen_values_.size());
    #pragma omp parallel for
    for(int i=0;i<size_z;++i){
        this->eigen_values_[i].resize(size_y);
        for(int j=0;j<size_y;++j){
            this->eigen_values_[i][j].resize(size_x);
            for(int k=0;k<size_x;++k){
                this->eigen_values_[i][j][k].resize(size_e,0.0);
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //read data
    progress = progressbar_new("Reading",this->eigen_values_.size());
    if(order == "xyz"){

        for(int i=0;i<this->eigen_values_.size();++i){
            for(int j=0;j<this->eigen_values_[i].size();++j){
                for(int k=0;k<this->eigen_values_[i][j].size();++k){
                    for(int m=0;m<this->eigen_values_[i][j][k].size();++m){
                        in_ev >> this->eigen_values_[i][j][k][m] ;
                        this->eigen_values_[i][j][k][m] *= normalized;
                    }
                }
            }
            progressbar_inc(progress);
        }
        progressbar_finish(progress);
    }
    else{
        cout << "ERROR : order " << order << " not handled" <<endl;
        exit(-1);
    }

    in_ev.close();

    return;
}

void tomo_super_tiff::load_eigen_values_separated(const char *prefix){
    cout << "reading " << prefix << "..." <<endl;
    cout << "changing directory to " << prefix <<endl;

    char original_directory[100] = {0};
    getcwd(original_directory,99);
    if( chdir(prefix) < 0 ){
        cerr << "ERROR: cannot change directory to " << prefix <<endl;
        exit(-1);
    }

    fstream in_info("info.txt", fstream::in);
    if(in_info.is_open() == false){
        cerr << "ERROR : cannot open info.txt" <<endl;
        exit(-1);
    }

    string buffer;
    string order;
    float normalized = 0.0;
    int size_e = 0;
    int size_x = 0;
    int size_y = 0;
    int size_z = 0;

    //exyz-size
    in_info >> buffer >> size_e >> size_x >> size_y >> size_z;
    //normalized
    in_info >> buffer >> normalized;
    //order
    in_info >> buffer >> order;

    in_info.close();

    //init
    this->eigen_values_.resize(size_z);
    progressbar *progress = progressbar_new("Initialize",this->eigen_values_.size());
    #pragma omp parallel for
    for(int i=0;i<size_z;++i){
        this->eigen_values_[i].resize(size_y);
        for(int j=0;j<size_y;++j){
            this->eigen_values_[i][j].resize(size_x);
            for(int k=0;k<size_x;++k){
                this->eigen_values_[i][j][k].resize(size_e,0.0);
            }
        }
        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    //read data
    progress = progressbar_new("Reading",this->eigen_values_.size());
    if(order == "xyz"){

        #pragma omp parallel for
        for(int i=0;i<this->eigen_values_.size();++i){
            for(int m=0;m<size_e;++m){
                char address_tif[100] = {0};
                sprintf(address_tif, "%d/%d.tiff", m, i);
                tomo_tiff tmp_tiff(address_tif);

                for(int j=0;j<this->eigen_values_[i].size();++j){
                    for(int k=0;k<this->eigen_values_[i][j].size();++k){
                        this->eigen_values_[i][j][k][m] = tmp_tiff[j][k] * normalized;
                    }
                }
            }

            #pragma omp critical
            progressbar_inc(progress);
        }
        progressbar_finish(progress);
    }
    else{
        cout << "ERROR : order " << order << " not handled" <<endl;
        exit(-1);
    }

    cout << "changing directory back to " << original_directory << endl;
    chdir(original_directory);

    return;
}

void create_experimental_data(const char *address){

    int size = 200;

    vector< vector< vector<float> > > volumes;

    //init
    volumes.resize(size);
    #pragma omp parallel for
    for(int i=0;i<size;++i){
        volumes[i].resize(size);
        for(int j=0;j<size;++j){
            volumes[i][j].resize(size);
        }
    }

    float rotation_r = 80.0;

    vector<float> A(3,0.0);
    vector<float> C(3,0.0);
    A[0] = 9.0; // x

    for(int t=0;t<15;++t){

        float theta = M_PI / 2.0 / 15.0 * (float)t;

        float r = 2.0;
        A[0] += (float)r;

        vector<float> B(3,0.0);
        B[0] = A[0]; // x

        C[0] = A[0]; // x
        C[1] = 100.0;
        C[2] = 100.0;

        A[1] = rotation_r * cos(theta) + C[1]; // y
        A[2] = rotation_r * sin(theta) + C[2]; // z

        B[1] = C[1]*2.0 - A[1]; // y
        B[2] = C[2]*2.0 - A[2]; // z


        vector<float> AB = B - A;

        #pragma omp parallel for
        for(int i=0;i<size;++i){
            for(int j=0;j<size;++j){
                for(int k=0;k<size;++k){

                    //create data: cylinder
                    vector<float> P(3,0.0);
                    P[0] = (float)k; // x
                    P[1] = (float)j; // y
                    P[2] = (float)i; // z

                    vector<float> AP = P - A;
                    vector<float> BP = P - B;
                    vector<float> BA = A - B;
                    float AC_length = vector_dot(AP,AB) / vector_length(AB);

                    //it's in the cylinder
                    if( vector_dot(AP,AB) > 0.0 && vector_dot(BP,BA) > 0.0 &&
                            (vector_dot(AP,AP) - AC_length * AC_length) <= ((float)r * (float)r) ){
                        volumes[i][j][k] = 1.0;
                    }
                }
            }
        }

        A[0] += (float)r + 9.0;
    }// for t

    //save volumes
    char original_directory[100];
    getcwd(original_directory,100);
    mkdir(address,0755);
    chdir(address);

    #pragma omp parallel for
    for(int i=0;i<volumes.size();++i){
        char filename[20];
        sprintf(filename, "%d.tif", i);

        tomo_tiff tmp(volumes[i]);
        tmp.save(filename);
    }

    //make filelist
    char filelist_directory[100];
    getcwd(filelist_directory,100);

    fstream out_filelist("exp.txt",fstream::out);

    out_filelist << volumes.size() <<endl;
    out_filelist << filelist_directory <<endl;

    for(int i=0;i<volumes.size();++i){
        char filename[20];
        sprintf(filename, "%d.tif", i);
        out_filelist << filename <<endl;
    }
    out_filelist.close();

    chdir(original_directory);
    return;
}

void merge_measurements(const char *address_filelist, const char *prefix_output){
    cout << "Merging measurements..." <<endl;

    vector< vector< vector<float> > > merge_measure;
    vector< vector< vector<float> > > tmp_measure;
    int size_filelist = 0;

    fstream in_filelist(address_filelist, fstream::in);
    if(!in_filelist.is_open()){
        cerr << "ERROR : cannot open " << address_filelist <<endl;
        exit(-1);
    }

    //read filelist & merge every measurement
    in_filelist >> size_filelist;
    for(int t=0;t<size_filelist;++t){

        char original_directory[100] = {0};
        string buffer_string;
        string address_measurement;
        int size_x = 0;
        int size_y = 0;
        int size_z = 0;
        float normalized = 0.0;
        string order;

        float enlarge_ratio_x = 0.0;
        float enlarge_ratio_y = 0.0;
        float enlarge_ratio_z = 0.0;

        in_filelist >> address_measurement;

        //change working directory
        cout << "changing directory to " << address_measurement << " ..." <<endl;
        getcwd(original_directory, 100);
        chdir(address_measurement.c_str());

        //read info.txt
        fstream in_info("info.txt", fstream::in);
        if(!in_info.is_open()){
            cerr << "ERROR : cannot open info.txt" <<endl;
            exit(-1);
        }

        //xyz-size
        in_info >> buffer_string >> size_x >> size_y >> size_z;
        //normalized
        in_info >> buffer_string >> normalized;
        //order
        in_info >> buffer_string >> order; // ignore for now
        in_info.close();

        //init merge_measure with the size of the first measurement
        if(t == 0){
            progressbar *progress = progressbar_new("Initialize", size_z);
            merge_measure.resize(size_z);
            tmp_measure.resize(size_z);
            #pragma omp parallel for
            for(int i=0;i<size_z;++i){
                merge_measure[i].resize(size_y);
                tmp_measure[i].resize(size_y);
                for(int j=0;j<size_y;++j){
                    merge_measure[i][j].resize(size_x, 0.0);
                    tmp_measure[i][j].resize(size_x, 0.0);
                }
                #pragma omp critical
                progressbar_inc(progress);
            }
            progressbar_finish(progress);

        }else{//init tmp_measure to zero
            progressbar *progress = progressbar_new("Zeroing", tmp_measure.size() );
            #pragma omp parallel for
            for(int i=0;i<tmp_measure.size();++i){
                for(int j=0;j<tmp_measure[i].size();++j){
                    for(int k=0;k<tmp_measure[i][j].size();++k){
                        tmp_measure[i][j][k] = 0.0;
                    }
                }
                #pragma omp critical
                progressbar_inc(progress);
            }
            progressbar_finish(progress);
        }

        //calculate enlarge ratio
        enlarge_ratio_z = (float)merge_measure.size() / (float)size_z;
        enlarge_ratio_y = (float)merge_measure[0].size() / (float)size_y;
        enlarge_ratio_x = (float)merge_measure[0][0].size() / (float)size_x;

        //loading .tifs to tmp_measure
        char progress_label[50];
        sprintf(progress_label, "Loading %d", t);
        progressbar *progress = progressbar_new(progress_label, size_z);

        #pragma omp parallel for
        for(int i=0;i<size_z;++i){

            //make address
            char address_tif[100] = {0};
            sprintf(address_tif, "%d.tif", i);

            tomo_tiff tif(address_tif);
            for(int j=0;j<size_y;++j){
                for(int k=0;k<size_x;++k){

                    int index_z = (int)( enlarge_ratio_z * (float)i );
                    int index_y = (int)( enlarge_ratio_y * (float)j );
                    int index_x = (int)( enlarge_ratio_x * (float)k );

                    for(int sz=0;sz<(int)enlarge_ratio_z;++sz){
                        for(int sy=0;sy<(int)enlarge_ratio_y;++sy){
                            for(int sx=0;sx<(int)enlarge_ratio_x;++sx){

                                //boundary check
                                if( index_z-(int)(enlarge_ratio_z/2.0)+sz < 0 || index_z-(int)(enlarge_ratio_z/2.0)+sz >= tmp_measure.size() ||
                                        index_y-(int)(enlarge_ratio_y/2.0)+sy < 0 || index_y-(int)(enlarge_ratio_y/2.0)+sy >= tmp_measure[0].size() ||
                                        index_x-(int)(enlarge_ratio_x/2.0)+sx < 0 || index_x-(int)(enlarge_ratio_x/2.0)+sx >= tmp_measure[0][0].size())
                                    continue;

                                tmp_measure[ index_z-(int)(enlarge_ratio_z/2.0)+sz ][ index_y-(int)(enlarge_ratio_y/2.0)+sy ][ index_x-(int)(enlarge_ratio_x/2.0)+sx ] += tif[j][k]*normalized;
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            progressbar_inc(progress);
        }
        progressbar_finish(progress);

        //maxing merge_measure with tmp_measure
        sprintf(progress_label, "Maxing %d", t);
        progress = progressbar_new(progress_label, merge_measure.size() );

        #pragma omp parallel for
        for(int i=0;i<merge_measure.size();++i){
            for(int j=0;j<merge_measure[i].size();++j){
                for(int k=0;k<merge_measure[i][j].size();++k){
                    merge_measure[i][j][k] = merge_measure[i][j][k] > tmp_measure[i][j][k] ? merge_measure[i][j][k] : tmp_measure[i][j][k];
                }
            }
            #pragma omp critical
            progressbar_inc(progress);
        }
        progressbar_finish(progress);

        //change directory back to the original one
        chdir(original_directory);
    }
    in_filelist.close();

    //normalize
    float max_merge = -1;

    for(int i=0;i<merge_measure.size();++i){
        for(int j=0;j<merge_measure[i].size();++j){
            for(int k=0;k<merge_measure[i][j].size();++k){
                max_merge = max_merge < merge_measure[i][j][k] ? merge_measure[i][j][k] : max_merge ;
            }
        }
    }

    #pragma omp parallel for
    for(int i=0;i<merge_measure.size();++i){
        for(int j=0;j<merge_measure[i].size();++j){
            for(int k=0;k<merge_measure[i][j].size();++k){
                merge_measure[i][j][k] /= max_merge;
            }
        }
    }

    //output merge_measure to prefix_output
    char original_directory[100] = {0};
    cout << "change directory to " << prefix_output <<endl;
    getcwd(original_directory, 100);
    mkdir(prefix_output, 0755);
    chdir(prefix_output);

    //save info.txt
    fstream out_info("info.txt", fstream::out);
    if(out_info.is_open() == false){
        cerr << "ERROR : cannot open info.txt" <<endl;
        exit(-1);
    }
    out_info << "xyz-size " << merge_measure[0][0].size() << " " << merge_measure[0].size() << " " << merge_measure.size() <<endl;
    out_info << "normalized " << fixed << setprecision(8) << max_merge <<endl;
    out_info << "order xyz"<<endl;
    out_info.close();

    progressbar *progress = progressbar_new("Save", merge_measure.size());
    #pragma omp parallel for
    for(int i=0;i<merge_measure.size();++i){
        char address_tif[100] = {0};
        sprintf(address_tif, "%d.tif", i);
        tomo_tiff tif( merge_measure[i] );
        tif.save( address_tif, 50000 );

        #pragma omp critical
        progressbar_inc(progress);
    }
    progressbar_finish(progress);

    cout << "change directory back to " << original_directory <<endl;
    chdir(original_directory);
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
