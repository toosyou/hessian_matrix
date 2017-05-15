#ifndef TOMO_TIFF
#define TOMO_TIFF

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <iomanip>
#include <sstream>

extern "C"{
#include "progressbar.h"
#include "statusbar.h"
}

#include <omp.h>

#define TIFF_IMAGE_MEDIUM_SIZE 500
#define TIFF_IMAGE_LARGE_SIZE 1000

using namespace std;

class volume;

void merge_measurements(const char* address_filelist, const char* prefix_output);

vector<float> operator -(vector<float> &a, vector<float> &b);
vector<float> operator +(vector<float> &a, vector<float> &b);
float vector_dot(vector<float> &a, vector<float> &b);
float vector_length(vector<float> &a);

void create_experimental_data(const char* address);

class slice{

    string address_;
    unsigned int height_;
    unsigned int width_;
    uint8_t bits_per_sample_;
    int samples_per_pixel_;

    vector< vector<float> > gray_scale_;

    public:

    slice(){
        this->height_ = -1;
        this->width_ = -1;
        this->bits_per_sample_ = 0;
        this->samples_per_pixel_ = -1;
        this->gray_scale_.clear();
    }
    slice(vector< vector<float> >& data){
        this->gray_scale_ = data;
        this->height_ = data.size();
        this->width_ = data[0].size();
        this->bits_per_sample_ = 16;
        this->samples_per_pixel_ = 1;
    }

    vector<float>& operator [](int index_y);
    int size(void){return this->gray_scale_.size();}
    void resize(int size){
        this->gray_scale_.resize(size);
        this->height_ = size;
    }
    void clear(){
        this->gray_scale_.clear();
    }

    friend class volume;

};

class matrix{

    vector< vector<float> > number_;

    public:

    matrix(const int size = 0, const float number = 0.0){
        this->resize(size,number);
    }

    vector<float>& operator [](const int index_i){
        return this->number_[index_i];
    }
    void resize(const int size,const float number = 0.0){
        this->number_.resize(size);
        for(int i=0;i<size;++i){
            this->number_[i].resize(size,number);
        }
    }
    int size(void){
        return this->number_.size();
    }

    const matrix& operator =(const matrix& b){
        this->number_ = b.number_;
        return *this;
    }

    matrix operator *(float ratio){
        matrix rtn;
        rtn.number_ = this->number_;
        for(int i=0;i<rtn.number_.size();++i){
            for(int j=0;j<rtn.number_[i].size();++j){
                rtn.number_[i][j] = this->number_[i][j] * ratio;
            }
        }
        return rtn;
    }

    void operator +=(matrix b){
        if(this->number_.size() != b.size()){
            cerr << "matrixes' size don't match" <<endl;
            return;
        }
        for(int i=0;i<this->number_.size();++i){
            for(int j=0;j<this->number_[i].size();++j){
                this->number_[i][j] += b.number_[i][j];
            }
        }
    }

    float det(){
        if(number_.size() == 3){
            /*      00      01      02
             *
             *      10      11      12
             *
             *      20      21      22
             */
            float ans = 0.0;
            ans += number_[0][0] * number_[1][1] * number_[2][2];
            ans += number_[0][1] * number_[1][2] * number_[2][0];
            ans += number_[1][0] * number_[2][1] * number_[0][2];

            ans -= number_[0][2] * number_[1][1] * number_[2][0];
            ans -= number_[0][1] * number_[1][0] * number_[2][2];
            ans -= number_[0][0] * number_[2][1] * number_[1][2];
            return ans;
        }
        else if(number_.size() == 2){
            float ans = 0.0;
            /*      00      01
             *
             *      10      11
             */
            ans += number_[0][0] * number_[1][1];
            ans -= number_[0][1] * number_[1][0];
            return ans;
        }
        else{
            cerr << "matrix size : " << number_.size() << " determine not handled!" <<endl;
            return 0.0;
        }
    }

    float trace(){
        float ans = 0.0;
        for(int i=0;i<number_.size();++i){
            ans += number_[i][i];
        }
        return ans;
    }
};

class volume{

    string prefix_;
    vector<slice> intensity_;//[z][y][x]
    vector< vector< vector<float> > > gaussian_window_;
    vector< vector< vector<matrix> > >differential_matrix_;
    vector< vector< vector<matrix> > >tensor_;
    vector< vector< vector<float> > >measure_;
    vector< vector< vector< vector<float> > > >eigen_values_;

    float normalized_measure_;

    // Noble's cornor measure :
    //      Mc = 2* det(tensor) / ( trace(tensor) + c )

    void make_gaussian_window_(const int size, const float standard_deviation);
    void make_differential_matrix_();
    void make_tensor_(const int window_size);
    void make_eigen_values_();

    //serial process
    void make_differential_matrix_(int start_z, int number_z);
    void make_tensor_(const int window_size, int index_z);
    void eigen_values_initialize_();
    void make_eigen_values_(int index_z);

    float Ix_(int x, int y, int z);
    float Iy_(int x, int y, int z);
    float Iz_(int x, int y, int z);

    float summation_within_window_gaussianed_(int x, int y, int z, int size);

    public:

    volume(vector<vector<vector<float> > >);
    volume(){}

    void calculate_hessian(const int window_size, float threshold = 0.0000015, const float standard_deviation=0.8);

    int size_original_data(void){return this->intensity_.size();}

    vector< vector< vector< vector<float> > > > get_eigen_value(){return this->eigen_values_;}

};

#endif // TOMO_TIFF
