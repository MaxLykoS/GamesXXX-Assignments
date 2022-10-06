#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

const float H = 0.017453293f;

int main()
{
    Eigen::MatrixXf m(3, 1);
    Eigen::Matrix3f rot, trans;
    float hu = 45 * H;
    m << 2.0, 1.0, 1.0;
    rot << std::cos(hu), -std::sin(hu), 0, 
           std::sin(hu), std::cos(hu),  0, 
           0, 0, 1;
    trans << 1, 0, 1,
             0, 1, 2,
             0, 0, 1;
    m = trans * rot * m;
    std::cout << m << std::endl;

    return 0;
}