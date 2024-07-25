#include "..\src\MeshViewer\MeshDefinition.h"
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
namespace CABT
{
    typedef double scalar;
    typedef Eigen::Matrix<scalar, -1, -1> MatrixXs;
    typedef Eigen::Matrix<scalar, 3, 20> mat3_10;
    typedef Eigen::Matrix<scalar, 3, 20> mat3_20;
    typedef Eigen::Matrix<scalar, 3, 27> mat3_27;
    typedef Eigen::Matrix<scalar, 3, 64> mat3_64;
    typedef Eigen::Matrix<scalar, 3, 1> vec3s;
    typedef Eigen::Matrix<scalar, 20, 4> mat20_4;
    typedef Eigen::Matrix<scalar, 20, 20> mat20_20;
    typedef Eigen::Matrix<scalar, 3, 3> mat3_3;
    typedef Eigen::Tensor<scalar, 3> tensor3;
}