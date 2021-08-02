#ifndef SS_MATRIX_HPP
#define SS_MATRIX_HPP

#include <vector>
#include <Eigen/Dense>

namespace SS
{
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;

    struct MatrixTools
    {
        static std::vector<double> toSTDVect(const Vector &v);
        static Matrix readMatrix(const std::string &, std::size_t, std::size_t, bool skipHead);
    };
}

#endif