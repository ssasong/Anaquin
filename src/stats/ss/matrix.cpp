#include <fstream>
#include <rapidcsv.h>
#include "stats/ss/matrix.hpp"

using namespace SS;

Matrix MatrixTools::readMatrix(const std::string &file, std::size_t nrows, std::size_t ncols, bool skipHead)
{
    Matrix x(nrows, ncols);
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        for (auto j = 0; j < row.size(); j++)
        {
            x(i, j) = stod(row[j]);
        }
    }

    return x;
}

std::vector<double> MatrixTools::toSTDVect(const Vector &v)
{
    std::vector<double> r;
    r.resize(v.size());
    Vector::Map(&r[0], v.size()) = v;
    return r;
}