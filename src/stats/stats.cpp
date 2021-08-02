#include <rapidcsv.h>
#include "tools/tools.hpp"
#include "stats/stats.hpp"
#include "stats/ss/stats.hpp"

using namespace Anaquin;

template <typename T> bool isMissing(const T& x) { return x == MISSING || x == "." || x == "-"; }

SS::Matrix Anaquin::TMM(const SS::Matrix &x)
{
    const auto total = x.colwise().sum();
    const auto totalA = total.mean();

    auto calcFactorWeighted = [&](std::size_t i, std::size_t j = 0, double logratioTrim = 0.3, double sumTrim = 0.05)
    {
        if (i == j) { return 1.0; }
        
        auto obs = x.col(i);
        auto ref = x.col(j);
        
        const auto nO = obs.sum();
        const auto nR = ref.sum();
        
        auto v    = std::vector<double>();
        auto logR = std::vector<double>();
        auto absE = std::vector<double>();

        for (auto k = 0; k < obs.size(); k++)
        {
            logR.push_back(log2((obs[k]/nO) / (ref[k]/nR)));
            absE.push_back((log2(obs[k]/nO) + log2(ref[k]/nR)) / 2.0);
            v.push_back(((nO - obs[k])/nO/obs[k]) + ((nR - ref[k])/nR/ref[k]));
        }
        
        v = removeNA(v);
        logR = removeNA(logR);
        absE = removeNA(absE);
        assert(logR.size() == absE.size());

        const auto n = logR.size();
        const auto loL = floor(n * logratioTrim) + 1.0;
        const auto hiL = n + 1.0 - loL;
        const auto loS = floor(n * sumTrim) + 1.0;
        const auto hiS = n + 1.0 - loS;

        auto rankLR = rank(logR);
        auto rankAE = rank(absE);

        auto logRK = std::vector<double>();
        auto absEK = std::vector<double>();
        
        auto sum1 = 0.0;
        auto sum2 = 0.0;
        
        for (auto i = 0; i < logR.size(); i++)
        {
            const auto keep = (rankLR[i] >= loL & rankLR[i] <= hiL) & (rankAE[i] >= loS & rankAE[i] <= hiS);
            
            if (keep)
            {
                sum1 += logR[i] / v[i];
                sum2 += 1.0 / v[i];
            }
        }
        
        return (std::pow(2.0, sum1 / sum2));
    };
    
    auto fk = std::vector<double>();
    for (auto i = 0; i < x.cols(); i++) {
        fk.push_back(calcFactorWeighted(i));
    }
    
    const auto tmp = exp(SS::mean(logV(fk)));
    for (auto i = 0; i < fk.size(); i++) {
        fk[i] /= tmp;
    }
    
    auto newTotal = std::vector<double>();
    for (auto i = 0; i < total.size(); i++) {
        newTotal.push_back(total[i] / totalA);
    }
    
    for (auto i = 0; i < fk.size(); i++) {
        fk[i] = fk[i] * newTotal[i];
    }

    auto xx = x;
    for (auto j = 0; j < xx.cols(); j++)
    {
        for (auto i = 0; i < xx.rows(); i++)
        {
            xx(i,j) /= fk[j];
        }
    }
    
    return xx;
}

static std::vector<Token> RGetColumns(rapidcsv::Document &doc)
{
    std::vector<Token> toks;
    
    for (auto i = 0; i < doc.GetColumnCount(); i++)
    {
        toks.push_back(doc.GetColumnName(i));
    }
    
    return toks;
}

static void RWriteHeader(std::ofstream &w, const std::vector<Token> &toks)
{
    for (auto i = 0; i < toks.size(); i++)
    {
        w << toks.at(i) << (i != toks.size() - 1 ? "\t" : "");
    }
}

Count Anaquin::RFilterR(const FileName &src, const FileName &dst, const std::set<std::string> &cols)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    RWriteHeader(w, heads);
    
    std::vector<Token> tmp;
    tmp.resize(heads.size());
    
    std::map<std::string, std::size_t> m;
    
    for (auto iter = cols.begin(); iter != cols.end(); iter++)
    {
        // Guarantee there's always an entry
        m[*iter] = std::find(heads.begin(), heads.end(), *iter) - heads.begin();
    }
    
    // Number of filtered
    auto k = 0;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        auto shouldWrite = [&]()
        {
            for (auto iter = m.begin(); iter != m.end(); iter++)
            {
                if (isMissing(row[iter->second]))
                {
                    return false;
                }
            }
            
            return true;
        };
        
        if (shouldWrite())
        {
            k++;
            w << std::endl;
            
            for (auto j = 0; j < row.size(); j++)
            {
                w << row[j] << (j != row.size() - 1 ? "\t" : "");
            }
        }
    }
    
    w.close();
    return k;
}

SequinStats Anaquin::RLinear(const FileName &src, const Label &name, const Label &x, const Label &y)
{
    SequinStats s;
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    const auto i = std::find(heads.begin(), heads.end(), x) - heads.begin();    // Index for independent
    const auto j = std::find(heads.begin(), heads.end(), y) - heads.begin();    // Index for dependent
    const auto n = std::find(heads.begin(), heads.end(), name) - heads.begin(); // Index for ID

    assert(i != heads.size());
    assert(j != heads.size());
    assert(n != heads.size());

    for (auto k = 0; k < doc.GetRowCount(); k++)
    {
        const auto row = doc.GetRow<std::string>(k);
        
        if (!isMissing(row[i]) && !isMissing(row[j]))
        {
            s.add(row[n], stod(row[i]), stod(row[j]));
        }
    }
    
    return s;
}

std::vector<std::string> Anaquin::RReadTSV(const FileName &file, const Label &key)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto cols = RGetColumns(doc);
    
    auto ki = std::find(cols.begin(), cols.end(), key) - cols.begin();
    assert(ki != cols.size());
    
    std::vector<std::string> r;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        r.push_back(doc.GetRow<std::string>(i)[ki]);
    }
    
    return r;
}

std::map<std::string, std::string> Anaquin::RReadTSV(const FileName &file, const Label &key, const Label &val)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto cols = RGetColumns(doc);
    
    auto ki = std::find(cols.begin(), cols.end(), key) - cols.begin();
    auto vi = std::find(cols.begin(), cols.end(), val) - cols.begin();
    assert(ki != cols.size());
    assert(vi != cols.size());
    
    std::map<std::string, std::string> r;

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        r[row[ki]] = row[vi];
    }
    
    return r;
}

std::map<double, double> Anaquin::RBinaryTSV(const FileName &file, const Label &key, const Label &val)
{
    std::map<double, double> r;
    
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto cols = RGetColumns(doc);
    assert(cols.size() == 2);
    
    auto ki = std::find(cols.begin(), cols.end(), key) - cols.begin();
    auto vi = std::find(cols.begin(), cols.end(), val) - cols.begin();

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        assert(row.size() == 2);
        r[stod(row[ki])] = (row[vi] == MISSING ? NAN : stod(row[vi]));
    }
    
    return r;
}

void Anaquin::RMeanCV(const FileName &src, const FileName &dst, const Label &key, const Label &val, Count nExp, Count nObs, Count nCV, float d)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp3 = tmpFile();
    const auto tmp4 = tmpFile();
    
    RFilterC(src,  tmp1, std::set<std::string> { key, val }, true);
    RFilterR(tmp1, tmp2, std::set<std::string> { key });
    RFilterR(tmp2, tmp1, std::set<std::string> { val });
    
    RAggregateMean(tmp1, tmp2, key, d);  // tmp2 holds the results
    RAggregateSD(tmp1, tmp3, key, d);    // tmp3 holds the results
    RAggregateCount(tmp1, tmp4, key, d); // tmp4 holds the results

    auto m = RBinaryTSV(tmp2, key, val); // Mean
    auto s = RBinaryTSV(tmp3, key, val); // SD
    auto n = RBinaryTSV(tmp4, key, val); // Count
    
    assert(m.size() == s.size());
    assert(m.size() == n.size());
    
    std::ofstream w(dst);
    w << "NAME\tCOUNT\tMEAN\tCV";

    for (const auto &i : m)
    {
        assert(s.count(i.first));
        assert(n.count(i.first));
        
        // Coefficient of variation
        const auto cv = !s[i.first] ? "NA" : toString((double) s[i.first] / m[i.first], nCV);

        w << std::endl << toString(i.first, nExp) << "\t" << n.at(i.first) << "\t" << toString(m[i.first], nObs) << "\t" << cv;
    }
    
    w.close();
}

void Anaquin::RNoLast(const FileName &src,
                      const FileName &dst,
                      const Column   &col,
                      const Label    &sep)
{
    std::ofstream w(dst);

    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    RWriteHeader(w, heads);

    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    std::map<Column, double> m;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        auto row = doc.GetRow<std::string>(i);        
        w << std::endl;
        
        for (auto j = 0; j < row.size(); j++)
        {
            auto val = row[j];
            
            if (j == k)
            {
                val = noLast(val, sep);
            }
            
            w << val << (j != row.size() - 1 ? "\t" : "");
        }
    }
}

void Anaquin::RAggregate(const FileName &src, const FileName &dst, const Label &col, Apply f, float d, Imputation impute)
{
    std::ofstream w(dst);
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    
    auto cols = RGetColumns(doc);
    assert(cols.size() >= 2);
    RWriteHeader(w, cols);
    
    // Column index
    auto k = std::find(cols.begin(), cols.end(), col) - cols.begin();
    
    std::set<Label>  k1;
    std::set<double> k2;

    std::map<Label, std::map<Label, std::vector<double>>> m;

    for (auto i = 0; i < cols.size(); i++)
    {
        if (i != k) { m[cols[i]]; }
    }
    
    auto tmp = std::map<double, std::string>();

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        assert(row.size() >= 2);
        
        for (auto j = 0; j < row.size(); j++)
        {
            if (j != k)
            {
                if (row[j] == MISSING)
                {
                    continue;
                }
                
                auto x = row[k];
                
                if (isNumber(x) && d > 0)
                {
                    auto found = false;
                    const auto tmp2 = stod(x);
                    
                    for (auto &i : tmp)
                    {
                        if (fabs(i.first - tmp2) < d)
                        {
                            found = true;
                            x = i.second;
                            break;
                        }
                    }
                    
                    if (!found) { tmp[stod(x)] = x; }
                }
                
                if (x == MISSING)
                {
                    switch (impute)
                    {
                        case Imputation::None:   { break; }
                        case Imputation::ToZero: { x = "0"; break; }
                        case Imputation::Remove: { continue; }
                    }
                }
                
                k1.insert(x); if (isNumber(x)) { k2.insert(stod(x)); }
                m[cols[j]][x].push_back(stod(row[j]));
            }
        }
    }
    
    std::vector<Label> keys(k1.begin(), k1.end());
    
    if (k1.size() == k2.size())
    {
        std::sort(keys.begin(), keys.end(), [] (const std::string &x, const std::string &y) {
            return std::stod(x) < std::stod(y);
        });
    }

    for (const auto &key : keys)
    {
        w << std::endl << key;

        for (const auto &col : cols)
        {
            if (m.count(col))
            {
                if (m[col].count(key))
                {
                    w << "\t" << f(m[col][key]);
                }
                else
                {
                    w << "\t" << MISSING;
                }
            }
        }
    }

    w.close();
}

LadderResult Anaquin::RLadder(const FileName &src)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    
    // Only ladders
    RGrep(src, tmp1, "NAME", "_LD_");
    
    // Writer for adding "Copy" column
    std::ofstream w1(tmp2);
    
    rapidcsv::Document doc(tmp1, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    w1 << "COPY\t"; RWriteHeader(w1, heads);
    
    // Column index for the ID
    auto k = std::find(heads.begin(), heads.end(), "NAME") - heads.begin();
    
    assert(k != heads.size());
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        w1 << std::endl;
        
        auto cn = -1;
        if (isSubstr(row[k], "_A")) { cn = 1; }
        if (isSubstr(row[k], "_B")) { cn = 2; }
        if (isSubstr(row[k], "_C")) { cn = 4; }
        if (isSubstr(row[k], "_D")) { cn = 8; }
        assert(cn != -1);
        
        w1 << cn;
        for (auto j = 0; j < row.size(); j++)
        {
            w1 << "\t" << row[j];
        }
    }
    
    // tmp2 has the results
    w1.close();
    
    const auto key = "COPY";
    const auto val = "Q50";
    
    RFilterC(tmp2, tmp1, std::set<Label> { key, val }, true);
    const auto nonZero = !RFilterR(tmp1, tmp2, std::set<Label> { val }); // tmp2 holds the results
    
    std::map<double, double> mu, sd, q0, q25, q50, q75, q100;
    
    if (!nonZero)
    {
        RAggregateMean(tmp2, tmp1, key); mu   = RBinaryTSV(tmp1, key, val);
        RAggregateSD  (tmp2, tmp1, key); sd   = RBinaryTSV(tmp1, key, val);
        RAggregateMin (tmp2, tmp1, key); q0   = RBinaryTSV(tmp1, key, val);
        RAggregateQ25 (tmp2, tmp1, key); q25  = RBinaryTSV(tmp1, key, val);
        RAggregateQ50 (tmp2, tmp1, key); q50  = RBinaryTSV(tmp1, key, val);
        RAggregateQ75 (tmp2, tmp1, key); q75  = RBinaryTSV(tmp1, key, val);
        RAggregateMax (tmp2, tmp1, key); q100 = RBinaryTSV(tmp1, key, val);
    }
    
    assert(mu.size() == sd.size());
    assert(sd.size() == q0.size());
    assert(q0.size() == q25.size());
    assert(q25.size() == q50.size());
    assert(q50.size() == q75.size());
    assert(q75.size() == q100.size());
    
    std::vector<double> keys;
    for (const auto &i : mu) { keys.push_back(i.first); }
    
    std::vector<double> ros;
    std::vector<double> aps;
    
    LadderResult x;

    for (auto i = 0; i < keys.size(); i++)
    {
        const auto cv = mu[keys[i]] == 0.0 ? std::numeric_limits<double>::quiet_NaN() : (double) sd[keys[i]] / mu[keys[i]];
        const auto ro = !i ? NAN : q50[keys[i]] / q50[keys[i-1]];
        const auto cn = (int) std::pow(2, i);
        
        x.cn[cn].cv    = cv;
        x.cn[cn].sd    = sd[keys[i]];
        x.cn[cn].mean  = mu[keys[i]];
        x.cn[cn].q0    = q0[keys[i]];
        x.cn[cn].q25   = q25[keys[i]];
        x.cn[cn].q50   = q50[keys[i]];
        x.cn[cn].q75   = q75[keys[i]];
        x.cn[cn].q100  = q100[keys[i]];
        x.cn[cn].ratio = ro;
        
        if (cn > 2)          { aps.push_back(mu[keys[i]]); }
        if (!std::isnan(ro)) { ros.push_back(ro); }
    }
    
    x.meanR = ros.empty() ? NAN : SS::mean(ros);
    
    return x;
}

double Anaquin::RLadTable(const FileName &src, const FileName &dst, const Label &c1, const Column &c2)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp3 = tmpFile();
    
    // Only ladders
    RGrep(src, tmp1, c1, "_LD_");

    // Writer for adding "Copy" column
    std::ofstream w1(tmp2);

    rapidcsv::Document doc(tmp1, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    w1 << "COPY\t"; RWriteHeader(w1, heads);
    
    // Column index for the ID
    auto k = std::find(heads.begin(), heads.end(), c1) - heads.begin();

    assert(k != heads.size());
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        w1 << std::endl;

        auto cn = -1;
        if (isSubstr(row[k], "_A")) { cn = 1; }
        if (isSubstr(row[k], "_B")) { cn = 2; }
        if (isSubstr(row[k], "_C")) { cn = 4; }
        if (isSubstr(row[k], "_D")) { cn = 8; }
        assert(cn != -1);
        
        w1 << cn;
        for (auto j = 0; j < row.size(); j++)
        {
            w1 << "\t" << row[j];
        }
    }
    
    // tmp2 has the results
    w1.close();
    
    const auto key = "COPY";
    const auto val = c2;
    
    RFilterC(tmp2, tmp1, std::set<Label> { key, val }, true);
    const auto nonZero = !RFilterR(tmp1, tmp2, std::set<Label> { val }); // tmp2 holds the results

    std::map<double, double> mu, sd, q0, q25, q50, q75, q100;
    
    if (!nonZero)
    {
        RAggregateSD (tmp2, tmp1, key); sd   = RBinaryTSV(tmp1, key, val);
        RAggregateMin(tmp2, tmp1, key); q0   = RBinaryTSV(tmp1, key, val);
        RAggregateQ25(tmp2, tmp1, key); q25  = RBinaryTSV(tmp1, key, val);
        RAggregateQ50(tmp2, tmp1, key); q50  = RBinaryTSV(tmp1, key, val);
        RAggregateQ75(tmp2, tmp1, key); q75  = RBinaryTSV(tmp1, key, val);
        RAggregateMax(tmp2, tmp1, key); q100 = RBinaryTSV(tmp1, key, val);
        RAggregateSum(tmp2, tmp3, "COPY");
        RAggregateMean(tmp3, tmp1, key); mu = RBinaryTSV(tmp1, key, val);
    }
    
    assert(mu.size()  == sd.size());
    assert(sd.size()  == q0.size());
    assert(q0.size()  == q25.size());
    assert(q25.size() == q50.size());
    assert(q50.size() == q75.size());
    assert(q75.size() == q100.size());

    std::ofstream w2(dst);
    w2 << "COPY\tSUM\tSD\tCV\tQ0\tQ25\tQ50\tQ75\tQ100\tRATIO";

    std::vector<double> keys;
    for (const auto &i : mu) { keys.push_back(i.first); }

    // List of ratios
    std::vector<double> ros;
    
    for (auto i = 0; i < keys.size(); i++)
    {
        const auto cv = mu[keys[i]] == 0.0 ? "NA" : toString((double) sd[keys[i]] / mu[keys[i]], 4.0);
        const auto ro = !i ? NAN : mu[keys[i]] / mu[keys[i-1]];
        const auto cn = std::to_string((int) std::pow(2, i));
        
        #define X(x) "\t" << x[keys[i]]
        w2 << std::endl << cn << X(mu) << X(sd) << "\t" << cv << X(q0) << X(q25) << X(q50) << X(q75) << X(q100) << "\t" << replaceNA(ro, 4);
        
        if (!std::isnan(ro))
        {
            ros.push_back(ro);
        }
    }
    
    w2.close();
    
    // Mean ratios
    return ros.empty() ? NAN : SS::mean(ros);
}

void Anaquin::RApply(const FileName &src, const FileName &dst, const Label &col, RApplyF f)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    RWriteHeader(w, heads);
    
    // Column index
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        auto row = doc.GetRow<std::string>(i);
        row[k] = f(row[k]);
        
        w << std::endl;
        
        for (auto j = 0; j < row.size(); j++)
        {
            w << row[j] << (j != row.size() - 1 ? "\t" : "");
        }
    }
    
    w.close();
}

bool Anaquin::RHead(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    return (std::find(heads.begin(), heads.end(), col) - heads.begin()) < heads.size();
}

Count Anaquin::RCountCols(const FileName &file)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    return RGetColumns(doc).size();
}

std::vector<double> Anaquin::RArith(const FileName &file, const Column &x, const Column &y, Arithmetic a, bool skipNA)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto c1 = std::find(heads.begin(), heads.end(), x) - heads.begin();
    auto c2 = std::find(heads.begin(), heads.end(), y) - heads.begin();
    assert(c1 < heads.size() && c2 < heads.size());
    
    std::vector<double> m;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        auto row = doc.GetRow<std::string>(i);
        
        if (row[c1] == "NA" || row[c2] == "NA")
        {
            if (!skipNA)
            {
                m.push_back(NAN);
            }

            continue;
        }
        
        switch (a)
        {
            case Arithmetic::Add:
            {
                m.push_back(stod(row[c2]) + stod(row[c1]));
                break;
            }

            case Arithmetic::Subtract:
            {
                m.push_back(stod(row[c2]) - stod(row[c1]));
                break;
            }

            case Arithmetic::Multipy:
            {
                m.push_back(stod(row[c2]) * stod(row[c1]));
                break;
            }

            case Arithmetic::Divide:
            {
                m.push_back(stod(row[c2]) / stod(row[c1]));
                break;
            }
        }        
    }
    
    return m;
}

Count Anaquin::RCount(const FileName &file, const Column &col, const std::string &x, bool exact)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(col.empty() || k < heads.size());

    Count n = 0;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto j = doc.GetRow<std::string>(i)[k];
        
        if (col.empty() || ((exact && j == x) || (!exact && isSubstr(j, x))))
        {
            n++;
        }
    }
    
    return n;
}

std::vector<std::string> Anaquin::RList(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
 
    std::vector<std::string> x;

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        x.push_back(doc.GetRow<std::string>(i)[k]);
    }

    return x;
}

double Anaquin::RSum(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    double n = 0;

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto x = doc.GetRow<std::string>(i)[k];
        
        if (x != MISSING)
        {
            n += stod(x);
        }
    }
    
    return n;
}

double Anaquin::RMean(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    std::vector<double> r;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto x = doc.GetRow<std::string>(i)[k];
        
        if (x != MISSING)
        {
            r.push_back(stod(x));
        }
    }
    
    return SS::mean(r);
}

void Anaquin::RGrep(const FileName &src, const FileName &dst, const Label &col, const std::set<std::string> &vals, bool keep, bool isNum, bool exact)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    RWriteHeader(w, heads);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        const auto match = std::any_of(vals.begin(), vals.end(), [&](const std::string &x)
        {
            return (!isNum && (exact ? row[k] == x : isSubstr(row[k], x))) || (isNum && stod(row[k]) == stod(x));
        });
        
        if ((keep && match) || (!keep && !match))
        {
            w << std::endl;
            
            for (auto j = 0; j < row.size(); j++)
            {
                w << row[j] << (j != row.size() - 1 ? "\t" : "");
            }
        }
    }
    
    w.close();
}

void Anaquin::RAggregateMean(const FileName &src, const FileName &dst, const Label &col, float d, Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return SS::mean(x);
    }, d, impute);
}

void Anaquin::RAggregateCount(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return std::count_if(x.begin(), x.end(), [&](double x)
        {
            return x > 0 && !std::isinf(x) && !std::isnan(x);
        });
    }, d);
}

void Anaquin::RAggregateSD(const FileName &src, const FileName &dst, const Label &col, float d,  Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return SS::SD(x);
    }, d, impute);
}

void Anaquin::RAggregateMin(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return min(x);
    }, d);
}

void Anaquin::RAggregateQ25(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return quant(x, 0.25);
    }, d);
}

void Anaquin::RAggregateQ50(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return med(x);
    }, d);
}

void Anaquin::RAggregateQ75(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return quant(x, 0.75);
    }, d);
}

void Anaquin::RAggregateMax(const FileName &src, const FileName &dst, const Label &col, float d)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return max(x);
    }, d);
}

void Anaquin::RAggregateSum(const FileName &src, const FileName &dst, const Label &col, float d, Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return std::accumulate(x.begin(), x.end(), 0.0);
    }, d, impute);
}

void Anaquin::RFilterC(const FileName &src, const FileName &dst, const std::set<std::size_t> &m, bool keep)
{
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);

    std::ofstream w(dst);
    std::vector<Token> tmp;
    
    for (auto i = 0; i < heads.size(); i++)
    {
        if ((!keep && !m.count(i)) || (keep && m.count(i)))
        {
            tmp.push_back(heads[i]);
        }
    }
    
    for (auto i = 0; i < tmp.size(); i++)
    {
        w << tmp[i] << (i != tmp.size() - 1 ? "\t" : "");
    }
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        auto started = false;
        for (auto j = 0; j < row.size(); j++)
        {
            if (!j)
            {
                w << std::endl;
            }
            
            if ((!keep && !m.count(j)) || (keep && m.count(j)))
            {
                w << (started ? "\t" : "") << row[j];
                started = true;
            }
        }
    }
    
    w.close();
}

void Anaquin::RRenameC(const FileName &src, const FileName &dst, const std::vector<Column> &cols)
{
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));    
    std::ofstream w(dst);

    for (auto j = 0; j < cols.size(); j++)
    {
        w << (j ? "\t" : "") << cols[j];
    }

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        w << std::endl;
        
        for (auto j = 0; j < row.size(); j++)
        {
            w << (j ? "\t" : "") << row[j];
        }
    }
    
    w.close();
}

void Anaquin::RFilterC(const FileName &src, const FileName &dst, const std::set<Label> &cols, bool keep)
{
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    std::set<std::size_t> m;
    
    for (auto iter = cols.begin(); iter != cols.end(); iter++)
    {
        m.insert(std::find(heads.begin(), heads.end(), *iter) - heads.begin());
    }
    
    std::ofstream w(dst);
    std::vector<Token> tmp;
    
    for (auto i = 0; i < heads.size(); i++)
    {
        if ((!keep && !m.count(i)) || (keep && m.count(i)))
        {
            tmp.push_back(heads[i]);
        }
    }
    
    for (auto i = 0; i < tmp.size(); i++)
    {
        w << tmp[i] << (i != tmp.size() - 1 ? "\t" : "");
    }
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        auto started = false;
        for (auto j = 0; j < row.size(); j++)
        {
            if (!j)
            {
                w << std::endl;
            }
            
            if ((!keep && !m.count(j)) || (keep && m.count(j)))
            {
                w << (started ? "\t" : "") << row[j];
                started = true;
            }
        }
    }
    
    w.close();
}
