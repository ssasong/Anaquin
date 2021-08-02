#include "sequins/genomics/g_calibrate.hpp"

using namespace Anaquin;

typedef GCalibrate::Options Options;

GCalibrate::Stats GCalibrate::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    GBroadBam::Options o_(cloneO(o)); o_.index = o.index; o_.debug = true; o_.showGen = false;
    o_.meth  = o.meth;
    o_.origW = o.work;
    o_.customSequinThreshold = o.customSequinThreshold;
    
    if (o.debug)
    {
        switch (o_.meth)
        {
            case CalibrateMethod::Custom:  { o.logInfo("Custom");  break; }
            case CalibrateMethod::None:    { o.logInfo("None");    break; }
            case CalibrateMethod::Mean:    { o.logInfo("Mean");    break; }
            case CalibrateMethod::Median:  { o.logInfo("Median");  break; }
            case CalibrateMethod::Percent: { o.logInfo("Percent"); break; }
        }
    }
    
    const auto isCancer = o.isCancer;

    // Generate outputs like BroadBAM (we're going to move files later)
    GBroadBam::report(f1, f2, o_);
    
    auto m1 = [&](const FileName &src, const FileName &dst)
    {
        mv(o_.work + "/" + src, o.work + "/" + dst);
    };

    auto m2 = [&](const FileName &src, const FileName &dst)
    {
        mv(o.work + "/" + src, o.work + "/" + dst);
    };

    if (o.writeS) { m1("sample.bam", "sample.bam"); }
    if (o.writeD) { m1("trimmed.bam", "sequin.bam"); } // Not from the untrimmed "sequin.bam"
    if (o.writeC) { m1("sequin_calibrated.bam", "calibrated.bam"); }

    const auto both = !f1.empty() && !f2.empty();
    
    /*
     * Merged alignments include both sample and sequin calibrated reads. However, it would make
     * only sense for a decoy protocol, where the sequin reads are not aligned to the human
     * chromosomes.
     */
    
    if (!both)
    {
        m1("merged.bam", "merged.bam");
    }
    
    const std::string prefix = o.isCancer ? "cancer" : "calibrate";
    m1("broad_bam.txt", prefix + "_report.txt");

    if (o.debug)
    {
        m1("broadBAM_regions.tsv",   prefix + "_regions.tsv");
        m1("broadBAM_variants.tsv",  prefix + "_variants.tsv");
        m1("broadBAM_synthetic.tsv", prefix + "_synthetic.tsv");
        m1("broadBAM_errors.tsv",    prefix + "_errors.tsv");
        m1("broadBAM_features.tsv",  prefix + "_features.tsv");
    }

    return GCalibrate::Stats();
}

void GCalibrate::report(const FileName &file, const Options &o)
{
    analyze(file, "", o);
}

void GCalibrate::report(const FileName &f1, const FileName &f2, const Options &o)
{
    analyze(f1, f2, o);
}
