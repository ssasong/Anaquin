#include "tools/tools.hpp"
#include "data/resources.hpp"

#include "resources/rna.txt"
#include "resources/meta.txt"
#include "resources/norm.txt"
#include "resources/split.txt"
#include "resources/cancer.txt"
#include "resources/anaquin.txt"
#include "resources/somatic.txt"
#include "resources/germline.txt"
#include "resources/calibrate.txt"
#include "resources/broad_bam.txt"
#include "resources/broad_vcf.txt"

#include "resources/hg382chrQ.py"
#include "resources/plotGene.R"
#include "resources/plotNorm.R"
#include "resources/plotAllele.R"
#include "resources/plotLinear.R"
#include "resources/plotInsert.R"
#include "resources/plotSomatic.R"
#include "resources/plotLDensity.R"
#include "resources/plotKSomatic.R"
#include "resources/plotLogistic.R"
#include "resources/plotKSynthetic.R"
#include "resources/plotQualFilter.R"
#include "resources/plotStrelkaROC.R"
#include "resources/plotSplitSomatic.R"

using namespace Anaquin;

static std::string clean(const std::string &x)
{
    return x.substr(0, x.find("<<@@@@>>"));
}

typedef std::string Scripts;
#define ToString(x) clean(std::string(reinterpret_cast<char*>(x)))

Scripts Manual() { return ToString(data_manuals_anaquin_txt); }

Scripts Anaquin::hg382chrQ()        { return ToString(scripts_hg382chrQ_py);   }
Scripts Anaquin::PlotGene()         { return ToString(src_r_plotGene_R);       }
Scripts Anaquin::PlotNorm()         { return ToString(src_r_plotNorm_R);       }
Scripts Anaquin::PlotLinear()       { return ToString(src_r_plotLinear_R);     }
Scripts Anaquin::PlotInsert()       { return ToString(src_r_plotInsert_R);     }
Scripts Anaquin::PlotAllele()       { return ToString(src_r_plotAllele_R);     }
Scripts Anaquin::PlotSomatic()      { return ToString(src_r_plotSomatic_R);    }
Scripts Anaquin::PlotLDensity()     { return ToString(src_r_plotLDensity_R);   }

Scripts Anaquin::PlotKSomatic()     { return ToString(src_r_plotKSomatic_R);   }
Scripts Anaquin::PlotKSynthetic()   { return ToString(src_r_plotKSynthetic_R); }
Scripts Anaquin::PlotStrelkaROC()   { return ToString(src_r_plotStrelkaROC_R); }
Scripts Anaquin::PlotQualFilter()   { return ToString(src_r_plotQualFilter_R); }
Scripts Anaquin::PlotSplitSomatic() { return ToString(src_r_plotSplitSomatic_R); }

Scripts Anaquin::rna()       { return ToString(data_manuals_rna_txt);       }
Scripts Anaquin::meta()      { return ToString(data_manuals_meta_txt);      }
Scripts Anaquin::norm()      { return ToString(data_manuals_norm_txt);      }
Scripts Anaquin::split()     { return ToString(data_manuals_split_txt);     }
Scripts Anaquin::cancer()    { return ToString(data_manuals_cancer_txt);    }
Scripts Anaquin::somatic()   { return ToString(data_manuals_somatic_txt);   }
Scripts Anaquin::germline()  { return ToString(data_manuals_germline_txt);  }
Scripts Anaquin::calibrate() { return ToString(data_manuals_calibrate_txt); }
Scripts Anaquin::broadBAM()  { return ToString(data_manuals_broad_bam_txt); }
Scripts Anaquin::broadVCF()  { return ToString(data_manuals_broad_vcf_txt); }
