#!/usr/bin/env python3
"""Exercise the native correction consumer with an explicit sparse ROOT payload."""
import array
import pathlib
import shlex
import subprocess
import tempfile
import unittest


class Run3CorrectorTest(unittest.TestCase):
    def test_native_inverse_fallback_and_missing_support(self):
        import ROOT
        repo = pathlib.Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory(prefix="run3-pair-loader-") as tmp:
            tmp = pathlib.Path(tmp)
            path = tmp / "table.root"
            f = ROOT.TFile(str(path), "RECREATE")
            ROOT.TObjString("run3-an751-v1").Write("Run3PairCorrectionVersion")
            d = f.mkdir("JCard")
            d.cd()
            for name, values in {"DataType": [10], "CentralityBinEdges": [0,100],
                                 "TrackPairPtBinEdges": [2,500], "JetPtBinEdgesEEC": [40,160]}.items():
                v = ROOT.TVectorT("float")(1, len(values))
                for i, value in enumerate(values, 1):
                    v[i] = value
                v.Write(name)
            f.cd()
            edges = array.array("d", [0,.1,.4,.8])
            for name, values in {
                "trackPairEfficiencyCorrection_C0T0A0": [1/3, float("nan"), 1/4],
                "trackPairEfficiencyCorrectionCloseToJet_C0T0A0J0": [.5,float("nan"),float("nan")]
            }.items():
                h = ROOT.TH1D(name, name, 3, edges)
                for i, value in enumerate(values, 1):
                    h.SetBinContent(i, value)
                h.Write()
            f.Close()
            source = tmp / "check.cxx"
            source.write_text(r'''#include "TrackPairEfficiencyCorrector.h"
#include <cmath>
#include <stdexcept>
#include <limits>
int main(int argc, char** argv){
 TrackPairEfficiencyCorrector c(argv[1],false,false);
 if(!c.IsRun3Table()) return 1;
 if(c.GetTrackPairEfficiencyCorrection(.05,20,3,4,60).first != 2) return 2;
 if(c.GetTrackPairEfficiencyCorrection(.05,20,4,3,-1).first != 3) return 3;
 if(c.GetTrackPairEfficiencyCorrection(.5,20,3,4,60).first != 4) return 4;
 for(int test=0; test<5; ++test){
   try {
     c.GetTrackPairEfficiencyCorrection(test==0?.2:test==1?.9:.05,
       test==2?101:20, test==3?1:3, 4, test==4?170:60);
     return 10+test;
   } catch(const std::out_of_range&) {}
 }
 try {c.GetTrackPairEfficiencyCorrection(.05,20,3,std::numeric_limits<double>::quiet_NaN(),60); return 20;}
 catch(const std::out_of_range&) {}
 c.SetDisableCorrection(true);
 if(c.GetTrackPairEfficiencyCorrection(.2,20,3,4,60).first != 1) return 21;
 return 0;
}
''')
            flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], text=True))
            command = ["g++", "-O2", "-I"+str(repo/"src"), str(source),
                       str(repo/"src/TrackPairEfficiencyCorrector.cxx"), str(repo/"src/BinningCard.cxx"),
                       *flags, "-o", str(tmp/"check")]
            subprocess.run(command, check=True, capture_output=True, text=True)
            subprocess.run([str(tmp/"check"), str(path)], check=True, capture_output=True, text=True)


if __name__ == "__main__":
    unittest.main()
