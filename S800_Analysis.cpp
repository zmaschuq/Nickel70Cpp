#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "H5Cpp.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TApplication.h"

using namespace std;

struct S800Data {
    double tof_corr;
    double de;
};

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    H5::Exception::dontPrint(); 

    ofstream S800Event("TrackEventsPID.txt", ios::trunc);

    TH2D *h2 = new TH2D("h2", "70Ni PID;TOF;Energy Loss", 1000, -2000, 2000, 1000, 0, 4000);
    
    string root_path = "/groups/tahn1/data/70Ni_NSCL/rootS800/cal/run-2020-00.root";
    string h5_path = "/groups/tahn1/data/70Ni_NSCL/h5/run_0020.h5";

    TFile *inputFile = TFile::Open(root_path.c_str(), "READ");
    TTree *tree = (TTree*)inputFile->Get("caltree");

    TLeaf *S800TS = tree->GetLeaf("fts");
    TLeaf *TotalDE = tree->GetLeaf("fIC.fsum");
    TLeaf *tacObj = tree->GetLeaf("fTOF.ftac_obj");
    TLeaf *multihit_rf = tree->GetLeaf("fMultiHitTOF.fRf");

    map<long long, S800Data> s800_map;
    int numEntries = tree->GetEntries();

    for (Long64_t i = 0; i < numEntries; i++) {
        tree->GetEntry(i);
        
        // Check to see if leaf exists, if it does get the value
        double MH_rf_val = (multihit_rf && multihit_rf->GetLen() > 0) ? multihit_rf->GetValue(0) : 0;
        double de_val  = (TotalDE) ? TotalDE->GetValue() : 0;
        double tacObj_val = (tacObj) ? tacObj->GetValue(0) : 0;

        if (tacObj_val != 0 && MH_rf_val != 0 && de_val > 500) {
            S800Data data;
            data.de = de_val;
            data.tof_corr = -MH_rf_val + tacObj_val; 
            
            s800_map[(long long)S800TS->GetValue()] = data;
        }
    }

    H5::H5File file(h5_path.c_str(), H5F_ACC_RDONLY);
    H5::Group getGroup = file.openGroup("/get");
    
    long long ts_offset = -1492; 
    int match_count = 0;

    TCutG *gate = new TCutG("gate1", 5);
    gate->SetPoint(0, 675, 1750);
    gate->SetPoint(1, 675, 2150);
    gate->SetPoint(2, 740, 2150);
    gate->SetPoint(3, 740, 1750);
    gate->SetPoint(4, 675, 1750);

    for (int i = 10759; i <= 32028; i++) {
        // Try-catch statement in the case event doesn't exist
        try {
            H5::DataSet dataset = getGroup.openDataSet("evt" + to_string(i) + "_header");
            vector<long long> hb(3);
            dataset.read(hb.data(), H5::PredType::NATIVE_LLONG);
            long long corrected_ts = hb[2] + ts_offset;

            if (s800_map.count(corrected_ts)) {
                
                h2->Fill(s800_map[corrected_ts].tof_corr, s800_map[corrected_ts].de);
                if (gate->IsInside(s800_map[corrected_ts].tof_corr, s800_map[corrected_ts].de)) {

                    S800Event << i << endl;
                }

                match_count++;
            }
        } catch (...) { continue; }
    }

    cout << "Matches Plotted: " << match_count << endl;

    TCanvas *c1 = new TCanvas("c1", "70Ni PID - Corrected", 1000, 800);
    c1->SetLogz(); 
    h2->SetStats(0);
    
    // Zoom to the data area
    if (match_count > 0) {
        h2->GetXaxis()->SetRangeUser(h2->GetMean(1)-500, h2->GetMean(1)+500);
        h2->GetYaxis()->SetRangeUser(h2->GetMean(2)-1000, h2->GetMean(2)+1000);
    }

    h2->Draw("colz");
    gate->Draw("SAME");
    app.Run();
    return 0;
}
