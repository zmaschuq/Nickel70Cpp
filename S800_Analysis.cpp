#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "H5Cpp.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TApplication.h"

using namespace std;

struct S800Data {
    float tof_corr;
    float de;
};

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    H5::Exception::dontPrint(); 

    ofstream S800Event("TrackEventsPID.txt", ios::trunc);

    TH2F *h2 = new TH2F("h2", "68Ni PID Setting: MultiHit;TOF (Obj - E1Up);Energy Loss", 500, -200, 0, 1000, 1000, 3000);

    string rootfile, h5filename;

    cout << "Input root file name: " << endl;
    cin >> rootfile;
    cout << "Input HDF5 File Name:" << endl;
    cin >> h5filename;
    
    string root_path = "/groups/tahn1/data/70Ni_NSCL/rootS800/cal/" + rootfile;
    string h5_path = "/groups/tahn1/data/70Ni_NSCL/h5/" + h5filename;

    TFile *inputFile = TFile::Open(root_path.c_str(), "READ");
    TTree *tree = (TTree*)inputFile->Get("caltree");

    TLeaf *S800TS = tree->GetLeaf("fts");
    TLeaf *TotalDE = tree->GetLeaf("fIC.fsum");
    TLeaf *ScintObj = tree->GetLeaf("fMultiHitTOF.fObj");
    TLeaf *ScintE1U = tree->GetLeaf("fMultiHitTOF.fE1Up");

    vector<long long> S800TSVec;
    map<long long, S800Data> s800_map;

    int numEntries = tree->GetEntries();
    S800TSVec.reserve(numEntries);

    for (Long64_t i = 0; i < numEntries; i++) {
        tree->GetEntry(i);
        
        float de_val  = TotalDE->GetValue();
        float ScintVal1 = ScintE1U->GetValue(0);
        float ScintVal2 = ScintObj->GetValue(0);
        long long S800TimeStamp = S800TS->GetValue();

        S800Data data;
        data.de = de_val;
        data.tof_corr = ScintVal2; 

        S800TSVec.push_back(S800TimeStamp);
            
        s800_map[(long long)S800TimeStamp] = data;
        h2->Fill(s800_map[S800TimeStamp].tof_corr, s800_map[S800TimeStamp].de);
    }

    sort(S800TSVec.begin(), S800TSVec.end());

    H5::H5File file(h5_path.c_str(), H5F_ACC_RDONLY);
    H5::Group getGroup = file.openGroup("/get");

    TCutG *gate = new TCutG("gate1", 5);
    gate->SetPoint(0, -47.5, 2200);
    gate->SetPoint(1, -39, 2200);
    gate->SetPoint(2, -39, 1850);
    gate->SetPoint(3, -47.5, 1950);
    gate->SetPoint(4, -47.5, 2200);

    for (int i = 10759; i <= 32028; i++) {
        try {
            H5::DataSet dataset = getGroup.openDataSet("evt" + to_string(i) + "_header");
            vector<long long> hb(3);
            dataset.read(hb.data(), H5::PredType::NATIVE_LLONG);
            long long ATTPCTS = hb[2];

            // BINARY SEARCH: Find the closest S800 timestamp
            auto it = lower_bound(S800TSVec.begin(), S800TSVec.end(), ATTPCTS);

            long long minDelta = -1;

            if (it == S800TSVec.begin()) {
                minDelta = abs(ATTPCTS - *it);
            } else if (it == S800TSVec.end()) {
                minDelta = abs(ATTPCTS - *(it - 1));
            } else {
                // It's between *it and *(it-1), find which is closer
                minDelta = min(abs(ATTPCTS - *it), abs(ATTPCTS - *(it - 1)));
            }
            
            long long TS_Offset = ATTPCTS - minDelta;
            if (gate->IsInside(s800_map[TS_Offset].tof_corr, s800_map[TS_Offset].de)) {
                S800Event << i << endl;
            }

        } catch (...) { continue; }
    }

    TCanvas *c1 = new TCanvas();
    c1->SetLogz(); 
    //h2->SetStats(0);
    
    // Zoom to the data area
    h2->GetXaxis()->SetRangeUser(h2->GetMean(1)-150, h2->GetMean(1)+150);
    h2->GetYaxis()->SetRangeUser(h2->GetMean(2)-1000, h2->GetMean(2)+1000);
    
    h2->Draw("colz");
    gate->Draw("SAME");
    app.Run();
    return 0;
}
