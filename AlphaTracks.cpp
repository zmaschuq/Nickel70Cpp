#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>
#include <unordered_set>
#include "H5Cpp.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"

using namespace std;
#define PI 3.14159

double *extrapolate(vector<double> &, double, const double);
vector<double> g(vector<float> &, vector<double> &, float, float &, float &);
double BasisSplines(int, int, double, const vector<double> &);
double EvaluateSpline(double, const vector<double> &, const vector<double> &, int);
vector<double> ComputingCoeff(vector<double> &, vector<double> &, vector<double> &, int);
vector<double> ConvolvedBraggModel(vector<double> &, vector<double> &, 
const vector<double> &, int, int, double);

int main() {
    
    int argc = 0;
    char *argv[] = {};
    TApplication theApp("App", &argc, argv);
    
    const int num_of_pads = 10240;
    vector<float> x_pad_mapping(num_of_pads, 0.0f);
    vector<float> y_pad_mapping(num_of_pads, 0.0f);
    vector<float> r_pad_mapping(num_of_pads, 0.0f);

    // Read Mapping File
    ifstream MappingFile("center_of_mass.txt");
    string line;

    while (getline(MappingFile, line)) {
        istringstream iss(line);
        int pad;
        float x, y;

        if (iss >> pad >> x >> y) {
            x_pad_mapping[pad] = x;
            y_pad_mapping[pad] = y;
            r_pad_mapping[pad] = sqrt(x * x + y * y);
        }
    }

    string RunFileName;
    cout << "Input run file name:" << endl;
    cin >> RunFileName; //HDF5 File Name: run_00XX.h5

    string H5FilePath;
    H5FilePath = "/groups/tahn1/data/70Ni_NSCL/h5/" + RunFileName;

    ifstream myfile(H5FilePath);
    if (!myfile) {
        cerr << "This File Does Not Exist" << endl;
        return 1;
    }
    myfile.close();

    // ofstream outputfile("AlphaTracksHTResults.txt", ios::trunc);
    // outputfile << "Event" << "  " << "Slope" << "  " << "Intercept" << "  " << "Intersections" << endl;
    ofstream Ni70_Results("70Ni_Results.txt", ios::app);
    // Ni70_Results << "Event" << "  " << "Angle" << "  " << "Recoil Energy" << "  " << "Vertex" << endl;
    Ni70_Results << "Event" << "  " << "Angle" << "  " << "Recoil Energy" << "  " << "Vertex" << "  " << 
    "Intersections" << endl;

    // Reading in Energy Loss File
    ifstream SRIMEnergyLossFile("energy_loss_He_11MeV.txt");
    double energyLoss, distanceTravel;

    vector<double> dEdx;
    vector<double> alphaDistanceTravel;
    vector<double> reverseEnergyLoss;
    vector<double> reverseAlphaDistance;

    while (SRIMEnergyLossFile >> distanceTravel >> energyLoss) {

        dEdx.push_back(energyLoss);
        alphaDistanceTravel.push_back(distanceTravel);
    }

    vector<double>::reverse_iterator rit;
    vector<double>::reverse_iterator rit1;

    // SRIM calculates energy loss FROM end of track, so want to reverse
    for (rit = dEdx.rbegin(); rit != dEdx.rend(); ++rit) {reverseEnergyLoss.push_back(*rit);}
    for (rit1 = alphaDistanceTravel.rbegin(); rit1 != alphaDistanceTravel.rend(); ++rit1) {
        reverseAlphaDistance.push_back(*rit1);
    }

    // Extrapolating the Bragg Curve with Linear Model
    double start = reverseAlphaDistance.back() + 1.0;
    double end   = reverseAlphaDistance.back() + 200.0;
    int nPoints  = 150;

    vector<double> extendAlphaDistance;
    extendAlphaDistance.reserve(nPoints);

    for (int i = 0; i < nPoints; i++) {
        double val = start + (end - start) * i / (nPoints - 1);
        val = std::round(val * 100.0) / 100.0; // round to 2 decimals
        extendAlphaDistance.push_back(val);
    }

    // Concatenate to original vector
    reverseAlphaDistance.insert(reverseAlphaDistance.end(),extendAlphaDistance.begin(),extendAlphaDistance.end());

    // Calculating new Energy Loss from Extrapolation
    double EPSlope = -6.6713232e-5; double EPIntrcpt = 0.248613273; // Calculated using Python
    double *dEdx_Extension = extrapolate(extendAlphaDistance, EPSlope, EPIntrcpt);
    
    for (int i = 0; i < 150; i++) {
        reverseEnergyLoss.push_back(dEdx_Extension[i]);
    }
    delete [] dEdx_Extension;
    dEdx_Extension = nullptr;

    // Splining the Original Bragg Curve
    int splineDeg = 3;
    vector<double> AlphaX(reverseAlphaDistance);
    vector<double> s_knots = {AlphaX[0], AlphaX[0], AlphaX[0], AlphaX[0], AlphaX[2], AlphaX[5], AlphaX[7], 
        AlphaX[10], AlphaX[13], AlphaX[16], AlphaX[23], AlphaX[28], AlphaX[33], AlphaX[36], AlphaX[50], 
        AlphaX[85], AlphaX[96], AlphaX[104], AlphaX[200], AlphaX.back(), AlphaX.back(), AlphaX.back(), AlphaX.back()};

    vector<double> splineCoeff = {1.72986300e-03,  1.79761991e-03,  2.02189255e-03,  2.35494371e-03, 
    2.97410001e-03,  3.47673751e-03,  4.27879205e-03,  5.15550036e-03, 
    6.40442249e-03,  8.30772500e-03,  9.06274315e-03,  7.81611792e-03, 
    4.52578335e-03,  1.53405384e-03, -2.41494076e-05,  2.15095570e-05, 
    -1.58384820e-05,  5.72638308e-06, -2.60945180e-06};

    auto splineFunc = [&](double *x, double *) {
        return EvaluateSpline(x[0], s_knots, splineCoeff, splineDeg);
    };

    double xmin = AlphaX.front(); double xmax = AlphaX.back();
    TF1 *fSpline = new TF1("Spline", splineFunc, xmin, xmax, 0);

    // Discretizing Bragg Curve for Convolution
    vector<double> xDiscretization;
    for (double x = 0.0; x <= *std::max_element(AlphaX.begin(), AlphaX.end()); x += 5.0) {xDiscretization.push_back(x);}
        
    vector<double> energyLossDiscret;
    for (double x : xDiscretization) {energyLossDiscret.push_back(fSpline->Eval(x));}
    delete fSpline;

    // Convolving the Bragg Curve
    vector<double> BraggCenters;
    vector<double> ConvBraggSum;
    vector<float> window_range;
    vector<float> Bragg_tau_range;
    float idxConv = 0.0; float idxConv1 = 0.0; float step_size = 5.0;

    while (idxConv < *std::max_element(xDiscretization.begin(), xDiscretization.end()) + 100.0) {
        window_range.push_back(idxConv);
        idxConv += 0.1;
    }

    while (idxConv1 < *std::max_element(window_range.begin(), window_range.end()) - 20.0) {
        Bragg_tau_range.push_back(idxConv1);
        idxConv1 += step_size;
    }

    for (float tau : Bragg_tau_range) {

        float BraggStart, BraggEnd;
        vector<double> Bragg_g_vector = g(window_range, energyLossDiscret, tau, BraggStart, BraggEnd);
        BraggCenters.push_back((BraggStart + BraggEnd) / 2.0);
        
        float BraggSum = 0;
        for (int i = 0; i < xDiscretization.size(); i++) {
            if (xDiscretization.at(i) >= BraggStart && xDiscretization.at(i) <= BraggEnd) {
                BraggSum += energyLossDiscret.at(i);
            }
        }
        ConvBraggSum.push_back(BraggSum);
    }

    // Splining Convolved Bragg Curve
    vector<double> s1_knots = {BraggCenters[0], BraggCenters[0], BraggCenters[0], BraggCenters[0], BraggCenters[175], 
    BraggCenters[450], BraggCenters[580], BraggCenters[660], BraggCenters[711], BraggCenters[717], BraggCenters[723], 
    BraggCenters[728], BraggCenters[735], BraggCenters[742], BraggCenters[748], BraggCenters[793], BraggCenters.back(), 
    BraggCenters.back(), BraggCenters.back(), BraggCenters.back()};

    vector<double> CoeffConvolvedBragg = ComputingCoeff(BraggCenters, ConvBraggSum, s1_knots, splineDeg);

    // Obtaining the Recoil Energy of Alphas
    ifstream EvXFile("EvsX_11MeV.txt");
    double EvX_Energy, EvX_Position;
    vector<double> x_R; vector<double> energy_R;
    vector<double> xR_rev; vector<double> energyR_rev;

    while (EvXFile >> EvX_Position >> EvX_Energy) {
        x_R.push_back(EvX_Position);
        energy_R.push_back(EvX_Energy);
    }

    vector<double>::reverse_iterator xRit;
    vector<double>::reverse_iterator eRit;

    for (xRit = x_R.rbegin(); xRit != x_R.rend(); ++xRit) {xR_rev.push_back(*xRit);}
    for (eRit = energy_R.rbegin(); eRit != energy_R.rend(); ++eRit) {energyR_rev.push_back(*eRit);}

    vector<double> EvsXKnots = {xR_rev[0], xR_rev[0], xR_rev[0], xR_rev[0], xR_rev[2], xR_rev[5], xR_rev[7], 
    xR_rev[10], xR_rev[13], xR_rev[16], xR_rev[23], xR_rev[28], xR_rev[33], xR_rev[36], xR_rev[50], xR_rev.back(), 
    xR_rev.back(), xR_rev.back(), xR_rev.back()};

    vector<double> EvXcoeff = ComputingCoeff(xR_rev, energyR_rev, EvsXKnots, splineDeg);
    
    auto EvXsplineFunc = [&](double *x, double *) {
        return EvaluateSpline(x[0], EvsXKnots, EvXcoeff, splineDeg);
    };

    TH2D *hEhT = new TH2D("hEhT", "Recoil Energy vs Angle", 181, 0, 180, 181, 0, 3);
    for (int i = 12000; i < 13500; i++) {
        // Reading in the HDF5 File
        int event = i; 
        H5::H5File file(H5FilePath, H5F_ACC_RDONLY);

        string datasetPath = "get/evt" + to_string(event) + "_data";

        // Check to see if event number exists
        if (!file.nameExists(datasetPath)) {
            cerr << "WARNING: Dataset not found for event " << event << ". Skipping." << endl;
            continue; 
        }

        H5::DataSet dataset = file.openDataSet(datasetPath);
        H5::DataSpace dataspace = dataset.getSpace();

        const int RANK = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[RANK];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        hsize_t n_rows = dims_out[0];
        hsize_t n_cols = dims_out[1];

        vector<double> raw_data(n_rows * n_cols);
        dataset.read(raw_data.data(), H5::PredType::NATIVE_DOUBLE);

        // Reshape into 2D vector
        vector<vector<double>> data(n_rows, vector<double>(n_cols));
        for (size_t i = 0; i < n_rows; ++i) {
            for (size_t j = 0; j < n_cols; ++j) {
                data[i][j] = raw_data[i * n_cols + j];
            }
        }

        vector<double> Q;
        vector<double> z;
        vector<double> x_list;
        vector<double> y_list;
        vector<double> r_list;

        for (const auto& row : data) {
            vector<double> filtered_values;
            
            for (const auto& val : row) {
                if (val < 5000) {filtered_values.push_back(val);}
            }

            if (filtered_values.size() < 25) {continue;}

            vector<double> first_twenty(filtered_values.begin() + 5, filtered_values.begin() + 25);
            double mean_baseline = accumulate(first_twenty.begin(), first_twenty.end(), 0.0) / first_twenty.size();

            vector<double> traces(row.begin() + 10, row.begin() + 500);
            if (traces.empty()) {continue;}

            double adc_max = *max_element(traces.begin(), traces.end());
            double adc = adc_max - mean_baseline;
            int tb_max = distance(traces.begin(), max_element(traces.begin(), traces.end()));

            const double drift_vel = 6.935e+6;
            const double frequency = 3.125e+6;
            double z_pos = drift_vel * tb_max / frequency;

            if (adc > 110 && row[4] > 0) {
                
                int cobo = static_cast<int>(row[0]);
                int asad = static_cast<int>(row[1]);
                int aget = static_cast<int>(row[2]);
                int channel = static_cast<int>(row[3]);
                int pad = static_cast<int>(row[4]);

                x_list.push_back(x_pad_mapping[pad]);
                y_list.push_back(y_pad_mapping[pad]);
                r_list.push_back(r_pad_mapping[pad]);
                Q.push_back(adc);
                z.push_back(z_pos);
            }
        }

        // Some events are empty, skip event
        if (r_list.empty()) {
            Ni70_Results << event << " " << "No Track!" << endl;
            continue;
        }
    // Plotting the original data R vs Z
    // TGraph *gr = new TGraph(x_list.size(), z.data(), r_list.data());
    // gr->SetTitle("AT-TPC XY Projection;Z (mm);R (mm)");
    // gr->SetMarkerColor(kBlue);
    // gr->SetMarkerStyle(20);
    // gr->SetMarkerSize(0.5);
    // gr->GetXaxis()->SetLimits(0.0, 1800.0);
    // gr->SetMinimum(0.0);
    // gr->SetMaximum(275.0);

    // TCanvas *c1 = new TCanvas();
    // gr->Draw("AP");

        // Here must use Hough Transform to isolate particle tracks and subsequently analyze them
        float theta_high = 90.0; float theta_low = -90.0; float theta_increment = 0.5;
        float r_increment = 2.5;
        float theta_diff;
    
        // cout << "Input lowest theta value, highest, and increment (DEGREES/FLOATS): " << endl;
        // cin >> theta_low >> theta_high >> theta_increment;
        // cout << "Input r increment (mm/Float):  " << endl;
        // cin >> r_increment;

        theta_diff = theta_high - theta_low;
        int thetabins = round(theta_diff/theta_increment) + 1.0;

        const float AT_TPC_Radius = 275.0;
        const float AT_TPC_Length = 1900.0;
        const float rhoBounds = sqrt(AT_TPC_Radius*AT_TPC_Radius + AT_TPC_Length*AT_TPC_Length);
        const float rhoMin = -rhoBounds, rhoMax = rhoBounds, rhoRange = rhoMax - rhoMin;
        int rbins = round(rhoRange/r_increment) + 1.0;
        
        vector<vector<int>> accumulator(rbins, vector<int> (thetabins));

        for (int i = 0; i < r_list.size(); i++) {
            for (float j = theta_low; j <= theta_high; j+=theta_increment) {

                float radians = j * PI/180.0;
                float rho = z[i]*cos(radians) + r_list[i]*sin(radians);
                
                int rIndex = round((rho - rhoMin)/r_increment);
                int thetaIndex = round((j - theta_low)/theta_increment);
                //if (rIndex >= 0 && rIndex < rbins) {accumulator[rIndex][thetaIndex]++;}
                accumulator[rIndex][thetaIndex]++;
            }
        }
    // Plotting 2D Histogram of Hough Transform
    // TH2D *h2 = new TH2D("h2", "Hough Space", thetabins, theta_low, theta_high, rbins, rhoMin, rhoMax);
    // for (int i = 0; i < rbins; i++) {
    //     for (int j = 0; j < thetabins; j++) {

    //         h2->SetBinContent(j, i, accumulator[i][j]);
    //     }
    // }
    // h2->GetXaxis()->SetTitle("theta");
    // h2->GetYaxis()->SetTitle("rho");

    // TCanvas *c2 = new TCanvas();
    // h2->Draw("colz");

        // Extracting the theta and r values that correspond to the highest voted cell
        int threshold = 20;
        vector<double> slopesHT; vector<double> interceptHT; vector<int> intersectionsHT;
        for (int i = 0; i < rbins; i++) {
            for (int j = 0; j < thetabins; j++) {

                if (accumulator[i][j] >= threshold) {

                    float r_peak = rhoMin + i*r_increment;
                    float theta_peak = (theta_low + j*theta_increment) * PI / 180.0;
                    float slope = -cos(theta_peak) / sin(theta_peak);
                    float intercept = r_peak / sin(theta_peak);
                    float roundSlope = std::round(slope * 100.0) / 100.0;
                    float roundIntercept = std::round(intercept * 100.0) / 100.0;

                    //if (roundSlope > 0.0) {
                    slopesHT.push_back(roundSlope);
                    interceptHT.push_back(roundIntercept);
                    intersectionsHT.push_back(accumulator[i][j]);
                    //} 
                } 
            }
        }

        // Vectors to store data from Hough Transform
        vector<double> isolated_r;
        vector<double> isolated_Q;
        vector<double> isolated_z;
        vector<double> isolated_x;
        vector<double> isolated_y;

        int max_intersections = -1;
        double best_slope = 0.0; double best_intercept = 0.0;

        for (int i = 0; i < slopesHT.size(); i++) {
            if (intersectionsHT[i] > max_intersections) {

                max_intersections = intersectionsHT[i];
                best_slope = slopesHT[i];
                best_intercept = interceptHT[i];
            }
        }
        
        double perpDist = 7.5;
        double slope = best_slope;
        double intercept = best_intercept;

        for (int i = 0; i < z.size(); i++) {
            double r_i = r_list[i];
            double z_i = z[i];

            double dist;
            if (isinf(slope)) {
                dist = fabs(z_i - intercept);  // vertical line
            } else {
                dist = fabs(r_i - slope * z_i - intercept) / sqrt(1 + slope*slope);
            }

            if (dist < perpDist) {
                isolated_r.push_back(r_i);
                isolated_Q.push_back(Q[i]);
                isolated_z.push_back(z_i);
                isolated_x.push_back(x_list[i]);
                isolated_y.push_back(y_list[i]);
            }
        }

        // Condition that does not meet definition of line (at least 20 points)
        if (isolated_r.empty()) {
            Ni70_Results << event << " " << "No Track!" << endl;
            continue;
        }

    // Plotting the Isolated Data
    // TGraph *gr1 = new TGraph(isolated_r.size(), isolated_z.data(), isolated_r.data());
    // gr1->SetTitle("Isolated Data;Z (mm);R (mm)");
    // gr1->SetMarkerColor(kBlue);
    // gr1->SetMarkerStyle(20);
    // gr1->SetMarkerSize(0.5);
    // gr1->GetXaxis()->SetLimits(0.0, 1800.0);
    // gr1->SetMinimum(0.0);
    // gr1->SetMaximum(275.0);

    // TCanvas *c3 = new TCanvas();
    // gr1->Draw("AP");

        vector<float> sigma_r(isolated_r.size(), 0.0f);
        for (int i = 0; i < isolated_r.size(); i++) {
            if (isolated_r[i] < 150) {sigma_r[i] = 2.5;}
            else {sigma_r[i] = 5.0;}
        }

        // Performin LSQ Fit to Extract Scattering Angle
        double A, B, C, D, E, F;
        A = 0.0; B = 0.0; C = 0.0; D = 0.0; E = 0.0; F = 0.0;

        for (int i = 0; i < isolated_r.size(); i++) {

            A += isolated_z.at(i) / (sigma_r[i]*sigma_r[i]);
            B += 1 / (sigma_r[i]*sigma_r[i]);
            C += isolated_r.at(i) / (sigma_r[i]*sigma_r[i]);
            D += (isolated_z.at(i)*isolated_z.at(i)) / (sigma_r[i]*sigma_r[i]);
            E += (isolated_z.at(i)*isolated_r.at(i)) / (sigma_r[i]*sigma_r[i]);
            F += (isolated_r.at(i)*isolated_r.at(i)) / (sigma_r[i]*sigma_r[i]);
        }

        double track_slope = (E*B - C*A) / (D*B - A*A);
        double track_intercept = (D*C - E*A) / (D*B - A*A);
        double slope_err = B / (B*D - A*A);
        double intercept_err = D / (B*D - A*A);
        double covar = -A / (B*D - A*A);

        double max_z = (*std::max_element(isolated_r.begin(), isolated_r.end()) - track_intercept) / track_slope;
        double interceptZ = -track_intercept / track_slope;

        double adjacent = max_z - interceptZ;
        double opposite = *std::max_element(isolated_r.begin(), isolated_r.end());
        double LabAngle = atan2(opposite, adjacent) * 180.0/PI;
        double uncertaintyAngle = 1 / (1.0 + track_slope*track_slope) * sqrt(slope_err) * 180.0/PI;

    // TF1 *LSQFit = new TF1("LSQ Fit", "[0]*x + [1]", 0, 1900.0);
    // LSQFit->SetParameters(track_slope, track_intercept);

    // LSQFit->SetLineStyle(1);
    // LSQFit->SetLineColor(2);
    // LSQFit->SetLineWidth(2);
    // LSQFit->Draw("L SAME"); // Draw on same plot as Isolated Data

    // cout << "Angle: " << LabAngle << " +/- " << uncertaintyAngle << endl;

    // Plot of original Bragg Curve
    // TGraph *gr2 = new TGraph(reverseEnergyLoss.size(), reverseAlphaDistance.data(), reverseEnergyLoss.data());
    // gr2->SetTitle("Energy Loss; Distance (mm); Energy Loss (MeV/mm)");
    // gr2->SetMarkerColor(kRed);
    // gr2->SetMarkerSize(400.0);

    // TCanvas *c4 = new TCanvas();
    // gr2->Draw("AP");
    
        // Convolving the ADC Data
        vector<double> centers;
        vector<double> convQsum;
        vector<float> step_func_range;
        vector<float> tau_range;

        float idx = 0.0; float idx1 = 21.75;
        while (idx < *std::max_element(isolated_r.begin(), isolated_r.end()) + 100.0) {
            step_func_range.push_back(idx);
            idx += 0.1;
        }

        while (idx1 < *std::max_element(step_func_range.begin(), step_func_range.end()) - 20.0) {
            tau_range.push_back(idx1);
            idx1 += step_size;
        }

        for (float tau : tau_range) {

            float g_start, g_end;
            vector<double> g_vector = g(step_func_range, isolated_Q, tau, g_start, g_end);
            centers.push_back(round((g_start + g_end) * 1000.0 / 2.0) / 1000.0);
            
            float sum_charge = 0;
            for (int i = 0; i < isolated_r.size(); i++) {
                if (isolated_r.at(i) >= g_start && isolated_r.at(i) <= g_end) {
                    sum_charge += isolated_Q.at(i);
                }
            }
            convQsum.push_back(sum_charge);
        }

        vector<float> SigmaR(centers.size(), 5.0f);
        vector<float> SigmaQ(centers.size(), 200.0f);

    // Plotting original ADC data and Convolved Data
    // TCanvas *c5 = new TCanvas();
    // c5->Divide(2, 1);
    
    // c5->cd(1);
    // TGraph *gr4 = new TGraph(isolated_r.size(), isolated_r.data(), isolated_Q.data());
    // gr4->SetTitle("ADC Data vs R; R (mm); Q (Arbitrary Units)");
    // gr4->SetMarkerColor(kRed);
    // gr4->SetMarkerStyle(8);
    // gr4->SetMarkerSize(1.0);
    // gr4->GetXaxis()->SetLimits(0, 300);
    // gr4->Draw("AP");

    // c5->cd(2);
    // TGraph *gr5 = new TGraph(centers.size(), centers.data(), convQsum.data());
    // gr5->SetTitle("Convolved Data; R (mm); Q (Arbitrary Units)");
    // gr5->SetMarkerColor(kBlue);
    // gr5->SetMarkerStyle(8);
    // gr5->SetMarkerSize(1.0);
    // gr5->GetXaxis()->SetLimits(0, 300);
    // gr5->Draw("AP");

        float scaleFactor = sin(LabAngle*PI / 180.0);

    // Plot of extrapolated Bragg Curve
    // fSpline->SetLineColor(kRed);
    // fSpline->SetLineWidth(1);
    // fSpline->SetNpx(1500);
    
    // TGraph *gr3 = new TGraph(reverseEnergyLoss.size(), reverseAlphaDistance.data(), reverseEnergyLoss.data());
    // gr3->SetTitle("11 MeV Bragg Curve Extrapolated; x (mm); dE/dx (MeV/mm)");
    // gr3->SetMarkerColor(kBlue);
    // gr3->SetMarkerSize(1.0);
    // gr3->SetMarkerStyle(8);
    // //gr3->GetXaxis()->SetLimits(3000, 3500);

    // TCanvas *c6 = new TCanvas();
    // gr3->Draw("AP");
    // fSpline->Draw("L");

    // TLegend *leg = new TLegend();
    // leg->AddEntry(gr3, "Bragg Curve", "p");
    // leg->AddEntry(fSpline, "Spline", "l");
    // leg->Draw();

    // TGraph *gr6 = new TGraph(ConvBraggSum.size(), BraggCenters.data(), ConvBraggSum.data());
    // gr6->SetMarkerColor(kRed);
    // gr6->SetMarkerSize(1.0);
    // gr6->SetMarkerStyle(8);

    // TCanvas *c7 = new TCanvas();
    // gr6->Draw("AP");

        vector<double> xData(centers);

        vector<int> shiftRange; vector<int> scaleYRange;
        shiftRange.clear(); scaleYRange.clear();

        for (int i = 3400; i < 3707; i++) {shiftRange.push_back(i);}
        for (int i = 100000; i < 500000; i += 10000) {scaleYRange.push_back(i);}

        vector<double> yData(convQsum);
        vector<float> QError; QError.reserve(yData.size());
                for (int i = 0; i < yData.size(); i++) {QError.push_back(SigmaQ[i]);}

        vector<vector<double>> chiSquaredMatrix(shiftRange.size(), vector<double>(scaleYRange.size(), NAN));

        for (int idxShift = 0; idxShift < shiftRange.size(); idxShift++) {
            
            vector<double> yModel = ConvolvedBraggModel(xData, s1_knots, CoeffConvolvedBragg, 
                    splineDeg, shiftRange[idxShift], scaleFactor);

            for (int idxScale = 0; idxScale < scaleYRange.size(); idxScale++) {

                double chi2 = 0.0;
                for (int k = 0; k < yData.size(); k++) {
                    
                    double diff = yData[k] - yModel[k]*scaleYRange[idxScale];
                    chi2 += diff*diff / (QError[k]*QError[k]);
                }
                chiSquaredMatrix[idxShift][idxScale] = chi2;
            }
        }

        // Obtaining best shift and scale value from minimum chi-squared
        double maxChi = std::numeric_limits<double>::infinity();
        int shiftIdx = -1; int scaleIdx = -1;

        for (int i = 0; i < chiSquaredMatrix.size(); i++) {
            for (int j = 0; j < chiSquaredMatrix[i].size(); j++) {

                double chiValue = chiSquaredMatrix[i][j];
                if (!std::isnan(chiValue) && chiValue < maxChi) {

                    maxChi = chiValue;
                    shiftIdx = i;
                    scaleIdx = j;
                }
            }
        }

        int bestShift = shiftRange[shiftIdx]; int bestScale = scaleYRange[scaleIdx];

        // cout << "Chi Value: " << maxChi << endl;
        // cout << "Shift: " << bestShift << " Scale: " << bestScale << endl;

    // Plotting fitted Convolved Bragg Curve to Convolved Data
    // vector<double> yBraggPlot = ConvolvedBraggModel(centers, s1_knots, CoeffConvolvedBragg, bestScale, splineDeg, 
    //     bestShift, scaleFactor);

    // TGraph *gr7 = new TGraph(centers.size(), centers.data(), yBraggPlot.data());
    // gr7->GetXaxis()->SetTitle("R (mm)");
    // gr7->GetYaxis()->SetTitle("Q (Arbitrary)");
    // gr7->GetXaxis()->SetLimits(0, 300);
    // gr7->SetLineStyle(10);
    // gr7->SetLineColor(2);
    // gr7->SetLineWidth(1);

    // TGraph *gr8 = new TGraph(convQsum.size(), centers.data(), convQsum.data());
    // gr8->SetMarkerStyle(20);
    // gr8->SetMarkerColor(kBlue);
    // gr8->SetMarkerSize(0.5);

    // TCanvas *c8 = new TCanvas();
    // gr7->Draw("AL");
    // gr8->Draw("P SAME");

        // Obtain Recoil Energy and Store in File
        double xRmin = xR_rev.front(); double xRmax = xR_rev.back();
        TF1 *EvXSpline = new TF1("Spline", EvXsplineFunc, xRmin, xRmax, 0);

        float EnergyofAlpha = EvXSpline->Eval(bestShift);
        // Ni70_Results << event << "  " << LabAngle << "  " << EnergyofAlpha << "  " << interceptZ << endl;
        Ni70_Results << event << "  " << LabAngle << "  " << EnergyofAlpha << "  " << interceptZ << "  " << 
        max_intersections << endl;

        delete EvXSpline;

        hEhT->Fill(LabAngle, EnergyofAlpha);
    
    // Plotting Splined Energy Curve
    // EvXSpline->SetLineColor(kRed);
    // EvXSpline->SetLineWidth(1);
    // EvXSpline->SetNpx(1500);
    
    // TGraph *gr9 = new TGraph(xR_rev.size(), xR_rev.data(), energyR_rev.data());
    // gr9->SetTitle("11 MeV Energy Curve ; x (mm); Energy (MeV)");
    // gr9->SetMarkerColor(kBlue);
    // gr9->SetMarkerSize(1.0);
    // gr9->SetMarkerStyle(8);
    // //gr3->GetXaxis()->SetLimits(3000, 3500);

    // TCanvas *c9 = new TCanvas();
    // gr9->Draw("AP");
    // EvXSpline->Draw("L SAME");
    }

    hEhT->GetXaxis()->SetTitle("Lab Angle");
    hEhT->GetYaxis()->SetTitle("Recoil Energy");
    TCanvas *c10 = new TCanvas();
    hEhT->Draw();

    theApp.Run();
    return 0;
}

double *extrapolate(vector<double> &x, double a, const double b) {
    
    double result;
    double *reverseAlphaExtension = new double [x.size()];

    for (int i = 0; i < x.size(); i++) {
        result = x.at(i)*a + b;
        result = std::max(result, 0.0); // y-values that are negative, keep at 0
        reverseAlphaExtension[i] = result;
    }

    return reverseAlphaExtension;
}

vector<double> g(vector<float> &t, vector<double> &Q, float tau, float &start_step, float &end_step) {

    int size = t.size();
    vector<double> convolutionVector (size, 0.0);

    float step_func_width = 20.0;
    start_step = tau;
    end_step = start_step + step_func_width;

    double maxQ = *max_element(Q.begin(), Q.end());
    for (int i = 0; i < size; i++) {
        if (t[i] >= start_step && t[i] <= end_step) {
            convolutionVector.at(i) = maxQ + 50.0;
        }
    }
    return convolutionVector;
}

double BasisSplines(int i, int deg, double BraggPos, const vector<double> &knots) {

    // Implementing Cox de Boor recursion
    if (deg == 0) {
        if (i == knots.size() - 2 && BraggPos == knots.back()) return 1.0;
        return (BraggPos >= knots[i] && BraggPos < knots[i+1]) ? 1.0 : 0.0;
    }

    double denom1 = knots[i + deg] - knots[i];
    double denom2 = knots[i + deg + 1] - knots[i + 1];
    double inv1 = (denom1 != 0.0) ? 1.0 / denom1 : 0.0;
    double inv2 = (denom2 != 0.0) ? 1.0 / denom2 : 0.0;
    
    double numerator1 = BraggPos - knots[i];
    double numerator2 = knots[i + deg + 1] - BraggPos;

    return numerator1 * inv1 * BasisSplines(i, deg - 1, BraggPos, knots) + 
    numerator2 * inv2 * BasisSplines(i + 1, deg - 1, BraggPos, knots);
}

double EvaluateSpline(double BraggPos, const vector<double> &knots, const vector<double> &coeff, int deg) {

    double sum = 0.0;
    BraggPos = std::max(knots.front(), std::min(BraggPos, knots.back()));

    int nBasis = knots.size() - deg - 1;
    for (int i = 0; i < nBasis; i++) {sum += coeff.at(i) * BasisSplines(i, deg, BraggPos, knots);}

    return sum;
}

vector<double> ComputingCoeff(vector<double> &xDisc, vector<double> &ELossDisc, vector<double> &knots, int deg) {
    
    // Constructing design matrix to compute coefficients (c): Ac = y
    int NumBasis = knots.size() - deg - 1;
    int NumData = xDisc.size();

    TMatrixD DesignMtx(NumData, NumBasis);
    for (int i = 0; i < NumData; i++) {
        for (int j = 0; j < NumBasis; j++) {

            DesignMtx(i, j) = BasisSplines(j, deg, xDisc[i], knots);
        }
    }

    TVectorD y(NumData);
    for (int i = 0; i < NumData; i++) {y(i) = ELossDisc[i];}

    TDecompSVD svd(DesignMtx);
    TVectorD Coeff(y);
    Bool_t ok = svd.Solve(Coeff);

    vector<double> Coefficients;
    for (int i = 0; i < NumBasis; i++) {Coefficients.push_back(Coeff(i));}

    return Coefficients;
}

vector<double> ConvolvedBraggModel(vector<double> &x,
vector<double> &knots, const vector<double> &SplineCoeff, int deg, int ParamShift, double angScale) {

    vector<double> evalSpline;
    for (int i = 0; i < x.size(); i++) {

        double xOriginal = x.at(i) / angScale + ParamShift;
        double evaluation = EvaluateSpline(xOriginal, knots, SplineCoeff, deg);
        evalSpline.push_back(evaluation);
    }

    return evalSpline;
}