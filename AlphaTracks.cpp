#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sstream>
#include "H5Cpp.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TSpline.h"

using namespace std;
#define PI 3.14159

double *extrapolate(vector<double> &, double, const double);
vector<double> g(vector<float> &, vector<double> &, float, float &, float &);

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

    ifstream myfile(RunFileName);
    if (!myfile) {
        cerr << "This File Does Not Exist" << endl;
        return 1;
    }
    myfile.close();

    ofstream outputfile("AlphaTracksHTResults.txt", ios::trunc);
    outputfile << "Event" << "  " << "Slope" << "  " << "Intercept" << "  " << "Intersections" << endl;

    // Reading in the HDF5 File
    int event = 18111; // Must change this at some point to read all events
    H5::H5File file(RunFileName, H5F_ACC_RDONLY);

    string datasetPath = "get/evt" + to_string(event) + "_data";
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

        const double drift_vel = 11.59e+6;
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
    // Plotting the original data R vs Z
    TGraph *gr = new TGraph(x_list.size(), z.data(), r_list.data());
    gr->SetTitle("AT-TPC XY Projection;Z (mm);R (mm)");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(4000.0);
    gr->GetXaxis()->SetLimits(0.0, 1800.0);
    gr->GetYaxis()->SetLimits(0.0, 275.0);

    TCanvas *c1 = new TCanvas();
    gr->Draw("AP");

    // Here must use Hough Transform to isolate particle tracks and subsequently analyze them
    float theta_high, theta_low, theta_diff, theta_increment;
    float r_increment;
    
    cout << "Input lowest theta value, highest, and increment (DEGREES/FLOATS): " << endl;
    cin >> theta_low >> theta_high >> theta_increment;
    cout << "Input r increment (mm/Float):  " << endl;
    cin >> r_increment;

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
    TH2D *h2 = new TH2D("h2", "Hough Space", thetabins, theta_low, theta_high, rbins, rhoMin, rhoMax);
    for (int i = 0; i < rbins; i++) {
        for (int j = 0; j < thetabins; j++) {

            h2->SetBinContent(j, i, accumulator[i][j]);
        }
    }
    h2->GetXaxis()->SetTitle("theta");
    h2->GetYaxis()->SetTitle("rho");

    TCanvas *c2 = new TCanvas();
    h2->Draw("colz");

    // Extracting the theta and r values that correspond to the highest voted cell

    for (int i = 0; i < rbins; i++) {
        for (int j = 0; j < thetabins; j++) {

            if (accumulator[i][j] >= 10) {

                float r_peak = rhoMin + i*r_increment;
                float theta_peak = (theta_low + j*theta_increment) * PI / 180.0;
                float slope = -cos(theta_peak) / sin(theta_peak);
                float intercept = r_peak / sin(theta_peak);
                float roundSlope = round(slope * 100.0) / 100.0;
                float roundIntercept = round(intercept * 100.0) / 100.0;

                outputfile << event << "  " << roundSlope << "  " << roundIntercept << "  " << accumulator[i][j] << endl;
            }
        }
    }
    outputfile.close();

    string inputLine;
    ifstream HTResults("AlphaTracksHTResults.txt");
    getline(HTResults, inputLine); // Skip the first line

    // Extracting data from the Hough Transform for analysis
    double HTFileSlope, HTFileIntercepts, yValueHT, xValueHT;
    int HTFileIntersections, HTFileEvent;

    // Vectors to store data from Hough Transform
    vector<double> isolated_r;
    vector<double> isolated_Q;
    vector<double> isolated_z;
    vector<double> isolated_x;
    vector<double> isolated_y;

    while (HTResults >> HTFileEvent >> HTFileSlope >> HTFileIntercepts >> HTFileIntersections) {
        double slope = HTFileSlope;
        double intercept = HTFileIntercepts;

        for (int i = 0; i < z.size(); i++) {
            double r_i = r_list[i];
            double z_i = z[i];

            double dist;
            if (isinf(slope)) {
                dist = fabs(z_i - intercept);  // vertical line
            } else {
                dist = fabs(r_i - slope * z_i - intercept) / sqrt(1 + slope*slope);
            }

            if (dist < 10.0) {
                isolated_r.push_back(r_i);
                isolated_Q.push_back(Q[i]);
                isolated_z.push_back(z_i);
                isolated_x.push_back(x_list[i]);
                isolated_y.push_back(y_list[i]);
            }
        }
    }

    // Plotting the Isolated Data
    TGraph *gr1 = new TGraph(isolated_r.size(), isolated_z.data(), isolated_r.data());
    gr1->SetTitle("Isolated Data;Z (mm);R (mm)");
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerSize(4000.0);
    gr1->GetXaxis()->SetLimits(0.0, 1800.0);
    gr1->SetMinimum(0.0);
    gr1->SetMaximum(275.0);

    TCanvas *c3 = new TCanvas();
    gr1->Draw("AP");

    float sigma_r[isolated_r.size()] = {0};
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

    double max_z = (*max_element(isolated_r.begin(), isolated_r.end()) - track_intercept) / track_slope;
    double interceptZ = -track_intercept / track_slope;

    double adjacent = max_z - interceptZ;
    double opposite = *max_element(isolated_r.begin(), isolated_r.end());
    double LabAngle = atan2(opposite, adjacent) * 180.0/PI;
    double uncertaintyAngle = 1 / (1.0 + track_slope*track_slope) * sqrt(slope_err) * 180.0/PI;

    // TF1 *LSQFit = new TF1("LSQ Fit", "[0]*x + [1]", 0, opposite);
    // LSQFit->SetParameters(track_slope, track_intercept);

    // LSQFit->SetLineStyle(1);
    // LSQFit->SetLineColor(2);
    // LSQFit->SetLineWidth(8);
    // LSQFit->Draw("SAME"); // Draw on same plot as Isolated Data

    cout << "Angle: " << LabAngle << " +/- " << uncertaintyAngle << endl;

    // Reading SRIM File for Energy Loss Data
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

    // Plot of original Bragg Curve
    // TGraph *gr2 = new TGraph(reverseEnergyLoss.size(), reverseAlphaDistance.data(), reverseEnergyLoss.data());
    // gr2->SetTitle("Energy Loss; Distance (mm); Energy Loss (MeV/mm)");
    // gr2->SetMarkerColor(kRed);
    // gr2->SetMarkerSize(400.0);

    // TCanvas *c4 = new TCanvas();
    // gr2->Draw("AP");

    // Defining knots for spline
    // double knots[] = {reverseAlphaDistance[2], reverseAlphaDistance[5], reverseAlphaDistance[7],
    // reverseAlphaDistance[10], reverseAlphaDistance[13], reverseAlphaDistance[16], reverseAlphaDistance[23], 
    // reverseAlphaDistance[28], reverseAlphaDistance[33], reverseAlphaDistance[36], reverseAlphaDistance[50]};

    // Extrapolating Bragg Curve to Zero with Linear Model
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

    // Plot of extrapolated Bragg Curve
    // TGraph *gr3 = new TGraph(reverseEnergyLoss.size(), reverseAlphaDistance.data(), reverseEnergyLoss.data());
    // gr3->SetTitle("11 MeV Bragg Curve Extrapolated; x (mm); dE/dx (MeV/mm)");
    // gr3->SetMarkerColor(kBlue);
    // gr3->SetMarkerSize(6000.0);

    // TCanvas *c5 = new TCanvas();
    // gr3->Draw("AP");
    
    // Convolving the ADC Data
    vector<double> centers;
    vector<double> convQsum;
    vector<float> step_func_range;
    vector<float> tau_range;
    float step_size = 5.0; float idx = 0.0; float idx1 = 21.75;

    while (idx < *max_element(isolated_r.begin(), isolated_r.end()) + 100.0) {
        step_func_range.push_back(idx);
        idx += 0.1;
    }

    while (idx1 < *max_element(step_func_range.begin(), step_func_range.end()) - 20.0) {
        tau_range.push_back(idx1);
        idx1 += step_size;
    }

    for (float tau : tau_range) {

        float g_start, g_end;
        vector<double> g_vector = g(step_func_range, isolated_Q, tau, g_start, g_end);
        centers.push_back((g_start + g_end) / 2.0);
        
        float sum_charge = 0;
        for (int i = 0; i < isolated_r.size(); i++) {
            if (isolated_r.at(i) >= g_start && isolated_r.at(i) <= g_end) {
                sum_charge += isolated_Q.at(i);
            }
        }
        convQsum.push_back(sum_charge);
    }

    float SigmaR[centers.size()] = {5.0};
    float SigmaQ[centers.size()] = {200.0};

    // Plotting original ADC data and Convolved Data
    // TCanvas *c6 = new TCanvas();
    // c6->Divide(2, 1);
    
    // c6->cd(1);
    // TGraph *gr4 = new TGraph(isolated_r.size(), isolated_r.data(), isolated_Q.data());
    // gr4->SetTitle("ADC Data vs R; R (mm); Q (Arbitrary Units)");
    // gr4->SetMarkerColor(kRed);
    // gr4->SetMarkerSize(7500.0);
    // gr4->GetXaxis()->SetLimits(0, 300);
    // gr4->Draw("AP");

    // c6->cd(2);
    // TGraph *gr5 = new TGraph(centers.size(), centers.data(), convQsum.data());
    // gr5->SetTitle("Convolved Data; R (mm); Q (Arbitrary Units)");
    // gr5->SetMarkerColor(kBlue);
    // gr5->SetMarkerSize(7500.0);
    // gr5->GetXaxis()->SetLimits(0, 300);
    // gr5->Draw("AP");

    float scaleFactor = sin(LabAngle*PI / 180.0);
    vector<double> AlphaX(reverseAlphaDistance);

    // Splining the Original Bragg Curve
    vector<double> s_knots = {AlphaX[0], AlphaX[2], AlphaX[5], AlphaX[7], AlphaX[10], AlphaX[13], AlphaX[16], AlphaX[23], 
    AlphaX[28], AlphaX[33], AlphaX[36], AlphaX[50], AlphaX[85], AlphaX[96], AlphaX[104], AlphaX[200], AlphaX.back()};

    vector<double> splineCoeff = {1.72986300e-03,  1.79761991e-03,  2.02189255e-03,  2.35494371e-03, 
    2.97410001e-03,  3.47673751e-03,  4.27879205e-03,  5.15550036e-03, 
    6.40442249e-03,  8.30772500e-03,  9.06274315e-03,  7.81611792e-03, 
    4.52578335e-03,  1.53405384e-03, -2.41494076e-05,  2.15095570e-05, 
    -1.58384820e-05,  5.72638308e-06, -2.60945180e-06};
    

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
