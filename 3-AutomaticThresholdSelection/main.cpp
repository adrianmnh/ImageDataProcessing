/* 

Created by Adrian Noa

Objective: 

To take an image's histogram from a text and find the best threshold using the bi-Gaussian automatic threshold selection method.

This process yields a Gaussian function which will be graphed along with the histogram, and the best threshold.

Usage:

g++ noa_adrian_main.cpp -o main.exe && ./main.exe BiGuass_data2.txt outFile1.txt outFile2.txt && rm ./main.exe

*/


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <string>
using namespace std;


class BiMean{

public:

    int numRows, numCols, minVal, maxVal;
    int maxHeight = 0; // largest histAry[i]
    int maxGval; // maximum calculated distribution value;
    int offSet; // one-tenth of the maxVal-minVal
/* 
        in bimodal histogram, the first modal occupies at least one-tenth of
        the histogram population from minVal to maxVal of the histogram
*/ 
    int dividePt; // Initialy set of offset, increases by 1 each iteration.
/*
        Selected treshold value is at the point at divedePt where the "distance"
        between the two bi-Gaussian curves and the histogram is the minimum
*/

    // All arrays need to be dynamically allocated at run-time
    int* histAry;       // 1D[maxVal+1] to store Histogram Array
                        
    int* GaussAry;      // 1D[maxVal+1] to store "modified" Gaussian function

    char** histGraph;   // 2D[maxVal+1 x maxHeight+1] initialize to "blank"
/*        for displaying the histogram curve.       */

    char** GaussGraph;  // 2D[maxVal+1 x maxHeight+1] initialize to "blank" 
/*        for displaying Gaussian curves in 2D.       */

    BiMean(ifstream &inFile){
    // BiMean(int numRows, int numCols,int minVal,int maxVal){

    inFile >> numRows >> numCols >> minVal >> maxVal;
        histAry = new int[maxVal+1]();
        cout << "Allocated 1D Histogram Array" << endl;
    }

    int loadHist(ifstream &inFile){ // add histogram to histAry from inFile and returns the max among hist[i]
        int numberOfPixels;
        int max = -999;
        for (int i = 0; i<maxVal+1; i++){
            inFile >> numberOfPixels >> numberOfPixels;
            histAry[i] = numberOfPixels;
            max = (numberOfPixels > max) ? max = numberOfPixels : max; 
        }
        return max;
    }

    void allocate(){

        this->GaussAry = new int[maxVal+1]{0};
        cout << "Allocated 1D Gauss array" << endl;

        this->histGraph = new char*[maxHeight+1];
        cout << "Allocated Histogram 2D char array" << endl;

        this->GaussGraph = new char*[maxHeight+1];
        cout << "Allocated Gaussian graph 2D char array" << endl;

        for(int i=0; i<maxHeight+1;i++) {
            histGraph[i] = new char[maxVal+1]();
            GaussGraph[i] = new char[maxVal+1]();
            for (int j = 0; j < maxVal+1; j++) {
                histGraph[i][j]=' ';
                GaussGraph[i][j]=' ';
            }            
        }
        cout << "Allocated histogram and gaussian subarrays" << endl;
    } 

    void plotGraph(int* ary, char** graph, char symbol){
        int j;
        for (int i = 0; i < maxVal+1; i++){
            graph[maxHeight][i] = '_';
            if(ary[i] >= 0){
                for (j = maxHeight-ary[i]; j < maxHeight; j++) {     // display
                    graph[j][i] = symbol;                            // as
                }                                                    // line graph
            }
        }
    }

    void addVertical(char** graph, int thr){
        for (int j = 0; j < maxHeight; j++) {
            graph[j][thr] = '|';
        }  
    }

    double computeMean(int leftIndex, int rightIndex, int maxHeight){
        int sum = 0;
        int numPixels = 0;
        maxHeight = 0;

        for (int i = leftIndex; i<rightIndex; i++){
            sum += histAry[i]*i;
            numPixels += histAry[i];
            if(histAry[i]>maxHeight){
                maxHeight = histAry[i];
            }
        }
        return double(sum)/(double)numPixels;
    }

    double computeVar(int leftIndex, int rightIndex, double mean){
        double sum = 0.0;
        int numPixels = 0;
        for (int i=leftIndex ; i < rightIndex ; i++){
            sum += (double)histAry[i] * pow((double)i-mean,2);
            numPixels += histAry[i];
        }
        return (double)sum/(double)numPixels;
    }

    double modifiedGauss(int x, double mean, double var, int maxHeight){
        //G(X) = maxHeight * exp( - ( (x-mean)^2 / (2* c2) )

        double G = (double)maxHeight * exp(-1 * pow(x-mean, 2) / (2*var));
        // cout << x << " " << mean << " " << var << " " << maxHeight << " " << G << endl;
        return G;
    }

    void setZero(int* ary){
        for (int i = 0; i < maxVal+1; i++) {
            ary[i]=0;
        }
    }

    int biMeanGauss(int dividePt, ofstream &outFile){
        outFile << "DividePt" << "\tLeftSum" << "\t\tRightSum" << "\tTotalSum" << "  PrevDiff " << " BestThr" << endl;
        int bestThr = dividePt;
        double sum1;
        double sum2;
        double total;
        double minSumDiff = 99999.9;
        while( dividePt < (maxVal - offSet)){
            setZero(GaussAry);
            sum1 = fitGauss(0, dividePt, GaussAry);
            sum2 = fitGauss(dividePt, maxVal, GaussAry);
            total = sum1 + sum2;
            if(total<minSumDiff) {
                minSumDiff = total;
                bestThr = dividePt;
            }
            toFile(outFile, dividePt, sum1, sum2, total, minSumDiff, bestThr);
            dividePt++;
        }
        return bestThr;
    }

    void toFile(ofstream &outFile, int dividePt, double sum1, double sum2, double total, double minSumDiff, int bestThr){
        outFile << "\t" << dividePt 
                << "\t\t" << sum1 
                << "\t\t" << sum2 
                << "\t\t" << total 
                << "\t\t" << minSumDiff 
                << "\t\t" << bestThr << endl;
    }

    double fitGauss(int leftIndex, int rightIndex, int* GaussAry){
        double mean, var, Gval, maxGval, sum = 0.0;
        mean = computeMean(leftIndex, rightIndex, maxHeight);
        var = computeVar(leftIndex, rightIndex, mean);
        for (int i = leftIndex; i <= rightIndex; i++){
            Gval = modifiedGauss(i, mean, var, maxHeight);
            sum += abs(Gval - (double)histAry[i]);
            GaussAry[i] = (int)Gval;
        }
        return sum;
    }

    void bestFitGauss(int bestThrVal){
        double sum1, sum2;
        setZero(GaussAry);
        sum1 = fitGauss(0, bestThrVal, GaussAry);
        sum2 = fitGauss(bestThrVal, maxVal, GaussAry);
    }

    void bestThr(ofstream &outFile, int bestThrVal){
        outFile << "Best Threshold value: " << bestThrVal << endl << endl;
    }

    void drawGraph(ofstream &outFile, char** graphAry, string str, int* ary){
        drawTitle(outFile, str);
        for (int i = 0; i < maxHeight+1; i++){
            for (int j = 0; j < maxVal+1; j++){
                if(i == maxHeight) outFile << graphAry[i][j] << graphAry[i][j];
                else outFile << graphAry[i][j] << " ";
            } outFile << endl;
        } outFile << endl;
    }

    void plotAll(ofstream &outFile, int bestThr){
        drawTitle(outFile, "Gaussian Curve '+' with Histogram overlay'o/-'");
        plotOverlay(outFile, bestThr);
    }

    void drawTitle(ofstream &outFile, string str){
        int l = ((maxVal*2+1)/2)-1;
        int sl = l - str.length()/2;
        for(int r=0; r<3; r++){
            for (int i = 0; i < (maxVal+1)*2; i++){
                if(r==0) outFile << "*";
                else if(r==1){
                    if(i==sl-1) outFile << "  " << str << "  ";
                    else if(i<=sl-1) outFile << "*";
                    else if (i>sl+str.length()+2) outFile << "*";
                } else outFile << "*";
            } outFile << endl;
        }
    }

    void plotOverlay(ofstream &outFile, int bestThr){

        char** overlay = allocateOverlay();

        for (int i = 0; i < maxVal+1; i++) {            // draw Gaussian Curve

            for (int j = maxHeight-GaussAry[i]; j < maxHeight; j++) overlay[j][i] = '+'; 

            for (int j = maxHeight-histAry[i]; j < maxHeight; j++){  // draw Histogram
                if(overlay[j][i] == '+') overlay[j][i] = 'O';
                else if (overlay[j][i] == ' ')overlay[j][i] = '-';
                // else overlay[j][i] = '^';
            }         
        }

        addVertical(overlay, bestThr);
        drawToFile(outFile, overlay);
        deleteOverlay(overlay);
    }

    char** allocateOverlay(){
        char **overlay = new char*[maxHeight+1];
        for(int j=0; j<maxHeight+1; j++) {
            overlay[j] = new char[maxVal+1];
            for (int i = 0; i < maxVal+1; i++){
                if(j==maxHeight) overlay[j][i] = '_';
                else overlay[j][i] = ' ';
            }
        }
        cout << "allocated 2d for overlay" << endl;
        return overlay;
    }

    void drawToFile(ofstream &outFile, char** graph){
        for (int i = 0; i < maxHeight+1; i++){
            for (int j = 0; j < maxVal+1; j++){
                if(i==maxHeight) outFile << graph[i][j] << graph[i][j];
                else outFile << graph[i][j] << " ";
            }
            outFile << endl;
        }
    }
    void deleteOverlay(char** graph){
        for (int i = 0; i < maxHeight+1; i++) delete[] graph[i];
        delete[] graph;
        cout << "freed overlay resources" << endl;
    }


    ~BiMean(){
        delete[] histAry;
        delete[] GaussAry;
        cout << "deleted histogram and gaussian dynamic arrays" << endl;
        for(int i=0; i<maxHeight+1; i++){
            delete[] histGraph[i];
            delete[] GaussGraph[i];
        }
        cout << "delete all subarrays in histgraph and gausGraph" << endl;
        delete[] histGraph;
        cout << "deleting histGraph memory allocation" << endl;
        delete[] GaussGraph;
        cout << "deleting gausGraph memory allocation";
    }

};

int main(int argc, const char* argv[]) {
    cout << endl;

    if (argc != 4){
        printf("Not enough arguments\n");
        return 1;
    }

    ifstream inFile(argv[1]);
    ofstream outFile1(argv[2]); 
    ofstream outFile2(argv[3]); 

    if (!inFile.is_open()) {
        cout << "Unable to open file" << endl;
        exit(2);
    }

    BiMean biMean = BiMean(inFile);

    // int numRows, numCols, minVal, maxVal;
    biMean.maxHeight = biMean.loadHist(inFile);

    biMean.allocate();

    // imageProcessing.shout();
    biMean.plotGraph(biMean.histAry, biMean.histGraph, '*');
    biMean.drawGraph(outFile1, biMean.histGraph, "Histogram Graph", biMean.histAry);

    biMean.offSet = (biMean.maxVal - biMean.minVal) / 10;

    biMean.dividePt = biMean.offSet;
    int bestThrVal = biMean.biMeanGauss(biMean.dividePt, outFile2);

    biMean.bestFitGauss(bestThrVal);
    biMean.plotGraph(biMean.GaussAry, biMean.GaussGraph, '+');
    biMean.drawGraph(outFile1, biMean.GaussGraph, "Gaussian Curve Graph", biMean.GaussAry);

    biMean.bestThr(outFile1, bestThrVal);

    biMean.addVertical(biMean.histGraph, bestThrVal);

    biMean.drawGraph(outFile1, biMean.histGraph, "Best Threshold Histogram", biMean.histAry);

    biMean.plotAll(outFile1, bestThrVal);


    inFile.close();
    outFile1.close();
    cout << "outFile1 graphs created successfully" << endl;
    outFile2.close();
    cout << "outFile2 debug created successfully" << endl;


    return 0;
}