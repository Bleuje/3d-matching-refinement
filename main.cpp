/*****
Etienne JACOB, 16/08/2016
*****
Off files and correspondences given.
Trying to find out bad correspondences.
And trying to correct them.
---
C++11 features are used
*****/


#include "mesh.h"
#include "utils.h"
using namespace std;
using namespace m3d;

///Mesh location
const string pathMesh = "C:/Users/New user/Downloads/fmaps_siggraph2012_code/shapes/tosca-blended/";
///Initial correspondences location
const string pathCorr = "C:/Users/New user/Downloads/learning/LearningFinal/html/";
///Mesh name (without the postfix index)
const string prefix = "gtfaustoff_all";

int main(){
    srand(time(0));
    ios_base::sync_with_stdio(0);

    ///Load the pair of meshes
    MeshCorrespondence MC(pathMesh,pathCorr,prefix,41,1);
    cout << "loading ok ... \n\n" << flush;

    ///Computation of geodesic distances
    cout << "Computing geodesic distances...\n" << flush;
    MC.dijkstra();
    cout << "Geodesic distances computed. \n\n" << flush;

    ///Initial evaluation from ground truth
    cout << endl << MC.globalEval() << endl << endl;
    MC.showGroundTruthError("before");

    ///Algorithm steps
    for(int k=1;k<=1000;k++){
        cout << "100 steps No." << k << "\n"<<flush;
        //MC.findCorrect(1000);
        //cout << MC.globalEval() << endl;
        MC.optimize(100);
        cout << MC.globalEval() << endl;
    }

    ///Error compared to ground-truth (colored mesh)
    MC.showGroundTruthError("after4");
    ///Estimated error using differences of geodesic distances (colored mesh)
    MC.showGeodesicError("after4");

    ///Final evaluation from ground truth
    cout << endl << MC.globalEval() << endl << endl;

    ///Saving the result
    MC.saveCorr();

    return 0;
}
