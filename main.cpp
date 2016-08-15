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
const string pathMesh = "C:/Users/Etienne/Desktop/before_mincut/meshes/";
///Initial correpondances location
const string pathCorr = "C:/Users/Etienne/Desktop/before_mincut/data/";
///Mesh name (wuthout the prefix index)
const string prefix = "gtfaustoff_all";

int main(){
    srand(time(0));
    ios_base::sync_with_stdio(0);

    ///Load the pair of meshes
    MeshCorrespondance MC(pathMesh,pathCorr,prefix,38,6);
    cout << "loading ok ... \n\n" << flush;

    ///Computation of geodesic distances
    cout << "Computing geodesic distances...\n" << flush;
    MC.dijkstra();
    cout << "Geodesic distances computed. \n\n" << flush;

    ///Initial evaluation from ground truth
    cout << endl << MC.globalEval() << endl << endl;
    MC.showGroundTruthError("before");

    ///Algorithm steps
    for(int k=1;k<=100;k++){
        cout << "1000 steps No." << k << "\n"<<flush;
        MC.findGood(15000);
        cout << MC.globalEval() << endl;
        MC.optimize(1000);
        cout << MC.globalEval() << endl;
    }

    ///Error compared to ground-truth (colored mesh)
    MC.showGroundTruthError("after");
    ///Estimated error using differences of geodesic distances (colored mesh)
    MC.showGeodesicError("after");

    ///Final evaluation from ground truth
    cout << endl << MC.globalEval() << endl << endl;

    ///Saving the result
    MC.saveCorr();

    return 0;
}
