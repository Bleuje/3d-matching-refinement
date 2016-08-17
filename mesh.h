#ifndef MESH_H
#define MESH_H

#pragma once

#include <bits/stdc++.h>
#include "utils.h"
using namespace std;

///"Matching 3d"
namespace m3d {

    const int MAX_SIZE_V = 7000;
    const int MAX_SIZE_T = 15000;
    const float INFINITE = 1e9;

    const int SAMPLE_SIZE = 100;
    const int SAMPLE_SIZE_LARGE = 1000;
    const int BFS_NB_HERE = 500;
    const int BFS_NB_HERE2 = 40;
    const int BFS_NB_THERE = 40;

    const float DIFF_TOLERANCE = 0.02;
    const float CORRECT_TOLERANCE = 0.01;

    const float DIFF_MAX = 0.2;
    const float PARAM_B = 0.1;

    const int CORRECT_WEIGHT = 0;

    const float LAMBDA = 0.0;
    const int MESH2_TRIES = 100;

    static float geoDistances[MAX_SIZE_V][MAX_SIZE_V][2];

    struct Mesh
    {
        int n;
        int p;
        int m;
        bool type;
        string header;
        float vert[MAX_SIZE_V][3];
        int triangles[MAX_SIZE_T][3];
        vector<int> graphMesh[MAX_SIZE_V];

        bool isNew(const int& v,const int& val);

        Mesh();
        float shortDist(const int& i1,const int& i2);
        void load(const string& path,const string& name);
        set<int> bfs(const int& start,const int& nb);
        void dijkstraFill(const int& start);
        void dijkstraFillEverything();
    };

    struct Environment{
        string pathMesh;
        string pathCorr;
        string prefix;
    };

    class MeshCorrespondence
    {
    private:
        Mesh mesh1;
        Mesh mesh2;
        int id1;
        int id2;
        Environment env;
        vector<int> samplingVec;
        map<int,int> indexS;
        map<int,int> indexSL;
        int sampleEval[SAMPLE_SIZE];
        vector<int> corr;

    public:
        MeshCorrespondence(const string& pm,const string& pc,const string& pf,const int& k1,const int& k2);
        MeshCorrespondence(const string &pm, const string &pc, const string &pf,const string& corrName, const int& k1, const int& k2);
        float globalEval();
        float globalEvalGeo();
        double rho(const double& x,const double& y);
        double evaluate(const int& p,const set<int>& here);
        double evaluateLarge(const int& p,const set<int>& here);
        void optimize(const int& nSteps);
        void findCorrect(const int& nSteps);
        void showGroundTruthError(const string& name);
        void showGeodesicError(const string& name);
        void dijkstra();
        void setSampling();
        void saveCorr();
        void loadCorr();
        void loadCorrName(const string& name);
    };

}

#endif // MESH_H
