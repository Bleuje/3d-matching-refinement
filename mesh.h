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
    const int SAMPLE_SIZE = 100;
    const int SAMPLE_SIZE_LARGE = 1000;
    const float INFINITE = 1e9;
    const float VOTE_LIM = 0.02;
    const int BFS_NB = 300;

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
        vector<bool> isCorrect;
        vector<int> samplingVec;
        vector<int> samplingVecLarge;
        map<int,int> indexS;
        map<int,int> indexSL;
        int sampleEval[SAMPLE_SIZE];
        vector<int> corr;

    public:
        MeshCorrespondence(const string& pm,const string& pc,const string& pf,const int& k1,const int& k2);
        float simpleGlobalEval();
        float globalEval();
        long long evaluate(const int& p,const set<int>& here);
        long long evaluateLarge(const int& p,const set<int>& here);
        void optimize(const int& nSteps);
        void findCorrect(const int& nSteps);
        void showGroundTruthError(const string& name);
        void showGeodesicError(const string& name);
        void dijkstra();
        void setSampling();
        void saveCorr();
        void loadCorr();
    };

}

#endif // MESH_H
