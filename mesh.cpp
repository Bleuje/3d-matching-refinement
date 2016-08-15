#include "mesh.h"

namespace m3d {

    Mesh::Mesh()
    {

    }

    float Mesh::shortDist(const int& i1,const int& i2){
        return sqrt((vert[i1][0]-vert[i2][0])*(vert[i1][0]-vert[i2][0]) + (vert[i1][1]-vert[i2][1])*(vert[i1][1]-vert[i2][1]) + (vert[i1][2]-vert[i2][2])*(vert[i1][2]-vert[i2][2]));
    }

    bool Mesh::isNew(const int& v,const int& val){
        for(int i=0;i<int(graphMesh[v].size());i++){
            if(graphMesh[v][i]==val){
                return false;
            }
        }
        return true;
    }


    void Mesh::load(const string& path,const string& name){
        string filepath = path + name + ".off";
        ifstream in (filepath);

        in>>header;

        in>>n>>p>>m;
        for(int i=0;i<n;i++){
            for(int j=0;j<3;j++){
                in>>vert[i][j];
            }
        }
        for(int i=0;i<p;i++){
            int waste;in>>waste;
            int v1,v2,v3;in>>v1>>v2>>v3;
            if(isNew(v1,v2)) graphMesh[v1].push_back(v2);
            if(isNew(v2,v1)) graphMesh[v2].push_back(v1);
            if(isNew(v3,v2)) graphMesh[v3].push_back(v2);
            if(isNew(v2,v3)) graphMesh[v2].push_back(v3);
            if(isNew(v1,v3)) graphMesh[v1].push_back(v3);
            if(isNew(v3,v1)) graphMesh[v3].push_back(v1);
            triangles[i][0]=v1;
            triangles[i][1]=v2;
            triangles[i][2]=v3;
        }
    }

    inline set<int> Mesh::bfs(const int& start,const int& nb){
        set<int> res;
        vector<bool> visited(n,0);

        queue<int> frontier;
        frontier.push(start);
        visited[start]=true;

        while(!frontier.empty() && int(res.size())<nb){
            int next = frontier.front();
            for(int i = 0;i < int(graphMesh[next].size());i++){
                int k = graphMesh[next][i];
                if(!visited[k]){
                    frontier.push(k);
                    visited[k]=true;
                }
            }
            res.insert(next);
            frontier.pop();
        }
        return res;
    }

    void Mesh::dijkstraFill(const int& start){
        for(int i=0;i<n;i++){
            geoDistances[start][i][type] = INFINITE;
        }
        geoDistances[start][start][type] = 0;

        vector<bool> alreadyComputed(n,0);

        priority_queue<pair<float,int> > frontier;
        frontier.push(make_pair(0,start));

        while(!frontier.empty()){
            float nextVal = -frontier.top().first;
            int nextInd = frontier.top().second;
            frontier.pop();
            if(!alreadyComputed[nextInd]){
                for(int i = 0;i < int(graphMesh[nextInd].size());i++){
                    int k = graphMesh[nextInd][i];
                    float newDist =  nextVal + shortDist(nextInd,k);
                    if(newDist<geoDistances[start][k][type]){
                        geoDistances[start][k][type] = newDist;
                        frontier.push(make_pair(-newDist,k));
                    }
                }
            }
            alreadyComputed[nextInd]=true;
        }
    }

    void Mesh::dijkstraFillEverything(){
        for(int i=0;i<n;i++){
            dijkstraFill(i);
        }
    }

    float MeshCorrespondance::simpleGlobalEval(){
        double res = 0;
        for(int i=0;i<mesh1.n;i++){
            res += mesh2.shortDist(i,corr[i])/mesh1.n;
        }
        return res;
    }

    float MeshCorrespondance::globalEval(){
        double res = 0;
        for(int i=0;i<mesh1.n;i++){
            res += geoDistances[i][corr[i]][1]/mesh1.n;
        }
        return res;
    }

    long long MeshCorrespondance::evaluate(const int& p,const set<int>& here){
        long long res = 0;
        for(int i=0;i<SAMPLE_SIZE;i++){
            int k = samplingVec[i];
            res+=100000*min((float)0.10,(float)abs(geoDistances[k][p][0] - geoDistances[corr[k]][corr[p]][1]))*(1.0-0.75*here.count(k));
        }
        long long res2 = 0;
        int cnt = 0;
        set<int>::iterator it;
        for(it=here.begin();it!=here.end();it++){
            int k = *it;
            float difference = abs(geoDistances[k][p][0] - geoDistances[corr[k]][corr[p]][1]);
            res2 -= 0.1*100000*(difference<0.01)*(1+5*isGood[k]);
            cnt++;
        }
        return res/SAMPLE_SIZE + 0.2*(cnt>0?res2/cnt:0);
    }

    long long MeshCorrespondance::evaluateLarge(const int& p,const set<int>& here){
        long long res = 0;
        for(int i=0;i<SAMPLE_SIZE;i++){
            int k = samplingVecLarge[i];
            res+=100000*min((float)0.10,(float)abs(geoDistances[k][p][0] - geoDistances[corr[k]][corr[p]][1]));
        }
        long long res2 = 0;
        int cnt = 0;
        set<int>::iterator it;
        for(it=here.begin();it!=here.end();it++){
            int k = *it;
            //float difference = abs(geoDistances[k][p][0] - geoDistances[corr[k]][corr[p]][1]);
            float difference = 0;
            res2 -= 0.1*100000*(difference<0.01)*(1+5*isGood[k]);
            cnt++;
        }
        return res/SAMPLE_SIZE + 0.1*(cnt>0?res2/cnt:0);
    }

    void MeshCorrespondance::optimize(const int& nSteps){
        for(int k=0;k<nSteps;k++){
            int i = rand()%(mesh1.n);

            for(int l=0;l<SAMPLE_SIZE;l++){
                samplingVec[l]=rand()%mesh1.n;
            }

            set<int> here = mesh1.bfs(i,BFS_NB);

            int theMin = evaluate(i,here);
            int jmin = corr[i];
            for(int jj=0;jj<300;jj++){
                int j = rand()%(mesh2.n);
                corr[i]=j;
                int aux = evaluate(i,here);
                if(aux<theMin){
                    theMin = aux;
                    jmin=j;
                }
            }
            corr[i]=jmin;

            set<int> there = mesh2.bfs(corr[i],BFS_NB);
            set<int>::iterator it;
            for(it=there.begin();it!=there.end();it++){
                int j = *it;
                corr[i]=j;
                int aux = evaluate(i,here);
                if(aux<theMin){
                    theMin = aux;
                    jmin=j;
                }
            }
            corr[i]=jmin;

            set<int> here2;

            isGood[i]=(evaluateLarge(i,here)<100000*0.01);
        }
    }

    void MeshCorrespondance::findGood(const int& nSteps){
        for(int k=0;k<nSteps;k++){
            int i = rand()%(mesh1.n);

            for(int l=0;l<SAMPLE_SIZE;l++){
                samplingVecLarge[l]=rand()%(mesh1.n);
            }

            set<int> emp;

            int theMin = evaluateLarge(i,emp);

            isGood[i]=(theMin<100000*0.01);
            isGood[i]=0;
        }
    }


    void MeshCorrespondance::showGroundTruthError(const string& name){
        string filepath = env.pathMesh + env.prefix + to_string(id1) + "_" + to_string(id2) + "_" + name + "_error.off";
        ofstream out (filepath);

        vector<bool> state(mesh1.n,0);

        vector<float> error(mesh1.n,0);

        map<int,int> newindex;
        for(int i=0;i<mesh1.n;i++){
            error[i]=geoDistances[i][corr[i]][1];
        }
        //sort(pointsVec.begin(),pointsVec.end());

        out<<"OFF\n";

        vector<float> colorTriangle(mesh1.p,0);
        for(int i=0;i<mesh1.p;i++){
            colorTriangle[i]=max(error[mesh1.triangles[i][0]],max(error[mesh1.triangles[i][1]],error[mesh1.triangles[i][2]]))/0.2;
            colorTriangle[i]=min((float)1.0,colorTriangle[i]);
        }

        out<<mesh1.n<<" "<<mesh1.p<<" "<<0<<"\n";
        for(int i=0;i<mesh1.n;i++){
            if(true){
                out<<mesh1.vert[i][0]<<" "<<mesh1.vert[i][1]<<" "<<mesh1.vert[i][2]<<"\n";
            }
        }
        for(int i=0;i<mesh1.p;i++){
            float r = colorTriangle[i];
            float g = 0.0;
            float b = 1.0 - colorTriangle[i];
            out<<3<<" "<<mesh1.triangles[i][0]<<" "<<mesh1.triangles[i][1]<<" "<<mesh1.triangles[i][2]<<" "<<strOfFloat(r)<<" "<<strOfFloat(g)<<" "<<strOfFloat(b)<<" "<<"1.0"<<"\n";
        }
    }

    void MeshCorrespondance::showGeodesicError(const string& name){
        string filepath = env.pathMesh + env.prefix + to_string(id1) + "_" + to_string(id2) + "_" + name + "_error2.off";
        ofstream out (filepath);

        vector<float> error(mesh1.n,0);

        set<int> emp;

        for(int i=0;i<mesh1.n;i++){
            error[i]=evaluate(i,emp);
        }

        out<<"OFF\n";

        vector<float> colorTriangle(mesh1.p,0);
        for(int i=0;i<mesh1.p;i++){
            colorTriangle[i]=max(error[mesh1.triangles[i][0]],max(error[mesh1.triangles[i][1]],error[mesh1.triangles[i][2]]))/(0.1*100000);
            colorTriangle[i]=min((float)1.0,colorTriangle[i]);
        }

        out<<mesh1.n<<" "<<mesh1.p<<" "<<0<<"\n";
        for(int i=0;i<mesh1.n;i++){
            if(true){
                out<<mesh1.vert[i][0]<<" "<<mesh1.vert[i][1]<<" "<<mesh1.vert[i][2]<<"\n";
            }
        }
        for(int i=0;i<mesh1.p;i++){
            float r = colorTriangle[i];
            float g = 0.0;
            float b = 1.0 - colorTriangle[i];
            out<<3<<" "<<mesh1.triangles[i][0]<<" "<<mesh1.triangles[i][1]<<" "<<mesh1.triangles[i][2]<<" "<<strOfFloat(r)<<" "<<strOfFloat(g)<<" "<<strOfFloat(b)<<" "<<"1.0"<<"\n";
        }
    }

    MeshCorrespondance::MeshCorrespondance(const string &pm, const string &pc, const string &pf, const int& k1, const int& k2){
        env.pathMesh = pm;
        env.pathCorr = pc;
        env.prefix = pf;

        mesh1.load(pm,pf + to_string(k1));
        mesh1.type = 0;
        id1 = k1;

        mesh2.load(pm,pf + to_string(k2));
        mesh2.type = 1;
        id2 = k2;

        if(true) loadCorr();

        setSampling();

        vector<bool> test(mesh1.n,0);
        isGood = test;
    }

    void MeshCorrespondance::dijkstra(){
        mesh1.dijkstraFillEverything();
        mesh2.dijkstraFillEverything();
    }

    void MeshCorrespondance::setSampling(){
        if(int(samplingVec.size())<=0){
            vector<int> listIndex;
            for(int i=0;i<mesh1.n;i++){
                listIndex.push_back(i);
            }

            for(int i=mesh1.n-1;i>=0;i--){
                int k = rand()%(i+1);
                int aux = listIndex[i];
                listIndex[i]=listIndex[k];
                listIndex[k]=aux;
            }

            for(int i=0;i<min(int(listIndex.size()),SAMPLE_SIZE);i++){
                samplingVec.push_back(listIndex[i]);
                indexS[listIndex[i]]=i;
            }
            sort(samplingVec.begin(),samplingVec.end());
        }

        if(int(samplingVecLarge.size())<=0){
            vector<int> listIndex;
            for(int i=0;i<mesh1.n;i++){
                listIndex.push_back(i);
            }

            for(int i=mesh1.n-1;i>=0;i--){
                int k = rand()%(i+1);
                int aux = listIndex[i];
                listIndex[i]=listIndex[k];
                listIndex[k]=aux;
            }

            for(int i=0;i<min(int(listIndex.size()),SAMPLE_SIZE_LARGE);i++){
                samplingVecLarge.push_back(listIndex[i]);
                indexSL[listIndex[i]]=i;
            }
            sort(samplingVecLarge.begin(),samplingVecLarge.end());
        }
    }

    void MeshCorrespondance::loadCorr(){
        string filepath = env.pathCorr + "corr" + to_string(id1) + "_" + to_string(id2) + ".txt";
        ifstream in (filepath);

        int sz;
        in>>sz;

        corr = vector<int>(0);

        for(int i=0;i<sz;i++){
            int j;in>>j;
            corr.push_back(j-1);
        }
    }

    void MeshCorrespondance::saveCorr(){
        string filepath = env.pathCorr + env.prefix + "corr" + to_string(id1) + "_" + to_string(id2) + "_" +  to_string(rand()%1000) + ".txt";
        ofstream out (filepath);

        int sz=corr.size();
        out<<sz<<"\n";

        for(int i=0;i<sz;i++){
            out<<corr[i]+1<<" ";
        }
    }

}
