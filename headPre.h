//
// Created by Xinjie ZHOU on 24/05/2022.
//
#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <assert.h>
#include <boost/thread/thread.hpp>
#include "Heap.h"

#define INF 99999999

using namespace std;



struct CompDouble{//maximal first
    int id;
    double value;
    CompDouble(int _id, double _cc){
        id = _id; value = _cc;
    }
    bool operator < (const CompDouble b) const{
        return (value > b.value) || (b.value <= value && (id < b.id));
    }
};

struct CompInt{//minimum value and id first
    int id;
    int value;
    CompInt(int _id, int _cc){
        id = _id; value = _cc;
    }
    bool operator < (const CompInt b) const{
        return (value < b.value) || (b.value >= value && (id < b.id));
    }
};

vector<int> _DD_;
struct DegComp1{//min-first
    int x;
    DegComp1(int _x){
        x=_x;
    }
    bool operator< (const DegComp1 d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
};

struct Edge {
    int s, t;
    int weight;
    Edge():s(0),t(0),weight(-1){};
    Edge(int s_,int t_,int w){
        s=s_,t=t_,weight=w;
    };
};

typedef int vertex;
typedef unsigned long long edge;

class Graph{
public:
    string dataset;
    string graph_path;
    vertex node_num = 0;    //the number of vertices
    edge edge_num = 0;    //the number of edges
    vector<vector<pair<vertex,int>>> Neighbor;//adjacency list of original graph
    vector<Edge> Edges;//edge list for Bellman-Ford algorithm
    vector<unordered_map<int,int>> NeighborMap; //adjacency list of original graph, map version
    int threadnum = 10;// thread number
    int task=1;//task
    int strategy=2;//1: MDE-based ordering; 2: Degree-based ordering
//    map<long long int,int> OldToNewID;//map from old vertex id to new vertex id

    vector<pair<int,int>> Coordinate;//coordinates of vertex
    vector<pair<int,int>> Coordinate2;//coordinates of vertex, after conversion
    double graphDiameter=0;//in kilometer

    /*-------------Preprocessing--------------*/
    void ReadGraph(string filename, bool ifCo);
    void ReadCoordinates(string filename);
    void ReadRoadNetworksPre();
    void ReadRoadNetworksPre2();

    template <class T>
    pair<vertex,edge> DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_B, vertex nodenum);
    void ConnectivityCheck();
    void Preprocessing();
    void GraphDataGen(string filename,int strategy);//strategy: 0: random; 1: degree-proportional
    void WriteGraph(string filename, set<int>& set_LCC);

    void ODpairGenerate(int times);    //used to generate random OD pairs
    void ODpairGenerateDistanceDijk(int times);//use the shortest path distance
    void ODpairGenerateDistanceSpatial(int times);//use the spatial distance
    void UpdateGenerate(int times); //used to generate update pairs
    void UpdateGenerate2(int times); //used to generate increase and decrease update pairs

    vector<int> NodeOrder;  //vertex order, (nodeID,order)
    vector<int> vNodeOrder; //inverted list, (order,nodeID)
    vector<int> DD;
    vector<map<int,pair<int,int>>> E;
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    void VertexOrdering(int strategy);//1: MDE; 2: Degree
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void MDEOrderGenerate();//MDE-based elimination for vertex ordering
    void DegreeOrderGenerate();

    void RoadNetworkPreprocessing();
    double EuclideanDis(pair<double,double> s, pair<double,double> t);//input longitude and latitude, output distance in kilometers
    double ComputeDiameter(vector<pair<double,double>> & points);
    void ComputeDiameterOfNetwork(string filename);
    void ComputeMaximalStepNumOfShortestPath();
    void ComputeLongestShortestPath(unordered_set<int> & sourceIDs);
    pair<unsigned long int, unsigned long int> DijkstraSSSP(int ID); //Dijkstra's algorithm for SSSP, return the longest shortest path length and maximal step number of SP

    void ExtractSubgraphForScalabilityTest(bool ifCo);
};

