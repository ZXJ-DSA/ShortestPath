//
// Created by Xinjie ZHOU on 24/05/2022.
//

#ifndef HEAD_H
#define HEAD_H

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

class Graph{
public:
    string dataset;
    string graph_path;
    int node_num = 0;                   //the number of vertices
    unsigned long long edge_num = 0;    //the number of edges
    vector<vector<pair<int,int>>> Neighbor;//adjacency list of original graph
    vector<Edge> Edges;//edge list for Bellman-Ford algorithm
    vector<unordered_map<int,int>> NeighborMap; //adjacency list of original graph, map version
    int threadnum = 10;// thread number
    int algo = 1;
    int orderStrategy = 2;// 0: read from disk; 1: MDE-based; 2:Degree-based

    /*-------------Preprocessing--------------*/
    void ReadGraph(string filename);
    void ReadCoordinates(string filename);
    template <class T>
    pair<int,unsigned long long> DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_B, int nodenum);
    void ConnectivityCheck();
    void Preprocess();
    void ReadOrder(string filename);
    void VertexOrdering(int strategy);//0: Read; 1: MDE; 2: Degree;


    /*-------------Index-free Algorithms--------------*/
    vector<pair<int,int>> Coordinate;//coordinates of vertex
    vector<pair<int,int>> Coordinate2;//coordinates of vertex, after conversion

    int Dijkstra(int ID1, int ID2); //Dijkstra's algorithm
    int BiDijkstra(int ID1, int ID2);   //Bidirectional Dijkstra
    int Astar(int ID1, int ID2);    //A* algorithm
    int BiAstar(int ID1, int ID2);  //Bidirectional A*
    int BellmanFord(int ID1, int ID2);//Bellman-Ford's algorithm

    int EuclideanDis(int s, int t); //use original data
    int EuclideanDis2(int s, int t);//use processed data

    /*-------------Index-based Algorithms--------------*/
    vector<int> NodeOrder;  //vertex order, (nodeID,order)
    vector<int> vNodeOrder; //inverted list, (order,nodeID)
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(threadnum);

    /// CH algorithm
    vector<vector<pair<int, int>>> vvNode;//the adjacent and shortcut during contraction, duplicate of Neighbor
    vector<vector<pair<int, int>>> vvpShortCut;
    vector<vector<pair<int, int>>> AdjaShort;//final adja+shorcut for query answering
//    vector<vector<pair<int, int>>> AdjaShortR;//used for 1-N re-contraction, only store the lower adjacent nodes
    vector<map<int, vector<int>>> SupportNodes;//all possible supportive nodes of shortcut (u,w), uID<wID
    vector<int> DD,DD2;
    vector<map<int,pair<int,int>>> E;
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<map<int, vector<int>>> SCconNodesMT;
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void MDEOrderGenerate(string orderfile);//MDE-based elimination for vertex ordering
    int writeShortCutorder(string filename);
    int ReadShortCut(string filename);
    int CHContractionSearchMT(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);
    void CHConstructionMT();
    void CHConstruction();      //CH Construction
    int CH(int ID1, int ID2);   //CH Query

    /// PLL algorithm
    vector<unordered_map<int,int>> Label;
//    vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}

    void PLLConstruction(int strategy);
    void PLLConstruct();
    void DegreeOrderGenerate(string orderfile);//degree-based vertex ordering
    void writePLL(string filename, string filenameP);
    void readPLL(string filename, string filenameP);
    int PLLDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);

    int PLL(int ID1, int ID2);

    /// PSL algorithm
    vector<unordered_map<int,int>> Dhop;
    vector<int> Dvectex;
    vector<bool> DvertexNew;

    void PSLConstruction(int strategy);
    void PSLConstruct();
    bool DhopLableRefreshStep(int step);
    void labelMultiThread2New(vector<unordered_map<int,int>>& newDhop, vector<int>& p,int step);
    void threadDistribute(vector<vector<int>>& processID);
    int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    int ShortestDisQuery2(int ID1,int ID2);

    /*-------------Query Processing--------------*/
    int Query(int ID1, int ID2);
    void CorrectnessCheck(int num);
    void EfficiencyTest(int num);




};



#endif //HEAD_H
