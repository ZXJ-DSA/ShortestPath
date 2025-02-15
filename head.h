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
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "Heap.h"
//#include "labeling.hpp"

#define INF 99999999

typedef int vertex;
typedef unsigned long long ull;

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
vector<int> NodeOrder_;//nodeID order


struct DegComp{//min-first
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator< (const DegComp d) const{
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

struct OrderCompp{//prior to return the vertex with smaller order
    int x;
    OrderCompp(int _x){
        x=_x;
    }
    bool operator< (const OrderCompp& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};

vector<int> NodeOrders;
struct OrderComp{
    int x;
    int y;//order(x)<order(y)
    OrderComp(int _x, int _y){
        x=_x; y=_y;
    }
    bool operator< (const OrderComp& d) const{
        if(x==d.x && y==d.y){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrders[x]<NodeOrders[d.x];
            if(y!=d.y)
                return NodeOrders[y]<NodeOrders[d.y];
        }
    }
};

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const{
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

struct tri{
    int u;
    int v;//v is the contracted node
    int w;
};

struct Node{//tree node
    vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
//    vector<pair<int,pair<int,int>>> vertNo;//neighID/weight/count(how many ways can lead to this super edge weight)
    vector<int> pos;
    vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
//    vector<int> disNo;//the distance value of post-boundary strategy
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
    //vector<set<int>> FromNode;
//	set<int> changedPos;
    vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert), i.e., one-hop. there still be a case that the dis is also supported by neighbor's label even it is true.
    set<int> DisRe;//record the vertex id that the distance label should be updated
    vector<int> ch;
    int height=0;//tree height of a tree node
    int hdepth=0;//hdepty is the deepest node that a vertex still exists
    int pa;//parent, the pa of root vertex is 0
    int uniqueVertex;//?vertex id of this tree node?
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
    Node(){
        vert.clear();
//		neighInf.clear();
        pos.clear();
        dis.clear();
        cnt.clear();
        vAncestor.clear();
        ch.clear();
        pa = -1;
        uniqueVertex = -1;
        height = 0;
        hdepth = 0;
//		changedPos.clear();
        FN.clear();
        DisRe.clear();
//		piv.clear();
//        treeroot=-1;
    }
};

class Graph{
public:
    string dataset;
    string graph_path;
    int node_num = 0;                   //the number of vertices
    unsigned long long edge_num = 0;    //the number of edges
    vector<vector<pair<vertex,int>>> Neighbor;//adjacency list of original graph
    vector<Edge> Edges;//edge list for Bellman-Ford algorithm
    vector<unordered_map<vertex,int>> NeighborMap; //adjacency list of original graph, map version
    int threadnum = 10;// thread number
    int algo = 1;
    int orderStrategy = 2;// 0: read from disk; 1: MDE-based; 2:Degree-based
    bool ifDynamic = false;//if used for dynamic graphs

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
    void RetrievePath(int ID1, int ID2, vector<int> & prece);
    int DijkstraPaths(int ID1, int ID2);
    void RetrievePaths(int ID1, int ID2, vector<unordered_set<int>> & prece);
    int BiDijkstra(int ID1, int ID2);   //Bidirectional Dijkstra
    int Astar(int ID1, int ID2);    //A* algorithm
    int BiAstar(int ID1, int ID2);  //Bidirectional A*
    int BellmanFord(int ID1, int ID2);//Bellman-Ford's algorithm

    void DijkstraV(vector<pair<int,int>>& p);

    int EuclideanDis(int s, int t); //use original data
    int EuclideanDis2(int s, int t);//use processed data
    int EuclideanDis3(int s, int t);//use processed data

    /*-------------Index-based Algorithms--------------*/
    vector<int> NodeOrder;  //vertex order, (nodeID,order)
    vector<int> vNodeOrder; //inverted list, (order,nodeID)
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(threadnum);

    /// CH algorithm
    //CH with pruning
    vector<map<int,int>> OutEdgesM;
    vector<vector<pair<vertex, int>>> vvNode;//the adjacent and shortcut during contraction, duplicate of Neighbor
    vector<vector<pair<vertex, int>>> vvpShortCut;
    vector<vector<pair<vertex, int>>> AdjaShort;//final adja+shorcut for query answering
//    vector<vector<pair<int, int>>> AdjaShortR;//used for 1-N re-contraction, only store the lower adjacent nodes
//    vector<set<pair<int,int>>> InvalidWP;//record those invalid witness paths
//    unordered_map<pair<int,int>, int, hash_pair> EdgeRe;//edge (a,b)->edgeID
//    vector<unordered_map<pair<int,int>, int, hash_pair>> PathInfor; //u, (v,w) end points triple, length of path u<w
//    vector<vector<tri>> EdgeOnPath;//(a,b) end points of edge, path (u,v,w)
//    vector<map<int, vector<int>>> SupportNodes;//all possible supportive nodes of shortcut (u,w), uID<wID
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
    void CHConstruction(int strategy, bool ifDynamic);      //CH Construction with pruning
    int QueryCH(int ID1, int ID2);   //CH Query

    void IndexsizeCH();  //Core-tree index size computation

    //CH without pruning
    void CHConstructionWP(int strategy, bool ifDynamic);    //CH Construction without pruning
    int	QueryCHWP(int ID1, int ID2);

    void IndexsizeCHWP();  //Core-tree index size computation


    /// H2H algorithm
    vector<vector<int>> VidtoTNid;// (nodeID,vector<tree node rank>), record the tree node id whose unique vertex involves this vertex as neighbor
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<Node> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//For LCA query

    void H2HIndexConstruct(); //H2H index construction
    // MDE contraction
    void MDEContraction(string orderfile);
    //for order generation
    void deleteEOrderGenerate(int u,int v);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    void insertEMTorder(int u,int v,int w);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    /// Tree building
    void makeTree();
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    /// Index building
    void makeIndex();
    void makeRMQ();
    void makeRMQDFS(int p, int height);
    void makeIndexDFS(int p, vector<int>& list);

    void IndexsizeH2H();  //Core-tree index size computation

    //// For H2H query processing
    void EffiCheckH2H(string filename,int runtimes);
    int QueryH2H(int ID1,int ID2);
    int QueryH2HDebug(int ID1,int ID2);
    int LCAQuery(int _p, int _q);

    /// PLL algorithm
    hl::Labeling LabelV;//distance labels
    hl::PPR PPRV;//pruning point records
    vector<unordered_map<vertex,int>> Label;
    vector<unordered_map<vertex,unordered_set<vertex>>> PruningPointSet;//{(v,c),u}
    vector<unordered_map<vertex,vertex>> PruningPointSet2;//{(v,u),c}

    vector<omp_lock_t> oml;

    void PLLConstruction(int strategy, bool ifDynamic);
    void PLLConstruct();
    void PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor);
    void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);
    void DijksPrune1V(int nodeID, vector<vector<pair<vertex,int>>>& Neighbor);
    int PLLDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    int PLLDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode);

    void DegreeOrderGenerate(string orderfile);//degree-based vertex ordering
    void writePLL(string filename, string filenameP);
    void readPLL(string filename, string filenameP);
    void IndexSizePLL();

    int QueryPLL(int ID1, int ID2);

    /// PSL algorithm
    vector<unordered_map<vertex,int>> Dhop;
    vector<bool> DvertexNew;
    vector<vector<pair<vertex,int>>> DhopV;
    void PSLConstruction(int strategy, bool ifDynamic);
    void PSLConstruct();
    void PSLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    bool DhopLableRefreshStep(int step);
    bool DhopLableRefreshStepV(int step, vector<vector<pair<vertex,int>>> &Neighbor);
    void labelMultiThread2New(vector<unordered_map<vertex,int>>& newDhop, vector<int>& p, int step);
    void labelMultiThread2NewV(vector<vector<pair<vertex,int>>>& newDhop, int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void labelMultiThread2NewV2(vector<vector<pair<vertex,int>>>& newDhop, vector<int>& p, vector<vector<pair<vertex,int>>> &Neighbor);
    void threadDistribute(vector<vector<int>>& processID);
    int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    int ShortestDisQuery2(int ID1,int ID2);
    int ShortestDisQueryPeak(int ID1,int ID2,vector<int>& SupNode, int& d);
    int ShortestDisQueryPeakV(int ID1,int ID2, unordered_map<vertex,int>& Lh, vector<int>& SupNode, int& d);

    /// PCL algorithm
    int batchsize = threadnum;
    void PCLConstruction(int strategy, bool ifDynamic);
    void PCLConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    double PCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijks(vector<int> & IDs, vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijk(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PCLDijkV(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);

    double BPCLConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    double BPCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk(int nodeID, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijkV(int nodeID, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk2(vector<int>& p, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchPCLDijk2V(vector<int>& p, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor);

    void PruningPointBuild(bool ifParallel, vector<vector<int>> & processID, vector<vector<pair<vertex,int>>> &Neighbor);
    void ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID);
    void PPRConstruction(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstructionV(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstructionV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor);
    void PPRConstruction2(vector<int> & p, vector<vector<pair<vertex,int>>> &Neighbor);

    int ShortestDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode, int& d);

    /// GLL algorithm
    void GLLConstruction(int strategy, bool ifDynamic);
    double GLLConstruct(vector<vector<pair<vertex,int>>> &Neighbor);
    double GLLConstructV(vector<vector<pair<vertex,int>>> &Neighbor);
    void DijksPrune4(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void DijksPrune4V(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor);
    void BatchClean(vector<int> &bNodes, unordered_set<int>& setNodes, bool ifParallel);
    void LabelClean(int ID1, unordered_set<int>& setNodes);




    /*-------------Query Processing--------------*/
    int Query(int ID1, int ID2);
    void CorrectnessCheck(int num);
    void EfficiencyTest(int num);


    /*-------------Index Maintenance--------------*/
    void IndexMaintenance(int updateType, bool ifBatch, int batchNum, int batchSize);
    void IndexMaintenanceRealUpdate(string updateFile);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);

    /// For PLL/PSL algorithm
    set<pair<int,int>> NoSupportedPair;
    vector<unordered_set<int>> ChangedLabels;
    void IndexMaintenancePLL(int updateType, int batchNum);
    void DecreasePSL(vertex a, vertex b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors);
    void IncreasePSL(vertex a, vertex b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors);
    void CoarseUpdate(vertex LID, vertex HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<vertex,vertex>>& AL1, vector<pair<vertex,vertex>>& AL2, vector<pair<vertex,vertex>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor, bool ifDebug, int lid, int hid);//queue version
    void RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<vertex,vertex>>& AL1, vector<pair<vertex,vertex>>& AL2, vector<pair<vertex,vertex>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<vector<pair<vertex,int>>> &Neighbor, bool ifDebug, int lid, int hid);//queue version
    bool PPRCheck(vertex curID, vertex hubID, vector<vector<pair<vertex,int>>> &Neighbor, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<vertex,vertex>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid);//queue version
    void PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid);

    int DisQueryLower1(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors);
    int DisQueryPeak(vertex ID1, vertex ID2);
    int DisQueryVally(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors);
    pair<int,vertex> DisQueryVally2(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors);
    pair<int,vertex> DisQueryPeak2(vertex ID1, vertex ID2);

    /// For CHWP algorithm
    void IndexMaintenanceCHWP(int updateType, bool ifBatch, int batchNum, int batchSize);
    void IndexMaintenanceCHWPRealUpdate(string updateFile);
    void CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW increase
    void CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW decrease

    void CHIndexStore(string filename);
    void CHIndexCompare(string filename);

    //// For H2H Index Maintenance
    void IndexMaintenanceH2H(int updateType, bool ifBatch, int batchNum, int batchSize);
    void IndexMaintenanceH2HRealUpdate(string updateFile);
    void H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
    void H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel);

    void H2HIndexStore(string filename);
    void H2HIndexCompare(string filename);
};

#endif