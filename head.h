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
#include "Timer.h"
#include "Semaphore.h"

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

vector<int> _DD,_DD2;
struct DegComp{
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator< (const DegComp d) const{
        if(_DD[x]!=_DD[d.x])
            return _DD[x]<_DD[d.x];
        if(_DD2[x]!=_DD2[x])
            return _DD2[x]<_DD2[d.x];
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

    void ReadGraph(string filename);
    void ReadCoordinates(string filename);

    template <class T>
    pair<int,int> DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_B, int nodenum);
    void ConnectivityCheck();
    void Preprocess();
    void GraphDataGen(string filename);
    void WriteGraph(string filename, set<int>& set_LCC);

    int algo=1;

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
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    void deleteE(int u,int v);

    void insertEMTOrderGenerate(int u,int v,int w);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void MDEOrderGenerate(string orderfile);//MDE-based elimination for vertex ordering
    int writeShortCutorder(string filename);
    int ReadShortCut(string filename);
    int CHContractionSearchMT(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);
    void CHConstructionMT();
    void CHConstruction();      //CH Construction
    int CH(int ID1, int ID2);   //CH Query

    /// PLL algorithm
    vector<unordered_map<int,int>> Label;
    vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}

    void DegreeOrderGenerate(string orderfile);//degree-based vertex ordering
    void writePLL(string filename, string filenameP);
    void readPLL(string filename, string filenameP);
    int PLLDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
    void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);
    void PLLConstruct();
    void PLLConstruction();
    int PLL(int ID1, int ID2);



    /*-------------Query Processing--------------*/
    int Query(int ID1, int ID2);
    void CorrectnessCheck(int num);
    void EfficiencyTest(int num);




};

//created by Mengxuan, modified by Xinjie
namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

    template<int log_k, typename id_t, typename k_t >//id,value
    class heap {
    public:

        // Expose types.
        typedef k_t key_t;
        typedef id_t node_t;

        // Some constants regarding the elements.
        //static const node_t NULLINDEX = 0xFFFFFFFF;
        static const node_t k = 1 << log_k;//equals k = 1*2^log_k, usually log_k = 2, so k = 4

        // A struct defining a heap element.
        struct element_t {
            key_t key;
            node_t element;

            element_t() : key(0), element(0) {}

            element_t(const key_t k, const node_t e) : key(k), element(e) {}
        };


        //public:

        // Constructor of the heap.
        heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {//n is the number of elements in current
            // state, max_n is the size of heap
        }

        heap(): n(0), max_n(0), elements(0), position(0, NULLINDEX) {}

        ~heap(){}

        // Risize the heap
        inline void resize(node_t a){
            n = 0; max_n = a;
            elements.resize(a);
            position.resize(a, NULLINDEX);
        }

        // Size of the heap.
        inline node_t size() const {
            return n;
        }

        // Heap empty?
        inline bool empty() const {
            return size() == 0;
        }

        // Extract min element.
        inline void extract_min(node_t &element, key_t &key) {
            assert(!empty());

            element_t &front = elements[0];

            // Assign element and key.
            element = front.element;
            key = front.key;

            // Replace elements[0] by last element.
            position[element] = NULLINDEX;
            --n;//n=n-1
            if (!empty()) {
                front = elements[n];//elements[n] is the top element
                position[front.element] = 0;//make its position valid, it is also the smallest one
                sift_down(0);
            }
        }

        inline key_t top_key() {//get the key, i.e. minimal cost
            assert(!empty());

            element_t &front = elements[0];

            return front.key;

        }

        inline node_t top_value() {//get the value, i.e. id number of minimal cost

            assert(!empty());

            element_t &front = elements[0];

            return front.element;
        }

        // Update an element of the heap.
        inline void update(const node_t element, const key_t key) {

            if (position[element] == NULLINDEX) {//if originally NULL
                element_t &back = elements[n];//add new element to position n
                back.key = key;
                back.element = element;
                position[element] = n;//set position id to n
                sift_up(n++);
            } else {//if already valid, update the value
                node_t el_pos = position[element];//position information
                element_t &el = elements[el_pos];//get the element
                if (key > el.key) {//update the element
//                if (key > el.key || (key <= el.key && element > el.element)) {//update the element || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element)
                    el.key = key;
                    sift_down(el_pos);
                } else {
                    el.key = key;
                    sift_up(el_pos);
                }
            }
        }

        // Clear the heap.
        inline void clear() {
            for (node_t i = 0; i < n; ++i) {
                position[elements[i].element] = NULLINDEX;
            }
            n = 0;
        }

        // Cheaper erase.
        inline void erase(node_t v) {
            position[v] = NULLINDEX;
        }

        inline void clear_n() {
            n = 0;
        }

        // Test whether an element is contained in the heap.
        inline bool contains(const node_t element) const {
            return position[element] != NULLINDEX;
        }

        //return current elements information
        void get_elements(std::vector<std::pair<int,int>> &e_vector){
            std::pair<int,int> temp_pair;

            for(int i=0;i<n;i++){
                temp_pair.first = elements[i].key;
                temp_pair.second = elements[i].element;
                e_vector.push_back(temp_pair);
            }
        }

    protected:

        // Sift up an element.
        inline void sift_up(node_t i) {
            assert(i < n);
            node_t cur_i = i;
            while (cur_i > 0) {
                node_t parent_i = (cur_i - 1) >> log_k;//equals (cur_i - 1)/(2^log_k)
                if (elements[parent_i].key > elements[cur_i].key)//compare with parent node, if smaller, then swap
//                if (elements[parent_i].key > elements[cur_i].key || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element))//compare with parent node, if smaller, then swap
                    swap(cur_i, parent_i);
                else
                    break;
                cur_i = parent_i;
            }
        }

        // Sift down an element.
        inline void sift_down(node_t i) {
            assert(i < n);

            while (true) {
                node_t min_ind = i;
                key_t min_key = elements[i].key;

                node_t child_ind_l = (i << log_k) + 1;//equals i*2^log_k + 1
                node_t child_ind_u = std::min(child_ind_l + k, n);//equals min(child_ind_l+4,n)

                for (node_t j = child_ind_l; j < child_ind_u; ++j) {
                    if (elements[j].key < min_key) {
//                    if (elements[j].key < min_key || (elements[j].key >= min_key && elements[j].element < elements[i].element)) {
                        min_ind = j;
                        min_key = elements[j].key;
                    }
                }

                // Exchange?
                if (min_ind != i) {
                    swap(i, min_ind);
                    i = min_ind;
                } else {
                    break;
                }
            }
        }

        // Swap two elements in the heap.
        inline void swap(const node_t i, const node_t j) {
            element_t &el_i = elements[i];
            element_t &el_j = elements[j];

            // Exchange positions
            position[el_i.element] = j;
            position[el_j.element] = i;

            // Exchange elements
            element_t temp = el_i;
            el_i = el_j;
            el_j = temp;
        }

    private:

        // Number of elements in the heap.
        node_t n;

        // Number of maximal elements.
        node_t max_n;

        // Array of length heap_elements.
        std::vector<element_t> elements;

        // An array of positions for all elements.
        std::vector<node_t> position;
    };
}

#endif //HEAD_H
