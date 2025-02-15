// Basic classes, Timer, Semaphore, Heap
// created by Mengxuan, modified by Xinjie

#pragma once

#include <mutex>
#include <condition_variable>
#include <vector>
#include <chrono>

struct Timer
{
    std::chrono::high_resolution_clock::time_point t1, t2;//varibles for time record
    std::chrono::duration<double> time_span;
    void start()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }
    void stop()
    {
        t2 = std::chrono::high_resolution_clock::now();
    }
    double GetRuntime()//return time in second
    {
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);//std::chrono::nanoseconds
        return time_span.count();
    }
};

class Semaphore
{
public:
    Semaphore (int count_ = 0)
            : count(count_) {}

    inline void notify()
    {
        std::unique_lock<std::mutex> lock(mtx);
        count++;
        cv.notify_one();
    }

    inline void wait()
    {
        std::unique_lock<std::mutex> lock(mtx);

        while(count == 0)
        {
            cv.wait(lock);

        }
        count--;
    }

private:
    std::mutex mtx;
    std::condition_variable cv;
    int count;//the number of thread waiting here?
};

namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

    template<int log_k, typename id_t, typename k_t >//id,value
    class heap {//binary heap
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

// Hub Labels are lists of hubs and distances to them attached to every vertex in a graph.
// This file contains the class to store labels.
// Available methods include making a query, write labels to file, read labels from file,
// merge label tables, clean label tables.
//
//  Author: Xinjie ZHOU

#include <omp.h>
namespace hl {

#define INF 99999999
//typedef unsigned int vertex;
    typedef int vertex;
//typedef unsigned long long ull;

// Class to store labels
    class Labeling {

    public:
        std::vector< std::vector< std::pair<vertex,int> > > Labels;     // Lists of forward/reverse hubs

        vertex n;
        bool *lck;

        Labeling(size_t n = 0) :
                Labels(n, std::vector< std::pair<vertex,int> >()),
                n(n),
                lck(new bool[n]()) {}

        void resize(size_t n_=0){
            Labels.assign(n_, std::vector< std::pair<vertex,int> >()),
                    n=n_,
                    lck=new bool[n]();
        }

        void clear(){
            clearLabel();
//        delete[] lck;
        }

        // Find u-v distance
        int query(vertex u, vertex v) {
            int r = INF;
            for (size_t i=0, j=0; i < Labels[u].size() && j < Labels[v].size();) {
                if (Labels[u][i].first == Labels[v][j].first) {
                    r = std::min(r, Labels[u][i++].second + Labels[v][j++].second);
                } else if (Labels[u][i].first < Labels[v][j].first) ++i;
                else ++j;
            }
            return r;
        }

        bool cover(vertex u, vertex v, uint d=INF) {
            for (size_t i=0, j=0; i < Labels[u].size() && j < Labels[v].size();) {
                if (Labels[u][i].first == Labels[v][j].first) {
                    if (d >= Labels[u][i++].second + Labels[v][j++].second) return true;
                } else if (Labels[u][i].first < Labels[v][j].first) ++i;
                else ++j;
            }
            return false;
        }


//    inline bool clean_cover(Vertex u, Vertex v, unsigned int f, Distance d=infty, size_t hub_order=0) {
//        for (size_t i=0, j=0; i < label_v[u][f].size() && j < label_v[v][!f].size();) {
//            if (label_v[u][f][i] >= hub_order || label_v[v][!f][j] >= hub_order) return false;
//            if (label_v[u][f][i] == label_v[v][!f][j]) {
//                if (d >= label_d[u][f][i++] + label_d[v][!f][j++]) return true;
//            }
//            else if (label_v[u][f][i] < label_v[v][!f][j]) ++i;
//            else ++j;
//        }
//        return false;
//    }
//
//    inline void clean_roots(Vertex v, std::vector<Vertex> &order, unsigned int side)
//    {
//        std::vector<Vertex> temp_v;
//        std::vector<Distance> temp_d;
//        for (size_t i=0; i<label_v[v][side].size(); i++)
//        {
//            size_t hub_order = label_v[v][side][i];
//            Vertex hub =  order[hub_order];
//            Distance hub_dist = label_d[v][side][i];
//            if (!clean_cover(hub, v, side, hub_dist, hub_order))//;
//            {
//                temp_v.push_back(hub_order);
//                temp_d.push_back(hub_dist);
//            }
//        }
//        temp_v.swap(label_v[v][side]);
//        temp_d.swap(label_d[v][side]);
//    }

        // Add hub (v,d) to forward or reverse label of u
        inline void add(vertex u, vertex v, uint d) {

            while (!__sync_bool_compare_and_swap(&lck[u], false, true)) {}
            Labels[u].emplace_back(v,d);
            lck[u]=false;
        }
        void add_lockfree(vertex u, vertex v, uint d) {
            Labels[u].emplace_back(v,d);
        }

        // Get labels
//    std::vector< std::vector<Vertex> > &get_label_hubs(Vertex u) { return label_v[u]; }
//    std::vector< std::vector<Distance> > &get_label_distances(Vertex u) { return label_d[u]; }

        // Get maximum label size
        size_t get_max() const {
            size_t max = 0;
            for (vertex v = 0; v < n; ++v){
                max = std::max(max, Labels[v].size());
            }

            return max;

            //size_t maxVal = 0;
            ////#pragma omp parallel for num_threads (NUM_THREAD) reduction (max: maxVal)
            ////{
            //    for (Vertex v=0; v<n; v++)
            //    {
            //        for (unsigned side = 0; side < 2; side++)
            //        {
            //            if (label_v[v][side].size() < maxVal)
            //                maxVal = label_v[v][side].size;
            //        }
            //    }
            //}
            //return maxVal;



        }

        // Get average label size
        double get_avg() const {
            long long total = 0;
            for (vertex v = 0; v < n; ++v)
                total += Labels[v].size() + Labels[v].size();
            return static_cast<double>(total)/n/2;
        }
        long long get_total() const {
            long long total = 0;
            for (vertex v = 0; v < n; ++v)
                total += Labels[v].size() + Labels[v].size();
            return total;
        }

        // Write labels to file
        bool write(char *filename) {
            std::ofstream file;
            file.open(filename);
            file << n << std::endl;
            for (vertex v = 0; v < n; ++v) {
                file << Labels[v].size();
                for (size_t i = 0; i < Labels[v].size(); ++i) {
                    file << " " << Labels[v][i].first;
                    file << " " << Labels[v][i].second;
                }
                file << std::endl;
            }
            file.close();
            return file.good();
        }

        // Read labels from file
        bool read(char *filename, vertex check_n = 0) {
            std::ifstream file;
            file.open(filename);
            file >> n;
            if (check_n && n != check_n) return false;
            Labels.resize(n, std::vector< std::pair<vertex,int> >());

            for (vertex v = 0; v < n; ++v) {
                size_t s;
                file >> s;
                Labels[v].resize(s);
                Labels[v].resize(s);
                for (size_t i = 0; i < s; ++i) {
                    file >> Labels[v][i].first;
                    file >> Labels[v][i].second;
                }
            }
            file >> std::ws;
            file.close();
            return file.eof() && !file.fail();
        }

        // Clear labels
        void clearLabel(unsigned NUM_THREAD=72) {
#pragma omp parallel for num_threads(NUM_THREAD)
            for (vertex v = 0; v < n; ++v) {
                Labels[v].clear();
            }
        }

        // Sort labels before making queries
        void sort(unsigned NUM_THREAD) {

            // std::vector<std::vector<std::pair<Vertex, Distance>>> label(NUM_THREAD)
            // maxSize = get_max();
            // for (int i=0; i<NUM_THREAD; i++)
            //     label[i] = std::vector<std::pair<Vertex, Distance>>.reserve(maxSize)
            //#pragma omp parallel for num_threads(NUM_THREAD) schedule (dynamic)
            //for (Vertex v = 0; v < n; ++v) {
            //    for (int side = 0; side < 2; ++side) {
            //        for (size_t i = 0; i < label_v[v][side].size(); ++i)
            //            label[i].push_back(std::make_pair(label_v[v][side][i], label_d[v][side][i]));
            //        std::sort(label.begin(),label.end());
            //        for (size_t i = 0; i < label_v[v][side].size(); ++i) {
            //            label_v[v][side][i] = label[i].first;
            //            label_d[v][side][i] = label[i].second;
            //        }
            //        label[i].clear();
            //    }
            //}
            Timer tt;
            tt.start();
#pragma omp parallel for num_threads(NUM_THREAD) schedule (dynamic, NUM_THREAD)
            for (vertex v = 0; v < n; ++v) {
                std::sort(Labels[v].begin(),Labels[v].end());
            }
            tt.stop();
            std::cout<<"Time for label sorting: "<<tt.GetRuntime()<<" s."<<std::endl;
        }

        void postProcess(std::vector<std::unordered_map<vertex,int>>& Label, int ThreadNum){
            Timer tt;
            tt.start();
            Label.clear();
            Label.assign(n,std::unordered_map<vertex,int>());

#pragma omp parallel for num_threads(ThreadNum) schedule (dynamic, ThreadNum)
            for(int id=0;id<n;++id){
                vertex id2;
                uint dis;
                for(auto it=Labels[id].begin();it!=Labels[id].end();++it){
                    id2=it->first; dis=it->second;
                    Label[id].insert({id2,dis});
                }
//                Labels[id].clear();
            }
//            Labels.clear();
            tt.stop();
            std::cout<<"Time for label post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
        }



    };

    struct Node{
        int ID;
        omp_lock_t lock;
    };

// Class to store pruning point records
    class PPR{

    public:
        std::vector< std::vector< std::pair<vertex,vertex> > > PPRs;     // Lists of pruning point records
//    std::vector< std::unordered_map< vertex, std::unordered_set<vertex> > > PPRSet; // Hash table-based
        std::vector<Semaphore*> vSm;//locks of boost

        vertex n;
        bool *lck;

        PPR(size_t n = 0) : n(n), lck(new bool[n]()) {
            PPRs.assign(n, std::vector< std::pair<vertex,vertex> >());
        }

        void resize(size_t n_=0){
            PPRs.assign(n_, std::vector< std::pair<vertex,vertex> >()),
                    n=n_,
                    lck=new bool[n]();
        }

        void clear(){
            PPRs.clear();
//        delete[] lck;
            vSm.clear();
        }

        // Add hub (v,d) to forward or reverse label of u
        inline void add(vertex u, vertex v, uint d) {

            while (!__sync_bool_compare_and_swap(&lck[u], false, true)) {}
            PPRs[u].emplace_back(v,d);
            lck[u]=false;
        }
        void add_lockfree(vertex u, vertex v, uint d) {
            PPRs[u].emplace_back(v,d);
        }

        void write(std::string filename){
            std::ofstream OF2(filename);
            if(!OF2){
                std::cout<<"Cannot open "<<filename<<std::endl;
                exit(1);
            }
            std::cout<<"Write PPR..."<<std::endl;
//        Timer tt;
//        tt.start();
            int hub,ID2;
            for(int ID1=0;ID1<n;ID1++){
                if(!PPRs[ID1].empty()){//Order
                    OF2<<ID1;
                    for(auto it=PPRs[ID1].begin();it!=PPRs[ID1].end();++it){//Order
                        hub=it->first; ID2=it->second;
                        OF2<<" "<<hub<<" "<<ID2;
                    }
                    OF2<<std::endl;
                }

            }
            OF2.close();
//        tt.stop();
//        std::cout<<"Done."<<std::endl;
        }

        void pprInsert(std::vector<vertex>& p,std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2){
            for(auto it=p.begin();it!=p.end();++it){
                vertex id=*it;
                vertex hub,id2;
                for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
                    hub=it->first; id2=it->second;

                    vSm[id]->wait();
                    PruningPointSet[id][hub].insert(id2);
                    PruningPointSet2[id][id2]=hub;
                    vSm[id]->notify();

                    vSm[id2]->wait();
                    PruningPointSet[id2][hub].insert(id);
                    PruningPointSet2[id2][id]=hub;
                    vSm[id2]->notify();
                }
                PPRs[id].clear();
            }

        }

        // post-processing
        void postProcess(std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2, std::vector<vertex>& vertices){
            Timer tt;
            tt.start();
            PruningPointSet.assign(n,std::unordered_map<vertex,std::unordered_set<vertex>>());
            PruningPointSet2.assign(n,std::unordered_map<vertex,vertex>());

            vertex id,hub,id2;
            for(int i=0;i<vertices.size();++i){
                id=vertices[i];

                for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
                    hub=it->first; id2=it->second;

                    PruningPointSet[id][hub].insert(id2);
                    PruningPointSet2[id][id2]=hub;

                    PruningPointSet[id2][hub].insert(id);
                    PruningPointSet2[id2][id]=hub;
                }
                PPRs[id].clear();
            }

            PPRs.clear();

            tt.stop();
            std::cout<<"Time for PPR post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
        }

        void postProcess(std::vector<std::unordered_map<vertex,std::unordered_set<vertex>>>& PruningPointSet, std::vector<std::unordered_map<vertex,vertex>>& PruningPointSet2, std::vector<vertex>& vertices, int ThreadNum){
            Timer tt;
            tt.start();
            PruningPointSet.assign(n,std::unordered_map<vertex,std::unordered_set<vertex>>());
            PruningPointSet2.assign(n,std::unordered_map<vertex,vertex>());

            /// boost-based implementation
            vSm.reserve(n);
            for(int i = 0; i < n; i++)
            {
                Semaphore* s = new Semaphore(1);
                vSm.push_back(s);
            }
            std::vector<std::vector<vertex>> processID(ThreadNum);
            int pid;
            for(int i=0;i<vertices.size();++i){
                pid=i%ThreadNum;
                processID[pid].push_back(vertices[i]);
            }
            boost::thread_group thread;
            for(int i=0;i<processID.size();i++){
                thread.add_thread(new boost::thread(&PPR::pprInsert, this, boost::ref(processID[i]), boost::ref(PruningPointSet), boost::ref(PruningPointSet2)));
            }
            thread.join_all();

//        /// multiple threads: openmp-based implementation
//        std::vector<Node> vs(n);
//        for(int i=0;i<n;++i){
//            vs[i].ID=i;
//            omp_init_lock(&vs[i].lock);
//        }
//        #pragma omp parallel for num_threads(ThreadNum) schedule (dynamic, ThreadNum)
//        for(int id=0;id<n;++id){
//            vertex hub,id2;
//            for(auto it=PPRs[id].begin();it!=PPRs[id].end();++it){
//                hub=it->first; id2=it->second;
//
//                omp_set_lock(&vs[id].lock);
//                PruningPointSet[id][hub].insert(id2);
//                PruningPointSet2[id][id2]=hub;
//                omp_unset_lock(&vs[id].lock);
//
//                omp_set_lock(&vs[id2].lock);
//                PruningPointSet[id2][hub].insert(id);
//                PruningPointSet2[id2][id]=hub;
//                omp_unset_lock(&vs[id2].lock);
//            }
//            PPRs[id].clear();
//        }
//        for (int i = 0; i < n; ++i) {
//            omp_destroy_lock(&vs[i].lock);// Destroy locks
//        }



            PPRs.clear();

            tt.stop();
            std::cout<<"Time for PPR post-processing: "<<tt.GetRuntime()<<" s."<<std::endl;
        }


    };

}



