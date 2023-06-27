//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef FOUNDA_HPP
#define FOUNDA_HPP

#include "head.h"

void Graph::GraphDataGen(string filename) {
    string r_edges = filename + ".original";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

    string line;
    int ID1,ID2,weight;
    int edgeNum=0;
    int wRange=10;

    getline(inFile,line);
    istringstream iss(line);
    char ch;
    if (!(iss >>  ch >> node_num >> edge_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    cout<<"Node number: "<<node_num<<" ; Edge number: "<<edge_num<<endl;
    NeighborMap.assign(node_num,unordered_map<int,int>());

    unordered_set<int> set_A; set<int> set_LCC;

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> ID1 >> ID2)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

        ID1-=1; ID2-=1;
        weight=rand()%wRange + 1;

        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
            if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//if not found
                NeighborMap[ID1].insert({ID2,weight});
                NeighborMap[ID2].insert({ID1,weight});
                edgeNum+=2;
                set_A.insert(ID1); set_A.insert(ID2);
            }

            if(weight<=0){
                cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
            }
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
        }
    }
    inFile.close();
    cout<<"Graph read done."<<endl;
    cout<<"Set_A: "<<set_A.size()<<"; edgeNum: "<<edgeNum<<endl;

    pair<int,int> nums = DFS_CC(NeighborMap, set_A, set_LCC, node_num);
    node_num=nums.first, edge_num=nums.second;

    WriteGraph(filename,set_LCC);

}

void Graph::WriteGraph(string filename, set<int>& set_LCC) {
    ofstream outFile(filename, ios::out);
    if (!outFile) {
        cout << "Write File opening failed. " << filename<<endl;
        assert(outFile);
        exit(1);
    }
    outFile << node_num <<" "<<edge_num<< endl;

    int ID=0;
    map<int,int> IDMap;//map from old ID to new ID

    for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
        IDMap.insert({*it,ID});
        ++ID;
    }
    assert(ID==node_num);

    for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
        ID=*it;
        for(auto it2=NeighborMap[ID].begin();it2!=NeighborMap[ID].end();++it2){
            if(IDMap.find(it2->first) != IDMap.end()){
                outFile<<IDMap[ID]<<" "<<IDMap[it2->first]<<" "<<it2->second<<endl;
            }else{
                cout<<"Wrong! "<<ID<<" "<<it2->first<<" "<<it2->second<<endl;
            }

        }
    }

    outFile.close();
    cout<<"Write Done."<<endl;

}

void Graph::ReadGraph(string filename){
    string r_edges = filename;
    string r_co = filename+".co";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

    string line;
    int ID1,ID2,weight;
    int edgeNum=0;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> node_num >> edge_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    cout<<"Node number: "<<node_num<<" ; Edge number: "<<edge_num<<endl;
    Neighbor.assign(node_num,vector<pair<int,int>>());

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> ID1 >> ID2 >> weight)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

//        ID1-=1; ID2-=1;

        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num){
            Neighbor[ID1].emplace_back(ID2,weight);
//            Neighbor[ID2].emplace_back(ID1,weight);
            ++edgeNum;
            if(weight<=0){
                cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
            }
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
        }
    }
    inFile.close();
    if(edgeNum != edge_num){
        cout<<"Nominate edge number: "<<edge_num<<" ; real edge number: "<<edgeNum<<endl;
        edge_num=edgeNum;
    }
    cout<<"Graph read done."<<endl;

    //read coordinates
    if(algo == 3 || algo == 4){
        ReadCoordinates(r_co);
    }

}

void Graph::ReadCoordinates(string filename) {

    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }

    string line;
    int ID,x,y;
    int nodeNum=0;
    bool flag= false;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> nodeNum)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    if(nodeNum != node_num){
        cout<<"Nominate node number of Coordinate file: "<<nodeNum<<" ; real node number: "<<node_num<<endl;
        flag=true;
    }
    Coordinate.assign(node_num,pair<int,int>());
    Coordinate2.assign(node_num,pair<int,int>());

    nodeNum=0;
    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> ID >> x >> y)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

        if(ID>=0 && ID<node_num){
            Coordinate[ID].first=x;
            Coordinate[ID].second=y;
            Coordinate2[ID].first=x*0.083907;
            Coordinate2[ID].second=y*0.111319;
            ++nodeNum;
        }else{
            cout<<"Graph data is wrong! "<<ID<<" "<<x<<" "<<y<<endl;
        }
    }
    inFile.close();
    if(nodeNum != node_num){
        cout<<"Nominate node number: "<<node_num<<" ; real node number of Coordinate file: "<<nodeNum<<endl;
        flag=true;
    }
    if(flag){
        exit(1);
    }
    cout<<"Coordinate read done."<<endl;
}


//function of checking the connectivity
template <class T>
pair<int,int> Graph::DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(nodenum,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());

    return make_pair(LCC.first.size(),LCC.second);
}


void Graph::ConnectivityCheck(){
    unordered_set<int> set_A; set<int> set_LCC;
    for(int i=0;i<node_num;++i){
        set_A.insert(i);
    }
    DFS_CC(Neighbor, set_A, set_LCC, node_num);
    cout<<"Connectivity check done."<<endl;
}

void Graph::Preprocess() {
    ConnectivityCheck();

    switch (algo) {
        case 1:{
            cout<<"Dijkstra's algorithm!"<<endl;
            break;
        }
        case 2:{
            cout<<"Bi-Dijkstra algorithm!"<<endl;
            break;
        }
        case 3:{
            cout<<"A* algorithm!"<<endl;
            break;
        }
        case 4:{
            cout<<"Bi-A* algorithm!"<<endl;
            break;
        }
        case 5:{
            cout<<"Bellman-Ford algorithm!"<<endl;
            for(int i=0;i<node_num;++i){
                for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
                    Edges.emplace_back(Edge(i,it->first,it->second));
                }
            }
            break;
        }
        case 6:{//CH
            cout<<"CH algorithm!"<<endl;
            CHConstruction();
            break;
        }
        case 7:{//PLL
            cout<<"PLL algorithm!"<<endl;
            PLLConstruction();
            break;
        }
        case 8:{//PSL
            cout<<"PSL algorithm!"<<endl;
            break;
        }

        default:
            break;
    }
}


int Graph::Query(int ID1, int ID2){
    int d=INF;
    switch (algo) {
        case 1:{
            d=Dijkstra(ID1,ID2);// Dijkstra
            break;
        }
        case 2:{
            d=BiDijkstra(ID1,ID2);// Bi-Dijkstra
            break;
        }
        case 3:{
            d=Astar(ID1,ID2);// A*
            break;
        }
        case 4:{
            d=BiAstar(ID1,ID2);//Bi-A*
            break;
        }
        case 5:{//Bellman-Ford
            d=BellmanFord(ID1,ID2);
            break;
        }
        case 6:{//CH
            d=CH(ID1,ID2);
            break;
        }
        case 7:{//PLL
            d=PLL(ID1,ID2);
            break;
        }
        case 8:{//PSL
            d=PLL(ID1,ID2);
            break;
        }
        default:
            break;
    }
    return d;
}

void Graph::CorrectnessCheck(int runtimes){
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ..."<<endl;
    Timer tt;
    double timeall=0;
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
        while(s==t){
            t=rand()%node_num;
        }
//        s=0, t=3;

        d1=Dijkstra(s,t);

        tt.start();
        d2=Query(s,t);
        tt.stop();
        timeall+=tt.GetRuntime();

//        cout<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl;

        }
    }

    cout<<"Average query time: "<<timeall*1000/runtimes<<" ms."<<endl;
}

void Graph::EfficiencyTest(int num) {
    srand (time(NULL));
    int s, t, d1, d2, d3;

    string r_query = graph_path+".query";

    ifstream inFile(r_query, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_query << endl;
        exit(1);
    }

    string line;
    int ID1,ID2;
    int queryNum=0;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> queryNum)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    vector<pair<int,int>> ODpairs;
//    ODpairs.assign(queryNum,pair<int,int>());

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> ID1 >> ID2)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

//        ID1-=1; ID2-=1;

        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num){
            ODpairs.emplace_back(ID1,ID2);
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<endl;
        }
    }
    inFile.close();

    if(num>queryNum){
        num=queryNum;
    }
//    num=1;

    cout<<"Efficiency test ("<<num<<" rounds) ..."<<endl;
    Timer tt;
    tt.start();
    for(int i=0;i<num;i++){
        s=ODpairs[i].first;
        t=ODpairs[i].second;
        d2=Query(s,t);
    }
    tt.stop();
    cout<<"Average query time: "<<tt.GetRuntime()*1000/num<<" ms."<<endl;
}

#endif //FOUNDA_HPP