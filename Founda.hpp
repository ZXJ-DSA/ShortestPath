//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef FOUNDA_HPP
#define FOUNDA_HPP

#include "head.h"


void Graph::ReadGraph(string filename){
    string r_edges = filename;
//    r_edges=filename+".new";
//    r_edges=filename+".new3";
//    r_edges=filename+".CH4";
//    r_edges=filename+".CH6";
    string r_co = filename+".co";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }
    cout<<"Graph file: "<<r_edges<<endl;
    Timer tt;
    tt.start();
    string line;
    int ID1,ID2,weight;
    unsigned long long edgeNum=0;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> node_num >> edge_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    cout<<"Node number: "<<node_num<<" ; Edge number: "<<edge_num<<endl;
    Neighbor.assign(node_num,vector<pair<vertex,int>>());
    vector<map<int,int>> NeighborMapT;
    NeighborMapT.assign(node_num,map<int,int>());
    bool flag_double=true;
    flag_double=false;
    if(edge_num>1000000000){
        flag_double= false;
    }

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
//            Neighbor[ID1].emplace_back(ID2,weight);
            if(NeighborMapT[ID1].find(ID2)==NeighborMapT[ID1].end()){//if not found
                NeighborMapT[ID1].insert({ID2,weight});
                NeighborMapT[ID2].insert({ID1,weight});
                edgeNum+=2;
            }else{
                if(NeighborMapT[ID1][ID2]!=weight){
                    cout<<"Inconsistent graph data!!! "<<ID1<<" "<<ID2<<" "<<NeighborMapT[ID1][ID2]<<" "<<weight<<endl; exit(1);
                }
            }

//            if(!flag_double){
////                Neighbor[ID2].emplace_back(ID1,weight);
//                NeighborMapT[ID2].insert({ID1,weight});
//                edgeNum+=2;
//            }else{
//                ++edgeNum;
//            }

            if(weight<=0){
                cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }
    inFile.close();
    if(edgeNum != edge_num){
        cout<<"Nominate edge number: "<<edge_num<<" ; real edge number: "<<edgeNum<<endl;
        edge_num=edgeNum;
        exit(1);
    }
    for(int i=0;i<NeighborMapT.size();++i){
        for(auto it=NeighborMapT[i].begin();it!=NeighborMapT[i].end();++it){
            Neighbor[i].emplace_back(it->first,it->second);
        }
    }
    tt.stop();
    int maxD = 0;
    for(auto it=Neighbor.begin();it!=Neighbor.end();++it){
        if(maxD < it->size()){
            maxD = it->size();
        }
    }
    cout<<"Maximal vertex degree: "<<maxD<<endl;
    cout<<"Graph read done. Time: "<<tt.GetRuntime()<<" s."<< endl;

//    ofstream outFile(graph_path+".2",ios::out);
//    if(!outFile.is_open()){
//        cout<<"Failed to open "<<graph_path+".2"<<endl;
//        exit(1);
//    }
//    for(int i=0;i<node_num;++i){
//        for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//            outFile<<i<<" "<<it->first<<endl;
//        }
//    }
//    outFile.close();

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
pair<int,unsigned long long> Graph::DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(nodenum,false);
    bool flag_finish = false;
    unsigned long long temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,unsigned long long> LCC;
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
            CHConstruction(orderStrategy,ifDynamic);
            break;
        }
        case 7:{//CHWP
            cout<<"CHWP algorithm!"<<endl;
            CHConstructionWP(orderStrategy,ifDynamic);
            break;
        }
        case 8:{//H2H
            cout<<"H2H algorithm!"<<endl;
            H2HIndexConstruct();
            break;
        }
        case 9:{//PLL
            cout<<"PLL algorithm!"<<endl;
            PLLConstruction(orderStrategy, ifDynamic);
            break;
        }
        case 10:{//PSL
            cout<<"PSL algorithm!"<<endl;
            PSLConstruction(orderStrategy, ifDynamic);
            break;
        }
        case 11:{//PCL
            cout<<"PCL algorithm!"<<endl;
            PCLConstruction(orderStrategy,ifDynamic);
            break;
        }
        case 12:{//GLL
            cout<<"GLL algorithm!"<<endl;
            GLLConstruction(orderStrategy,ifDynamic);
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
//            d=DijkstraPaths(ID1,ID2);// Dijkstra for multi shortest paths
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
            d=QueryCH(ID1,ID2);
            break;
        }
        case 7:{//CH
            d=QueryCHWP(ID1,ID2);
            break;
        }
        case 8:{
            d = QueryH2H(ID1,ID2);
            break;
        }
        case 9:{//PLL
            d=QueryPLL(ID1,ID2);//hash table-based implementation
//            d=LabelV.query(ID1,ID2);//vector-based implementation
            break;
        }
        case 10:{//PSL
            d=QueryPLL(ID1,ID2);//hash table-based implementation
//            d=LabelV.query(ID1,ID2);//vector-based implementation
            break;
        }
        case 11:{//PCL
            d=QueryPLL(ID1,ID2);//hash table-based implementation
//            d=LabelV.query(ID1,ID2);//vector-based implementation
            break;
        }
        case 12:{//GLL
            d=QueryPLL(ID1,ID2);//hash table-based implementation
//            d=LabelV.query(ID1,ID2);//vector-based implementation
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
    cout<<"Correctness check ("<<runtimes<<" rounds) ... "<<endl;
    Timer tt;
    double timeall=0;
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
        while(s==t){
            t=rand()%node_num;
        }
//        s=545, t=1000;//CH
//        s=1074, t=1839;//CHWP
//        s=115113,t=178597;//beijing
//        s=178597,t=88793;//beijing
//        s=243096,t=131349;//beijing
//        s=256545,t=161475;//beijing
//        s=140515,t=100026;//beijing
//        s=54815, t=193672;//CH

        d1=Dijkstra(s,t);
        tt.start();
        d2=Query(s,t);
        tt.stop();
        timeall+=tt.GetRuntime();

//        cout<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl;
//        QueryH2HDebug(s,t);
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<endl;
            QueryH2HDebug(s,t);
            exit(1);
        }
    }

//    cout<<"Average query time: "<<timeall*1000/runtimes<<" ms."<<endl;
}

void Graph::EfficiencyTest(int num) {
    srand (time(NULL));
    int s, t, d1, d2, d3;

    string r_query = graph_path+".query";
//    r_query="/Users/zhouxj/Documents/1-Research/Datasets/PopeElection2013/ave/PopeElection2013_OD_ShortDis_256.txt";

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
    int temp_i=0;

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
            temp_i++;
            if(temp_i>num){
                break;
            }
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<endl; exit(1);
        }
    }
    inFile.close();

    if(num>queryNum){
        num=queryNum;
    }
//    num=1;
//    ofstream outFile(r_query+"Result", ios::out);
    cout<<"Efficiency test ("<<num<<" rounds) ..."<<endl;
    vector<int> results(num,INF);
    Timer tt;
    tt.start();
//    DijkstraV(ODpairs);
    for(int i=0;i<num;i++){
//    for(int i=0;i<1;i++){
        s=ODpairs[i].first;
        t=ODpairs[i].second;
//        s=85250, t= 176775;
        d2=Query(s,t);
        results[i]=d2;
//        cout<<i<<": "<<s<<" "<<t<<" "<<d2<<endl;
//        outFile<<i<<": "<<s<<" "<<t<<" "<<d2<<endl;
//        cout<<d2<<" "<<tt.GetRuntime()*1000<<" ms"<<endl;
    }
    tt.stop();
    cout<<"Average query time: "<<tt.GetRuntime()*1000/num<<" ms."<<endl;
//    outFile.close();
}

void Graph::IndexMaintenance(int updateType, bool ifBatch, int batchNum, int batchSize){
    switch (algo) {
        case 7:{//CHWP
            IndexMaintenanceCHWP( updateType,ifBatch,batchNum,batchSize);
            break;
        }
        case 8:{//H2H
            IndexMaintenanceH2H(updateType,ifBatch,batchNum,batchSize);
            break;
        }
        case 9:{//PLL
            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 10:{//PSL
            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 11:{//PCL
            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 12:{//GLL
            IndexMaintenancePLL(updateType,batchNum);
            break;
        }

        default:
            break;
    }
}

void Graph::IndexMaintenanceRealUpdate(string updateFile){
    switch (algo) {
        case 7:{//CHWP
//            IndexMaintenanceCHWP( updateType,ifBatch,batchNum,batchSize);
            IndexMaintenanceCHWPRealUpdate(updateFile);
            break;
        }
        case 8:{//H2H
            IndexMaintenanceH2HRealUpdate(updateFile);
            break;
        }
        case 9:{//PLL
//            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 10:{//PSL
//            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 11:{//PCL
//            IndexMaintenancePLL(updateType,batchNum);
            break;
        }
        case 12:{//GLL
//            IndexMaintenancePLL(updateType,batchNum);
            break;
        }

        default:
            break;
    }
}

void Graph::VertexOrdering(int strategy) {
    string filename=graph_path+".order";
    filename=graph_path+".vertex_orderMDE2";
    switch (strategy) {
        case 0:{
            ReadOrder(filename);
            break;
        }
        case 1:{
            MDEOrderGenerate(filename);// for road networks
            break;
        }
        case 2:{
            DegreeOrderGenerate(filename);// for small-world networks
            break;
        }
        default:{
            cout<<"Wrong vertex ordering strategy! "<<strategy<<endl;
            exit(0);
        }
    }

    string orderfile = graph_path + ".order";
    int ifWrite=false;
//    ifWrite=true;
    if(ifWrite){
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << "\t" << NodeOrder[i] << endl;//ID, order
        ofile.close();
    }
}

void Graph::ReadOrder(string filename){
    string r_query = filename;

    ifstream inFile(r_query, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_query << endl;
        exit(1);
    }
    cout<<"Read vertex ordering... "<< filename<<endl;

    string line;
    int ID,order;
    int nodeNum=0;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> nodeNum)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }

    if(nodeNum != node_num){
        cout<<"Inconsistent vertex number! "<<nodeNum<<" "<<node_num<<endl;
        exit(0);
    }

    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> ID >> order)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }
        NodeOrder[ID] = order;
        vNodeOrder[order]= ID;
    }
    inFile.close();
    cout<<"Done."<<endl;
}

void Graph::ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldw;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>oldw;
        TestData.push_back(make_pair(make_pair(ID1, ID2), oldw));
    }
    IF.close();
}



#endif //FOUNDA_HPP