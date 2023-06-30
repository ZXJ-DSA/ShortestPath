//
// Created by Xinjie ZHOU on 24/05/2023.
//

#include "headPre.h"

void Graph::GraphDataGen(string filename, int strategy) {
    string r_edges = filename + ".original";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

    Timer tt;
    tt.start();
    string line;
    int ID1,ID2,weight,temp;
    unsigned long long edgeNum=0;
    int wRange=100;

    getline(inFile,line);
    istringstream iss(line);
    char ch;

    if (!(iss >> node_num >> edge_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
//    if (!(iss >>  ch >> node_num >> edge_num)){
//        cout<<"Wrong input syntax!"<<endl;
//        exit(1);
//    }
    cout<<"Node number: "<<node_num<<" ; Edge number: "<<edge_num<<endl;
    NeighborMap.assign(node_num,unordered_map<int,int>());

    unordered_set<int> set_A; set<int> set_LCC;

    if(strategy==0){//random strategy
        cout<<"Random strategy for edge weight. range: "<<wRange<<endl;
        while(getline(inFile,line)){
            if(line.empty())
                continue;
            istringstream iss(line);
        if (!(iss >> ID1 >> ID2 >> weight)){
//            if (!(iss >> ID1 >> ID2)){
//            if (!(iss >> ID1 >> ID2 >> weight >> temp)){
                cout<<"Wrong input syntax!"<<endl;
                exit(1);
            }

            ID1-=1; ID2-=1;
//            weight=rand()%wRange + 1;

            weight= max(weight,1);

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
    }
    else if(strategy==1){//degree-proportional strategy
        cout<<"Degree-proportional strategy for edge weight. range: "<<wRange<<endl;
        vector<int> degrees(node_num,0);

        while(getline(inFile,line)){
            if(line.empty())
                continue;
            istringstream iss(line);
//        if (!(iss >> ID1 >> ID2 >> weight)){
//            if (!(iss >> ID1 >> ID2)){
            if (!(iss >> ID1 >> ID2 >> weight >> temp)){
                cout<<"Wrong input syntax!"<<endl;
                exit(1);
            }

            ID1-=1; ID2-=1;
            weight = 0.1*(rand()%wRange) + 1;

            if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && ID1!=ID2){
                if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//if not found
                    NeighborMap[ID1].insert({ID2,weight});
                    NeighborMap[ID2].insert({ID1,weight});
                    edgeNum+=2;
                    set_A.insert(ID1); set_A.insert(ID2);
                    degrees[ID1]++; degrees[ID2]++;
                }

                if(weight<=0){
                    cout<<"Negative edge weight "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                }
            }else{
                cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
            }
        }

        int maxD = *max_element(degrees.begin(),degrees.end());
        double wei=0.0;
        cout<<"Maximum degree: "<<maxD<<endl;
        for(ID1=0;ID1<node_num;++ID1){
            for(auto it=NeighborMap[ID1].begin();it!=NeighborMap[ID1].end();++it){
                ID2 = it->first;
                wei = degrees[ID1]+degrees[ID2];
                wei = 1 - wei/(2*maxD);
                weight = wei * wRange + it->second;
                if(weight>0){
                    NeighborMap[ID1][ID2] = weight;
                }
            }
        }

    }


    inFile.close();
    tt.stop();
    cout<<"Graph read done. Time: "<<tt.GetRuntime()<<" s."<< endl;
    cout<<"Set_A: "<<set_A.size()<<"; edgeNum: "<<edgeNum<<endl;

    pair<int,unsigned long long> nums = DFS_CC(NeighborMap, set_A, set_LCC, node_num);
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
    cout<<"Writing graph into disk..."<<endl;
    Timer tt;
    tt.start();
    outFile << node_num <<" "<<edge_num<< endl;

    bool flag_double = true;
    if(edge_num>1000000000){
        flag_double=false;
        cout<<"Write edge once!"<<endl;
    }

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
                if(flag_double || IDMap[ID]<IDMap[it2->first]){
                    outFile<<IDMap[ID]<<" "<<IDMap[it2->first]<<" "<<it2->second<<endl;
                }

            }else{
                cout<<"Wrong! "<<ID<<" "<<it2->first<<" "<<it2->second<<endl;
            }

        }
    }

    outFile.close();
    tt.stop();
    cout<<"Write Done. Time: "<<tt.GetRuntime()<<" s."<< endl;

}

void Graph::ReadGraph(string filename){
    string r_edges = filename;
    string r_co = filename+".co";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

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

    bool flag_double=true;
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
            Neighbor[ID1].emplace_back(ID2,weight);
            if(!flag_double){
                Neighbor[ID2].emplace_back(ID1,weight);
                edgeNum+=2;
            }else{
                ++edgeNum;
            }

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
    tt.stop();
    cout<<"Graph read done. Time: "<<tt.GetRuntime()<<" s."<< endl;

    //read coordinates
//    ReadCoordinates(r_co);

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
pair<vertex, edge> Graph::DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
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

void Graph::Preprocessing() {

    switch (algo) {
        case 2:{
            cout<<"Query OD Pair Generation."<<endl;
            ODpairGenerate(10000);
            break;
        }
        case 3:{
            cout<<"Update OD Pair Generation!"<<endl;
            UpdateGenerate(1000);
            break;
        }
        case 4:{
            cout<<"Vertex Ordering Generation!"<<endl;
            VertexOrdering(strategy);//1: MDE-based ordering; 2: Degree-based ordering
            break;
        }


        default:
            break;
    }
}

//Function for generating random queries
void Graph::ODpairGenerate(int times){
    string RandomDis = graph_path + ".query";

    if(node_num == 0){

        ifstream inFile(graph_path, ios::in);
        if (!inFile) { // if not exist
            cout << "Fail to open file" << graph_path << endl;
            exit(1);
        }

        /// read graph and recording the degree of vertices
        inFile >> node_num >> edge_num;
        inFile.close();

    }

    /*---OD pairs generation---*/
    int pairs = 0;
    int node_start, node_end;
    int temp = 0;
    ofstream outFile(RandomDis, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Query OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    while (pairs < times) {
        node_start = rand() % node_num;
        node_end = rand() % node_num;
        while(node_end == node_start){
            node_end = rand() % node_num;
        }
//        outFile << node_start+1 << ' ' << node_end+1 << endl;
        outFile << node_start << ' ' << node_end << endl;
        ++pairs;
    }

    outFile.close();
    cout << "Finished." << endl;
}

//Function for generating update edges
void Graph::UpdateGenerate(int times){
    string filename = graph_path + ".update";
    /// read graph
    ifstream inFile(graph_path, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << graph_path << endl;
        exit(1);
    }
    inFile >> node_num >> edge_num;
    int ID1,ID2,weight;
    vector<pair<pair<int,int>,int>> edges;
    while(inFile){
        inFile >> ID1 >> ID2 >> weight;
        if(ID1 < ID2){
            edges.emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    inFile.close();

    int size = edges.size();
    assert(2*size == edge_num);
    /*---OD pairs generation---*/
    int pairs = 0;
    int id = 0;

    ofstream outFile(filename, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Update OD pairs file generating..." << endl;
    outFile << times << endl;
    //generate random OD pairs
    pairs = 0;
    unordered_set<int> edgeIdSet;
    edgeIdSet.clear();
    while (pairs < times) {
        id = rand() % size;
        if(edges[id].second >= 2 && edgeIdSet.find(id)==edgeIdSet.end()){//if edge weight is no smaller than 2, and it has not been added.
            edgeIdSet.insert(id);
            outFile << edges[id].first.first << ' ' << edges[id].first.second << ' '<<edges[id].second<< endl;
        }else{
            continue;
        }

        ++pairs;
    }
    outFile.close();
    cout << "Finished." << endl;
}
//function for vertex ordering
void Graph::VertexOrdering(int strategy) {
    ReadGraph(graph_path);
    ConnectivityCheck();
    switch (strategy) {
        case 1:{
            MDEOrderGenerate();// for road networks
            break;
        }
        case 2:{
            DegreeOrderGenerate();// for small-world networks
            break;
        }
        default:{
            cout<<"Wrong vertex ordering strategy! "<<strategy<<endl;
            exit(0);
        }
    }

    string orderfile = graph_path + ".order";
    int ifWrite=false;
    ifWrite=true;
    if(ifWrite){
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << "\t" << NodeOrder[i] << endl;//ID, order
        ofile.close();
    }
}

//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteECore(int u,int v){
//	if(Emap[u].find(v)!=Emap[u].end()){
//		Emap[u].erase(Emap[u].find(v));
//		DD[u]--;
//	}

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}
//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//		DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}
//function of generating CH shortcuts: MDE-based Elimination
void Graph::MDEOrderGenerate() {
    cout<<"MDE-based vertex ordering."<<endl;
    Timer tt;
    tt.start();
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }

    _DD_.assign(node_num,0); //_DD2_.assign(nodenum,0);
    DD.assign(node_num,0); //DD2.assign(nodenum,0);

    set<DegComp1> Deg;//min first
//    vector<bool> active(nodenum,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    unsigned int Twidth=0;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp1(i));
//            active[i] = true;
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }

    vNodeOrder.clear();
    vector<bool> existCore;
    existCore.assign(node_num,true);//if in the core, all vertices is originally in core
    vector<bool> change;
    change.assign(node_num,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);//temporal graph to store Neighbors in the core, for graph contraction

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first

        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp1(x));
            _DD_[x]=DD[x];
            Deg.insert(DegComp1(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }

        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(existCore[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        if(Twidth<Neigh.size()){
            Twidth = Neigh.size();
        }


        /// if still need to contract
        existCore[x]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(x,y);//delete x from y's adjacency list
            change[y]=true;
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();j++){
                ID2=Neigh[j].first;
                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);

                change[ID1]=true;
                change[ID2]=true;
            }
        }

    }


    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    tt.stop();


    cout<<"Time for node ordering: "<<tt.GetRuntime()<<" s. treewidth "<<Twidth<<endl;
}

void Graph::DegreeOrderGenerate(){
    cout<<"Degree-based vertex ordering."<<endl;
    Timer tt;
    tt.start();
    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);

    set<CompInt> Deg;
    int degree;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            Deg.insert(CompInt(i,degree));
        }else{
            cout<<"Degree 0! "<<i<<endl;
            exit(1);
        }
    }
    int id, rank=0;
    while(!Deg.empty()){
        id=Deg.begin()->id;
        NodeOrder[id]=rank;
        vNodeOrder[rank]=id;
        ++rank;
        Deg.erase(Deg.begin());
    }
    if(rank!=node_num){
        cout<<"Not all vertices are ranked! "<<rank<<" "<<node_num<<endl;
        exit(1);
    }
    tt.stop();

//    for(int i=node_num-1;i>node_num-5;--i){
//        cout<<"rank "<<i<<": "<<vNodeOrder[i]<<" "<<Neighbor[vNodeOrder[i]].size()<<endl;
//    }


    cout<<"Time for node ordering: "<<tt.GetRuntime()<<" s."<<endl;
}