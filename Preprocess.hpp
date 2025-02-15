//
// Created by Xinjie ZHOU on 24/05/2023.
//

#include <boost/algorithm/string.hpp>
#include "headPre.h"

void Graph::GraphDataGen(string filename, int strategy) {
    string r_edges = filename + ".original";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

    cout<<"Graph edge weight generation."<<endl;

    Timer tt;
    tt.start();
    string line;
    int ID1,ID2,weight,temp;
    unsigned long long edgeNum=0;
    int wRange=100;
    weight=1;

    getline(inFile,line);
    vector<string> re1;

    if(node_num==0){
        boost::split(re1,line,boost::is_any_of(" \t"));
        if(re1.size()==2){
            node_num=stoi(re1[0]), edge_num=stoul(re1[1]);
        }else if(re1.size()>=3 && re1[0]=="#"){
            node_num=stoi(re1[1]), edge_num=stoul(re1[2]);
        }else{
            cout<<"Wrong input! re.size: "<<re1.size()<<". " <<line<<endl;
            exit(1);
        }
    }


    cout<<"Node number: "<<node_num<<endl;
    NeighborMap.assign(node_num,unordered_map<int,int>());

    unordered_set<int> set_A; set<int> set_LCC;

    switch (strategy) {
        case 1:{
            cout<<"Random strategy for edge weight. range: "<<wRange<<endl;
            while(getline(inFile,line)){
                if(line.empty())
                    continue;
                vector<string> re1;
                boost::split(re1,line,boost::is_any_of(" \t"));
                if(re1.size()>=2){
                    ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                }
                else{
                    cout<<"Wrong line. "<< line<<endl;
                    exit(1);
                }
                istringstream iss(line);

                ID1-=1; ID2-=1;
                weight=rand()%wRange + 1;
//                weight= max(weight,1);

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
                    if(ID1!=ID2){
                        cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                        exit(1);
                    }
                }
            }
            break;
        }
        case 2:{
            cout<<"Degree-proportional strategy for edge weight. range: "<<wRange<<endl;
            vector<int> degrees(node_num,0);

            while(getline(inFile,line)){
                if(line.empty())
                    continue;
                vector<string> re1;
                boost::split(re1,line,boost::is_any_of(" \t"));
                if(re1.size()>=2){
                    ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                }
                else{
                    cout<<"Wrong line. "<< line<<endl;
                    exit(1);
                }

                ID1-=1; ID2-=1;
//                weight = 0.05*(rand()%wRange) + 1;

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
                    if(ID1!=ID2){
                        cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                        exit(1);
                    }
                }
            }

            int maxD = *max_element(degrees.begin(),degrees.end());
            double wei=0.0;
            cout<<"Maximum degree: "<<maxD<<endl;
            for(ID1=0;ID1<node_num;++ID1){
                for(auto it=NeighborMap[ID1].begin();it!=NeighborMap[ID1].end();++it){
                    ID2 = it->first;
                    wei = degrees[ID1]+degrees[ID2];
                    weight = 2*maxD/wei;
                    if(weight>0){
                        NeighborMap[ID1][ID2] = weight;
                    }else{
                        cout<<"Wrong edge weight. "<<weight<<endl;
                    }
                }
            }
            break;
        }
        case 3:{
            cout<<"Degree as numerator strategy for edge weight. range: "<<wRange<<endl;
            vector<int> degrees(node_num,0);

            while(getline(inFile,line)){
                if(line.empty())
                    continue;

                re1.clear();
                boost::split(re1,line,boost::is_any_of(" \t"));
                if(re1.size()>=2){
                    ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                }
                else{
                    cout<<"Wrong line. "<< line<<endl;
                    exit(1);
                }

                ID1-=1; ID2-=1;
                weight = 0.05*(rand()%wRange) + 1;

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
                    if(ID1!=ID2){
                        cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                        exit(1);
                    }
                }
            }

            int maxD = *max_element(degrees.begin(),degrees.end());
            double wei=0.0;
            cout<<"Maximum degree: "<<maxD<<endl;
            for(ID1=0;ID1<node_num;++ID1){
                for(auto it=NeighborMap[ID1].begin();it!=NeighborMap[ID1].end();++it){
                    ID2 = it->first;
                    wei = 1.0/degrees[ID1]+1.0/degrees[ID2];
                    weight = wei * wRange + it->second;
                    if(weight>0){
                        NeighborMap[ID1][ID2] = weight;
                    }else{
                        cout<<"Wrong edge weight! "<<endl;
                        exit(1);
                    }
                }
            }
            break;
        }
        case 4:{
            cout<<"Degree-based ranged random."<<endl;
            vector<int> degrees(node_num,0);
            weight = 1;
            while(getline(inFile,line)){
                if(line.empty())
                    continue;
                vector<string> re1;
                boost::split(re1,line,boost::is_any_of(" \t"));
                if(re1.size()>=2){
                    ID1=stoi(re1[0]), ID2=stoi(re1[1]);
                }
                else{
                    cout<<"Wrong line. "<< line<<endl;
                    exit(1);
                }


                ID1-=1; ID2-=1;
//                weight = 0.05*(rand()%wRange) + 1;

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
                    if(ID1!=ID2){
                        cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
                        exit(1);
                    }
                }
            }

            int maxD = *max_element(degrees.begin(),degrees.end());
//            int range1=maxD/5;
            int deg=0;
            cout<<"Maximum degree: "<<maxD<<endl;
            for(ID1=0;ID1<node_num;++ID1){
                for(auto it=NeighborMap[ID1].begin();it!=NeighborMap[ID1].end();++it){
                    ID2 = it->first;
                    if(degrees[ID1]>=degrees[ID2]){
                        deg = degrees[ID1];
                        if(deg>0.8*maxD){
                            weight = rand()%1000 + 1000;
                        }else if(deg>0.6*maxD){
                            weight = rand()%1000 + 2000;
                        }else if(deg>0.4*maxD){
                            weight = rand()%1000 + 3000;
                        }else if(deg>0.2*maxD){
                            weight = rand()%1000 + 4000;
                        }else{
                            weight = rand()%1000 + 5000;
                        }

                        if(weight>0){
                            NeighborMap[ID1][ID2] = weight;
                            NeighborMap[ID2][ID1] = weight;
                        }else{
                            cout<<"Wrong edge weight! "<<endl;
                            exit(1);
                        }
                    }

                }
            }
            break;
        }
        default:
            break;
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

void Graph::ReadGraph(string filename, bool ifCo){
    string r_edges = filename;
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
    if(ifCo){
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

void Graph::ReadRoadNetworksPre(){
    Timer tt;
    tt.start();
    string gFile=graph_path+"_Distance.gr";
    gFile=graph_path+"_Time.gr";
    string coFile=graph_path+"_Coordinate.co";

    /// read edges
    ifstream inGraph(gFile, ios::in);
    if (!inGraph) { // if not exist
        cout << "Fail to open file " << gFile << endl;
        exit(1);
    }
    cout<<"Reading graph edges. "<<gFile<<endl;
    string line;
    getline(inGraph,line);
    vector<string> vs;
    int ID1,ID2,weight;


    while(!line.empty()){
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs[0]=="c"){
            getline(inGraph,line);
            continue;
        }else if(vs[0]=="p"){
            node_num=stoul(vs[2]), edge_num=stoull(vs[3]);
            cout<<"Node number: "<<node_num<<" ; edge number: "<<edge_num<<endl;
            NeighborMap.assign(node_num,unordered_map<int,int>());
        }else if(vs[0]=="a"){
            if(vs.size()!=4){
                cout<<"Wrong syntax! "<<line<<endl; exit(1);
            }
            ID1=stoi(vs[1]), ID2=stoi(vs[2]), weight=stoi(vs[3]);
            ID1-=1, ID2-=1;
            if(ID1>=0&&ID1<node_num && ID2>=0&&ID2<node_num && weight>0){
                if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//if not found
                    NeighborMap[ID1].insert({ID2,weight});
                    NeighborMap[ID2].insert({ID1,weight});
                }

            }
            else if(weight!=0){
                cout<<"Wrong syntax! "<<line<<endl; exit(1);
            }
        }else{
            cout<<"Wrong syntax! "<<line<<endl; exit(1);
        }
        if(inGraph.eof()) break;
        getline(inGraph,line);
    }
    inGraph.close();
    cout<<"Finished."<<endl;


    /// read coordinates
    ifstream inCoord(coFile, ios::in);
    if (!inCoord) { // if not exist
        cout << "Fail to open file " << coFile << endl;
        exit(1);
    }
    cout<<"Reading graph coordinates. "<<coFile<<endl;

    getline(inCoord,line);
    int co1,co2;

    while(!line.empty()){
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs[0]=="c"){
            getline(inCoord,line);
            continue;
        }else if(vs[0]=="p"){
            vertex nNum=stoul(vs[4]);
            if(nNum!=node_num){
                cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
            }
            Coordinate.assign(node_num,pair<int,int>());
        }else if(vs[0]=="v"){
            if(vs.size()!=4){
                cout<<"Wrong syntax! "<<line<<endl; exit(1);
            }
            ID1=stoi(vs[1]), co1=stoi(vs[2]), co2=stoi(vs[3]);
            ID1-=1;
            if(ID1>=0&&ID1<node_num){
                Coordinate[ID1].first=co1, Coordinate[ID1].second=co2;

            }
            else{
                cout<<"Wrong syntax! "<<line<<endl; exit(1);
            }
        }else{
            cout<<"Wrong syntax! "<<line<<endl; exit(1);
        }
        if(inCoord.eof()) break;
        getline(inCoord,line);
    }
    inCoord.close();
    cout<<"Finished."<<endl;
    tt.stop();
    cout<<"Time for graph reading: "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::ReadRoadNetworksPre2(){
    Timer tt;
    tt.start();
    string gFile=graph_path+"_Distance.gr";
    gFile=graph_path+"_Time.gr";
    string coFile=graph_path+"_Coordinate.co";

    /// read edges
    ifstream inGraph(gFile, ios::in);
    if (!inGraph) { // if not exist
        cout << "Fail to open file " << gFile << endl;
        exit(1);
    }
    cout<<"Reading graph edges."<<endl;
    string line;
    vector<string> vs;
    int ID1,ID2,weight;


    // time graph
    getline(inGraph,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    node_num=stoi(vs[0]), edge_num=stoi(vs[1]);
    cout<<"Node number: "<<node_num<<" , edge number: "<<edge_num<<endl;
    NeighborMap.assign(node_num,unordered_map<int,int>());
    unsigned long long int edgeNum=0;
    while(getline(inGraph,line)){
        if(line.empty()) continue;
        istringstream iss(line);
        if (!(iss >> ID1 >> ID2 >> weight)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

        if(ID1>=0 && ID1<node_num && ID2>=0 && ID2<node_num && weight>0){
//            NeighborMap[ID1].insert({ID2,weight});
            if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//if not found
                NeighborMap[ID1].insert({ID2,weight});
                NeighborMap[ID2].insert({ID1,weight});
                edgeNum+=2;
            }
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<ID2<<" "<<weight<<endl;
        }
    }
    inGraph.close();
    cout<<"Finished."<<endl;
    cout<<"New edge number: "<<edgeNum<<endl;

    /// read coordinates
    ifstream inCoord(coFile, ios::in);
    if (!inCoord) { // if not exist
        cout << "Fail to open file " << coFile << endl;
        exit(1);
    }
    cout<<"Reading graph coordinates."<<endl;


    int co1,co2;

    getline(inCoord,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    int nodeNum=stoi(vs[0]);
    if(nodeNum != node_num){
        cout<<"Inconsistent node number! "<<node_num<<" "<<nodeNum<<endl; exit(1);
    }
    Coordinate.assign(node_num,pair<int,int>());
    while(getline(inCoord,line)){
        if(line.empty()) continue;
        istringstream iss(line);
        if (!(iss >> ID1 >> co1 >> co2)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

        if(ID1>=0 && ID1<node_num ){
            Coordinate[ID1].first=co1, Coordinate[ID1].second=co2;
        }else{
            cout<<"Graph data is wrong! "<<ID1<<" "<<co1<<" "<<co2<<endl;
        }
    }

    inCoord.close();
    cout<<"Finished."<<endl;
    tt.stop();
    cout<<"Time for graph reading: "<<tt.GetRuntime()<<" s."<<endl;

//    ifstream IF(graph_path+"_NodeIDMap", ios::in);
//    if (!IF) { // if not exist
//        cout << "Fail to open file " << endl;
//        exit(1);
//    }
////    vector<long long int> NewToOldID;
//    getline(IF,line);
//    vs.clear();
//    boost::split(vs,line,boost::is_any_of(" "));
//    NewToOldID.assign(stoi(vs[0]),-1);
//    for(int i=0;i<nodeNum;++i){
//        getline(IF,line);
//        vs.clear();
//        boost::split(vs,line,boost::is_any_of(" "));
//        NewToOldID[stoi(vs[1])]= stoll(vs[0]);
//    }
//    IF.close();
//    int tempID=154902;
//    cout<<tempID<<"("<<NewToOldID[tempID]<<")";
//    for(auto it=NeighborMap[tempID].begin();it!=NeighborMap[tempID].end();++it){
//        cout<<" "<<it->first<<"("<<NewToOldID[it->first]<<")";
//    }
//    cout<<endl;
}

//function of checking the connectivity
template <class T>
pair<vertex, edge> Graph::DFS_CC(T & Edges, unordered_set<int> & set_A, set<int> & set_LCC, vertex nodenum) {
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
//    std::sort(CCs.begin(), CCs.end());
//    cout<<"Degree-1 vertex in LCC: ";
//    for(auto it=set_LCC.begin();it!=set_LCC.end();++it){
//        if(NeighborMap[*it].size()==1){
//            cout<<*it<<"("<<NewToOldID[*it]<<") ";
//        }
//    }
//    cout<<endl;
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

    switch (task) {
        case 2:{
            cout<<"Query OD Pair Generation."<<endl;
            ODpairGenerate(100000);
            break;
        }
        case 3:{
            cout<<"Update OD Pair Generation!"<<endl;
            UpdateGenerate(100000);
//            UpdateGenerate2(1000);
            break;
        }
        case 4:{
            cout<<"Vertex Ordering Generation!"<<endl;
            VertexOrdering(strategy);//1: MDE-based ordering; 2: Degree-based ordering
            break;
        }
        case 5:{
            cout<<"Preprocessing for road network!"<<endl;
            RoadNetworkPreprocessing();
            break;
        }
        case 6:{
            cout<<"Compute the diameters of road network!"<<endl;
            ComputeDiameterOfNetwork(graph_path+".co");
            break;
        }
        case 7:{
            cout<<"Extract subgraph from the network!"<<endl;
            if(dataset=="NY" || dataset=="FLA" || dataset=="W" || dataset=="USA"){
                ReadGraph(graph_path,true);
                ExtractSubgraphForScalabilityTest(true);
            }
            else{
                ReadGraph(graph_path,false);
                ExtractSubgraphForScalabilityTest(false);
            }

            break;
        }
        case 8:{
            cout<<"Generation of queries with varying distance."<<endl;

            if(dataset=="NY" || dataset=="FLA" || dataset=="W" ||dataset=="USA"){
                cout<<"Spatial distance."<<endl;
                ODpairGenerateDistanceSpatial(10000);
            }else{
                cout<<"Shortest-path distance."<<endl;
                ODpairGenerateDistanceDijk(10000);
            }

            break;
        }
        case 9:{
            cout<<"Maximal step number of all shortest paths."<<endl;
            ComputeMaximalStepNumOfShortestPath();

            break;
        }
        default:
            break;
    }
}



//Function for generating random queries
void Graph::ODpairGenerate(int times){
    string RandomDis = graph_path + ".query";

    ifstream inFile(graph_path, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << graph_path << endl;
        exit(1);
    }

    /// read graph and recording the degree of vertices
    inFile >> node_num >> edge_num;
    inFile.close();
    cout<<"Node number: "<<node_num<<" ;  edge number: "<<edge_num<<endl;
    /*---OD pairs generation---*/
    int pairs = 0;
    int node_start, node_end;
    int temp = 0;
    ofstream outFile(RandomDis, ios::out);
    if (!outFile) {
        cout << "File opening failed." << endl;
        assert(outFile);
    }
    cout << "Query OD pairs file generating... "<<times << endl;
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
//Function for generating random queries with varying distance
void Graph::ODpairGenerateDistanceDijk(int times){

    ReadGraph(graph_path, false);

    if(node_num>1000000){
        times=1000;
    }
    int ID1, ID2;
    vector<vector<pair<int,int>>> Query;
    Query.assign(10,vector<pair<int,int>>());
    vector<pair<double,double>> disRange;
    int maxDis=0;//maximal distance
    cout<<"Query number: "<<times<<endl;
    Timer tt;
    tt.start();
    // Step 1: compute the maximal shortest path length
    for(int i=0;i<10;++i){
        ID1=rand()%node_num;

        benchmark::heap<2, int, int> pqueue(node_num);
        pqueue.update(ID1,0);
        vector<bool> closed(node_num, false);
        vector<int> distance(node_num, INF);
        distance[ID1]=0;
        int topNodeID, topNodeDis;

        int NNodeID,NWeigh;

        int d=0;//maximal distance

        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            if(d<topNodeDis){
                d=topNodeDis;
            }
            closed[topNodeID]=true;

            for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
                NNodeID=(*it).first;
                NWeigh=(*it).second+topNodeDis;
                if(!closed[NNodeID]){
                    if(distance[NNodeID]>NWeigh){
                        distance[NNodeID]=NWeigh;
                        pqueue.update(NNodeID, NWeigh);
                    }
                }
            }
        }

        if(maxDis<d){
            maxDis=d;
        }
    }

    cout<<"Simulated maximal shortest path length: "<<maxDis<<endl;

    double disUnite=pow(maxDis/1000,0.1);
    cout<<"distance range unite: "<<disUnite<<endl;


    for(int i=0;i<10;++i){
        disRange.emplace_back(1000*pow(disUnite,i), 1000*pow(disUnite,i+1));
    }
    vector<vector<int>> candidates;
    // Step 2: Obtain the queries
    for(int i=0;i<times;++i){
        ID1=rand()%node_num;
        candidates.assign(10,vector<int>());
        benchmark::heap<2, int, int> pqueue(node_num);
        pqueue.update(ID1,0);
        vector<bool> closed(node_num, false);
        vector<int> distance(node_num, INF);
        distance[ID1]=0;
        int topNodeID, topNodeDis;

        int NNodeID,NWeigh;

        int d=0;//maximal distance
        int graph_i=0;
        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            if(topNodeDis>disRange[graph_i].first && topNodeDis<=disRange[graph_i].second){
                candidates[graph_i].push_back(topNodeID);
            }else if(topNodeDis>disRange[graph_i].second){
                graph_i++;
                if(graph_i>9){
                    graph_i=9;
                }
            }
            closed[topNodeID]=true;

            for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
                NNodeID=(*it).first;
                NWeigh=(*it).second+topNodeDis;
                if(!closed[NNodeID]){
                    if(distance[NNodeID]>NWeigh){
                        distance[NNodeID]=NWeigh;
                        pqueue.update(NNodeID, NWeigh);
                    }
                }
            }
        }

        for(int j=0;j<10;++j){
            if(!candidates[j].empty()) {
                ID2=candidates[j][rand()%candidates[j].size()];
                Query[j].emplace_back(ID1,ID2);
            }
            else{
//                cout<<"Q"<<j<<" does not have candidate vertex! "<<endl; exit(1);
            }

        }
    }
    for(int i=0;i<10;++i){
        if(Query[i].size()<times){
            cout<<"Size Q"<<i+1<<": "<<Query[i].size()<<endl;
            while(Query[i].size()<times){
                ID1=rand()%node_num;
                candidates.assign(10,vector<int>());
                benchmark::heap<2, int, int> pqueue(node_num);
                pqueue.update(ID1,0);
                vector<bool> closed(node_num, false);
                vector<int> distance(node_num, INF);
                distance[ID1]=0;
                int topNodeID, topNodeDis;

                int NNodeID,NWeigh;

                int d=0;//maximal distance
                int graph_i=0;
                while(!pqueue.empty()){
                    pqueue.extract_min(topNodeID, topNodeDis);
                    if(topNodeDis>disRange[graph_i].first && topNodeDis<=disRange[graph_i].second){
                        candidates[graph_i].push_back(topNodeID);
                    }else if(topNodeDis>disRange[graph_i].second){
                        graph_i++;
                        if(graph_i>9){
                            graph_i=9;
                        }
                    }
                    closed[topNodeID]=true;

                    for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
                        NNodeID=(*it).first;
                        NWeigh=(*it).second+topNodeDis;
                        if(!closed[NNodeID]){
                            if(distance[NNodeID]>NWeigh){
                                distance[NNodeID]=NWeigh;
                                pqueue.update(NNodeID, NWeigh);
                            }
                        }
                    }
                }


                if(!candidates[i].empty()) {
                    ID2=candidates[i][rand()%candidates[i].size()];
                    Query[i].emplace_back(ID1,ID2);
                }


            }
        }
    }
    // Step 3: write results
    cout<<"Writing results..."<<endl;
    for(int i=0;i<10;++i){
        string filename = graph_path + "_SP.Q"+ to_string(i+1);
        ofstream OF(filename);
        if(!OF.is_open()){
            cout<<"Cannot open "<<filename<<endl; exit(1);
        }
//        cout<<"Query number of "<<i<<" : "<<Query[i].size()<<endl;
        OF<<Query[i].size()<<endl;
        for(auto it=Query[i].begin();it!=Query[i].end();++it){
            OF<<it->first<<" "<<it->second<<endl;
        }
        OF.close();
    }


    cout << "Finished." << endl;
    tt.stop();
    cout<<"Overall time: "<<tt.GetRuntime()<<" s."<<endl;
}
void Graph::ODpairGenerateDistanceSpatial(int times){
    ReadGraph(graph_path, false);

    int ID1, ID2;
    vector<vector<pair<int,int>>> Query;
    Query.assign(10,vector<pair<int,int>>());
    vector<pair<double,double>> disRange;
    cout<<"Query number: "<<times<<endl;
    Timer tt;
    tt.start();
    // Step 1: compute the maximal shortest path length
    ComputeDiameterOfNetwork(graph_path+".co");

    double disUnite=pow(graphDiameter/1,0.1);
    cout<<"distance range unite: "<<disUnite<< endl;

    for(int i=0;i<10;++i){
        disRange.emplace_back(1*pow(disUnite,i), 1*pow(disUnite,i+1));
        cout<<"Distance range of Q"<<i+1<<": ("<<disRange[i].first<<","<<disRange[i].second<<"]"<<endl;
    }
    vector<vector<int>> candidates;
    // Step 2: Obtain the queries
    double dis;
    for(int i=0;i<times;++i){
        ID1=rand()%node_num;
        candidates.assign(10,vector<int>());
        for(ID2=0;ID2<node_num;++ID2){
            if(ID2!=ID1){
                dis=EuclideanDis(make_pair((double)Coordinate[ID1].first/1000000, (double)Coordinate[ID1].second/1000000) , make_pair((double)Coordinate[ID2].first/1000000,(double)Coordinate[ID2].second/1000000));
                bool ifFind=false;
                for(int j=0;j<10;++j){
                    if(dis>disRange[j].first && dis<=disRange[j].second){
                        candidates[j].push_back(ID2);
                        ifFind=true;
                        break;
                    }
                }
                if(!ifFind){
//                    cout<<"Unsuccessful OD pair. "<<ID1<<" "<<ID2<<" "<<dis<<endl;
//                    exit(1);
                }
            }
        }

        for(int j=0;j<10;++j){
            if(!candidates[j].empty()) {
                ID2=candidates[j][rand()%candidates[j].size()];
                Query[j].emplace_back(ID1,ID2);
            }
            else{
//                cout<<"Q"<<j<<" does not have candidate vertex! "<<endl; exit(1);
            }

        }

    }

    for(int i=0;i<10;++i){
        if(Query[i].size()<times){
            cout<<"Size Q"<<i+1<<": "<<Query[i].size()<<endl;
            while(Query[i].size()<times){
                ID1=rand()%node_num, ID2=rand()%node_num;
                dis=EuclideanDis(make_pair((double)Coordinate[ID1].first/1000000, (double)Coordinate[ID1].second/1000000) , make_pair((double)Coordinate[ID2].first/1000000,(double)Coordinate[ID2].second/1000000));
                if(dis>disRange[i].first && dis<=disRange[i].second){
                    Query[i].emplace_back(ID1,ID2);
                }

            }
        }
    }

    // Step 3: write results
    cout<<"Writing results..."<<endl;
    for(int i=0;i<10;++i){
        string filename = graph_path + "_spatial.Q"+ to_string(i+1);
        ofstream OF(filename);
        if(!OF.is_open()){
            cout<<"Cannot open "<<filename<<endl; exit(1);
        }
        cout<<"Query number of "<<i<<" : "<<Query[i].size()<<endl;
        OF<<Query[i].size()<<endl;
        for(auto it=Query[i].begin();it!=Query[i].end();++it){
            ID1=it->first, ID2=it->second;
            dis= EuclideanDis(make_pair((double)Coordinate[ID1].first/1000000, (double)Coordinate[ID1].second/1000000) , make_pair((double)Coordinate[ID2].first/1000000,(double)Coordinate[ID2].second/1000000));
            if(dis<disRange[i].first){
                cout<<"Wrong query of Q"<<i<<": "<<ID1<<" "<<ID2<<" "<<dis<<" "<<disRange[i].first<<endl;
            }
            OF<<it->first<<" "<<it->second<<endl;
        }
        OF.close();
    }


    cout << "Finished." << endl;
    tt.stop();
    cout<<"Overall time: "<<tt.GetRuntime()<<" s."<<endl;
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
    cout<<node_num<<" "<<edge_num<<endl;
    if(edge_num<200000){
        times=1000;
    }
    cout<<"update number: "<<times<<endl;
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
//Function for generating update edges
void Graph::UpdateGenerate2(int times){
    string filename = graph_path + ".update2";
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
    int oldW,newW1,newW2;
    double prop;
    vector<pair<pair<int,int>,int>> uEdge;
    while (pairs < times) {
        id = rand() % size;
        prop = 0.01*(rand()%100+1);//[0.01,1]
        oldW = edges[id].second;
        newW1 = oldW * (1.0 - 0.5*prop);//[0.5,0.995]
        newW2 = oldW * (1.0 + 0.5*prop);//[1.005,1.5]
        if(newW1 >= 2 && edgeIdSet.find(id)==edgeIdSet.end()){//if edge weight is no smaller than 2, and it has not been added.
            edgeIdSet.insert(id);
            outFile << edges[id].first.first << ' ' << edges[id].first.second << ' '<<oldW<<' '<<newW1<<' '<<newW2<< endl;
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
    ReadGraph(graph_path, false);
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
        for(int i = 0; i < NodeOrder.size(); i++){
            ofile << i << "\t" << NodeOrder[i] << endl;//ID, order
        }

        ofile.close();
    }
}

void Graph::RoadNetworkPreprocessing(){

    ReadRoadNetworksPre();//for DIMACS data
//    ReadRoadNetworksPre2();//for private data

    unordered_set<int> vertices; vertices.clear();
    for(int i=0;i<node_num;++i){
        vertices.insert(i);
    }
    set<int> verticesFinal;
    pair<vertex, edge> LCC;
    LCC = DFS_CC(NeighborMap,vertices,verticesFinal,node_num);
    map<int,int> IDMap;
    int ID=0;
    for(auto it=verticesFinal.begin();it!=verticesFinal.end();++it){
        IDMap.insert({*it,ID});
        ++ID;
    }
    if(IDMap.size()!=LCC.first){
        cout<<"Wrong! Inconsistent! "<<IDMap.size()<<" "<<LCC.first<<endl; exit(1);
    }
    /// Write ID map
//    ifstream IF(graph_path+"_NodeIDMap", ios::in);
//    if (!IF) { // if not exist
//        cout << "Fail to open file " << endl;
//        exit(1);
//    }
//    string line;
//    vector<string> vs;
//    getline(IF,line);
//    vs.clear();
//    boost::split(vs,line,boost::is_any_of(" "));
//    int tempNum=stoi(vs[0]);
//    OldToNewID.clear();
//    long long int ID1;
//    int ID2;
//    for(int i=0;i<tempNum;++i){
//        getline(IF,line);
//        vs.clear();
//        boost::split(vs,line,boost::is_any_of(" "));
//        ID1=stoll(vs[0]);
//        ID2=stoi(vs[1]);
//        if(OldToNewID.find(ID1)==OldToNewID.end()){
//            OldToNewID.insert({ID1,ID2});
//        }else{
//            cout<<"Wrong ID. exist! "<<ID1<<endl; exit(1);
//        }
//    }
//    IF.close();
    ofstream OF(graph_path+".IDMap", ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<endl; exit(1);
    }
    OF<<IDMap.size()<<endl;
    for(auto it=IDMap.begin();it!=IDMap.end();++it){
        OF<<it->first<<" "<<it->second<<endl;//from old id to new id of LCC
    }
    OF.close();
    /// Write graph edges
    string wGraph=graph_path;
    wGraph=graph_path+".dis";//+"2";
    wGraph=graph_path+".time";//+"2";
    ofstream outGraph(wGraph, ios::out);
    if(!outGraph.is_open()){
        cout<<"Cannot open file "<<wGraph<<endl; exit(1);
    }
    cout<<"Writing graph edges. "<<wGraph<<endl;
    outGraph<<LCC.first<<" "<<LCC.second<<endl;
    int ID1,ID2,weight;
    for(auto it=verticesFinal.begin();it!=verticesFinal.end();++it){
        ID1=*it;
        for(auto it2=NeighborMap[ID1].begin();it2!=NeighborMap[ID1].end();++it2){
            ID2=it2->first, weight=it2->second;
            outGraph<<IDMap[ID1]<<" "<<IDMap[ID2]<<" "<<weight<<endl;
        }
    }
    outGraph.close();
    cout<<"Finished."<<endl;

    /// Write graph coordinates
    string wCoord=graph_path+".dis.co";
    wCoord=graph_path+".time.co";
    ofstream outCoord(wCoord, ios::out);
    if(!outCoord.is_open()){
        cout<<"Cannot open file "<<wCoord<<endl; exit(1);
    }
    cout<<"Writing graph coordinates. "<<wCoord <<endl;
    outCoord<<LCC.first<<endl;
    pair<double,double> maxLon=make_pair(-999999999,0);
    pair<double,double> maxLat=make_pair(0,-999999999);
    pair<double,double> minLon=make_pair(999999999,0);
    pair<double,double> minLat=make_pair(0,999999999);
    for(auto it=verticesFinal.begin();it!=verticesFinal.end();++it){
        ID1=*it;
        outCoord<<IDMap[ID1]<<" "<<Coordinate[ID1].first<<" "<<Coordinate[ID1].second<<endl;
        if(Coordinate[ID1].first>maxLon.first) {
            maxLon.first=Coordinate[ID1].first; maxLon.second=Coordinate[ID1].second;
        }
        if(Coordinate[ID1].second>maxLat.second) {
            maxLat.second=Coordinate[ID1].second; maxLat.first=Coordinate[ID1].first;
        }
        if(Coordinate[ID1].first<minLon.first) {
            minLon.first=Coordinate[ID1].first; minLon.second=Coordinate[ID1].second;
        }
        if(Coordinate[ID1].second<minLat.second) {
            minLat.second=Coordinate[ID1].second; minLat.first=Coordinate[ID1].first;
        }
    }
    outCoord.close();
    minLon.first/=1000000, minLon.second/=1000000, maxLon.first/=1000000, maxLon.second/=1000000;
    minLat.first/=1000000, minLat.second/=1000000, maxLat.first/=1000000, maxLat.second/=1000000;
    cout<<"Longitude range: "<<minLon.first<<" "<<maxLon.first<<" ; Latitude range: "<<minLat.second<<" "<<maxLat.second<<endl;
    vector<pair<double,double>> points;
    points.push_back(minLon); points.push_back(maxLon); points.push_back(minLat); points.push_back(maxLat);
    cout << "Maximal diameter: "<< ComputeDiameter(points)<< " km."<<endl;
    cout<<"Finished."<<endl;
}

void Graph::ComputeMaximalStepNumOfShortestPath(){
    string filename=graph_path+".co";
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    cout<<"Read file "<<filename<<endl;
    string line;
    int ID,x,y;
    int nodeNum=0;
    bool flag= false;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> node_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
//    if(nodeNum != node_num){
//        cout<<"Nominate node number of Coordinate file: "<<nodeNum<<" ; real node number: "<<node_num<<endl;
//        flag=true;
//    }
    Coordinate.assign(node_num,pair<int,int>());
    pair<double,double> maxLon=make_pair(-999999999,0);
    pair<double,double> maxLat=make_pair(0,-999999999);
    pair<double,double> minLon=make_pair(999999999,0);
    pair<double,double> minLat=make_pair(0,999999999);
    unordered_set<int> maxLonV;
    unordered_set<int> maxLatV;
    unordered_set<int> minLonV;
    unordered_set<int> minLatV;
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
            ++nodeNum;

            if(Coordinate[ID].first>maxLon.first) {
                maxLon.first=Coordinate[ID].first; maxLon.second=Coordinate[ID].second;
                maxLonV.clear(); maxLonV.insert(ID);
            }else if(Coordinate[ID].first==maxLon.first){
                maxLonV.insert(ID);
            }
            if(Coordinate[ID].second>maxLat.second) {
                maxLat.second=Coordinate[ID].second; maxLat.first=Coordinate[ID].first;
                maxLatV.clear(); maxLatV.insert(ID);
            }else if(Coordinate[ID].second==maxLat.second){
                maxLatV.insert(ID);
            }
            if(Coordinate[ID].first<minLon.first) {
                minLon.first=Coordinate[ID].first; minLon.second=Coordinate[ID].second;
                minLonV.clear(); minLonV.insert(ID);
            }else if(Coordinate[ID].first==minLon.first){
                minLonV.insert(ID);
            }
            if(Coordinate[ID].second<minLat.second) {
                minLat.second=Coordinate[ID].second; minLat.first=Coordinate[ID].first;
                minLatV.clear(); minLatV.insert(ID);
            }else if(Coordinate[ID].second==minLat.second){
                minLatV.insert(ID);
            }
        }else{
            cout<<"Graph data is wrong! "<<ID<<" "<<x<<" "<<y<<endl; exit(1);
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

    minLon.first/=1000000, minLon.second/=1000000, maxLon.first/=1000000, maxLon.second/=1000000;
    minLat.first/=1000000, minLat.second/=1000000, maxLat.first/=1000000, maxLat.second/=1000000;
    cout<<"Longitude range: "<<minLon.first<<" "<<maxLon.first<<" ; Latitude range: "<<minLat.second<<" "<<maxLat.second<<endl;
    vector<pair<double,double>> points;
    points.push_back(minLon); points.push_back(maxLon); points.push_back(minLat); points.push_back(maxLat);
    graphDiameter=ComputeDiameter(points);
    cout << "Maximal diameter: "<< graphDiameter<< " km."<<endl;


    ReadGraph(graph_path, false);


    maxLatV.insert(maxLonV.begin(),maxLonV.end());
    maxLatV.insert(minLonV.begin(),minLonV.end());
    maxLatV.insert(minLatV.begin(),minLatV.end());

    for(int i=0;i<50;++i){
        int ID=rand()%node_num;
        if(maxLatV.find(ID)==maxLatV.end()){
            maxLatV.insert(ID);
        }

    }
    ComputeLongestShortestPath(maxLatV);

}

pair<unsigned long int, unsigned long int> Graph::DijkstraSSSP(int ID){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> pre(node_num, 0);
    vector<int> step(node_num, 0);
    distance[ID]=0;
    pre[ID]=ID;
    step[ID]=-1;
    int topNodeID, topNodeDis;

    int NNodeID,NWeigh;

    unsigned long int d=0;//initialize d to infinite for the unreachable case
    unsigned long int stepNum=0;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        step[topNodeID]=step[pre[topNodeID]]+1;
        if(stepNum<step[topNodeID]){
            stepNum=step[topNodeID];
        }
        if(d<topNodeDis){
            d=topNodeDis;
        }
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    pre[NNodeID]=topNodeID;
                }
            }
        }
    }

//    if(stepNum!=step[topNodeID]){
//        cout<<"maximal step number is not the last vertex. "<<topNodeID<<"("<<d<<") "<<step[topNodeID]<<" "<<stepNum<<endl;
//    }

    return make_pair(d,stepNum);
}

void Graph::ComputeLongestShortestPath(unordered_set<int> & sourceIDs){
    cout<<"source vertex number: "<<sourceIDs.size()<<endl;
    unsigned long int maxDis=0;
    unsigned long int maxStep=0;
    Timer tt;
    tt.start();
    for(auto it=sourceIDs.begin();it!=sourceIDs.end();++it){
        auto re1 = DijkstraSSSP(*it);
        if(maxDis < re1.first){
            maxDis = re1.first;
        }
        if(maxStep < re1.second){
            maxStep = re1.second;
        }
    }
    tt.stop();
    cout<<"For all shortest paths, maximal distance: "<<(double)maxDis/1000<<" km ; maximal step number: "<<maxStep<<" ; average query time: "<< 1000*tt.GetRuntime()/sourceIDs.size()<<" ms."<<endl;

}

void Graph::ComputeDiameterOfNetwork(string filename){
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    cout<<"Read file "<<filename<<endl;
    string line;
    int ID,x,y;
    int nodeNum=0;
    bool flag= false;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> node_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
//    if(nodeNum != node_num){
//        cout<<"Nominate node number of Coordinate file: "<<nodeNum<<" ; real node number: "<<node_num<<endl;
//        flag=true;
//    }
    Coordinate.assign(node_num,pair<int,int>());
    pair<double,double> maxLon=make_pair(-999999999,0);
    pair<double,double> maxLat=make_pair(0,-999999999);
    pair<double,double> minLon=make_pair(999999999,0);
    pair<double,double> minLat=make_pair(0,999999999);
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
            ++nodeNum;

            if(Coordinate[ID].first>maxLon.first) {
                maxLon.first=Coordinate[ID].first; maxLon.second=Coordinate[ID].second;
            }
            if(Coordinate[ID].second>maxLat.second) {
                maxLat.second=Coordinate[ID].second; maxLat.first=Coordinate[ID].first;
            }
            if(Coordinate[ID].first<minLon.first) {
                minLon.first=Coordinate[ID].first; minLon.second=Coordinate[ID].second;
            }
            if(Coordinate[ID].second<minLat.second) {
                minLat.second=Coordinate[ID].second; minLat.first=Coordinate[ID].first;
            }
        }else{
            cout<<"Graph data is wrong! "<<ID<<" "<<x<<" "<<y<<endl; exit(1);
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

    minLon.first/=1000000, minLon.second/=1000000, maxLon.first/=1000000, maxLon.second/=1000000;
    minLat.first/=1000000, minLat.second/=1000000, maxLat.first/=1000000, maxLat.second/=1000000;
    cout<<"Longitude range: "<<minLon.first<<" "<<maxLon.first<<" ; Latitude range: "<<minLat.second<<" "<<maxLat.second<<endl;
    vector<pair<double,double>> points;
    points.push_back(minLon); points.push_back(maxLon); points.push_back(minLat); points.push_back(maxLat);
    graphDiameter=ComputeDiameter(points);
    cout << "Maximal diameter: "<< graphDiameter<< " km."<<endl;

}

double Graph::ComputeDiameter(vector<pair<double,double>> & points){
    double maxDis=0;
    for(int i=0;i<points.size();++i){
        for(int j=i+1;j<points.size();++j){
            double temp = EuclideanDis(points[i],points[j]);
            if(temp>maxDis) maxDis=temp;
        }

    }
    return maxDis;
}

void Graph::ExtractSubgraphForScalabilityTest(bool ifCo){
    int sourceID=rand()%node_num;

    vector<vector<int>> targetVertex(4);
    vector<bool> visited(node_num,false);
    int targetNum=node_num*0.2;

    queue<int> BFSQueue;
    BFSQueue.push(sourceID);
    visited[sourceID]=true;
    int graph_i=0;
    while(!BFSQueue.empty()){
        int topID=BFSQueue.front();
        BFSQueue.pop();
        targetVertex[graph_i].push_back(topID);
        if(targetVertex[graph_i].size()>=targetNum){
            graph_i++;
            if(graph_i==4){
                break;
            }
        }
        for(auto it=Neighbor[topID].begin();it!=Neighbor[topID].end();++it){
            if(!visited[it->first]){
                BFSQueue.push(it->first);
                visited[it->first]=true;
            }
        }
    }
    vector<set<int>> graphVertex(4);
    graphVertex[0].insert(targetVertex[0].begin(),targetVertex[0].end());
    for(int i=1;i<4;i++){
        graphVertex[i]=graphVertex[i-1];
        graphVertex[i].insert(targetVertex[i].begin(),targetVertex[i].end());
    }
    for(int i=0;i<4;++i){
        string filename=graph_path+"_"+to_string(20*(i+1));
        ofstream OF(filename);
        if(!OF.is_open()){
            cout<<"Cannot open "<<filename<<endl; exit(1);
        }
        int graphVNum=graphVertex[i].size();

        map<int,int> IDMap;
        vector<map<int,int>> NeighborTemp(graphVNum,map<int,int>());
        int ID1, ID2, weight;
        int newID=0;
        for(auto it=graphVertex[i].begin();it!=graphVertex[i].end();++it) {
            ID1 = *it;
            IDMap.insert({ID1,newID});
            newID++;
        }
        unsigned long long int eNum=0;
        for(auto it=graphVertex[i].begin();it!=graphVertex[i].end();++it){
            ID1=*it;
            for(auto it2=Neighbor[ID1].begin();it2!=Neighbor[ID1].end();++it2){
                ID2=it2->first, weight=it2->second;
                if(graphVertex[i].find(ID2)!=graphVertex[i].end()){//if found
                    NeighborTemp[IDMap[ID1]][IDMap[ID2]]=weight;
                    eNum++;
                }
            }
        }
        cout<<"Graph "<<i<<" : vertex number "<<graphVNum<<" ; edge number "<<eNum<<endl;
        OF<<graphVNum<<" "<<eNum <<endl;
        for(int i=0;i<NeighborTemp.size();++i){
            for(auto it=NeighborTemp[i].begin();it!=NeighborTemp[i].end();++it){
                OF<<i<<" "<<it->first<<" "<<it->second<<endl;
            }
        }
        OF.close();

        if(ifCo){
            string filename2=filename+".co";
            ofstream OF2(filename2);
            if(!OF2.is_open()){
                cout<<"Cannot open "<<filename2<<endl; exit(1);
            }
            OF2<<graphVNum<<endl;
            for(auto it=graphVertex[i].begin();it!=graphVertex[i].end();++it){
                ID1=*it;
                OF2<<IDMap[ID1]<<" "<<Coordinate[ID1].first<<" "<<Coordinate[ID1].second<<endl;
            }
            OF2.close();
        }
        cout<<"Write done."<<endl;
    }

}

//function of computing Euclidean distance
double Graph::EuclideanDis(pair<double,double> s, pair<double,double> t)
{
    double d=111.319;//distance in equator per degree, kilometer
    double coe;
    double temp=(s.second+t.second)/2;
    temp=temp*3.1415926/180;
//    temp=temp*3.1415926/(180*1000000);
    coe=cos(temp);
    double y1 = s.second * d, y2 = t.second * d;
    double x1 = s.first * d * coe, x2 = t.first * d * coe;
    double xx = x1 - x2, yy = y1 - y2;
    return sqrt(xx * xx + yy * yy);//Euclidean Distance in meter
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