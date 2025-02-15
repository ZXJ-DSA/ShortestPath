//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef PLL_HPP
#define PLL_HPP

#include "head.h"

/*-------------Index-based Algorithms--------------*/
/// PLL
//PLL construction entry
void Graph::PLLConstruction(int strategy, bool ifDynamic){
    string indexfile=graph_path+".PLL";
    string PPRfile=graph_path+".PPR";

    VertexOrdering(strategy);


    if(ifDynamic){
        cout<<"For dynamic graph."<<endl;
    }else{
        cout<<"For static graph."<<endl;
    }


    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;
    bool ifHash= true;
    ifHash=false;

    if(ifHash){
        t1=std::chrono::high_resolution_clock::now();
        PLLConstruct();
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"PLL construction time: "<<runT<<" s."<<endl;
    }else{
        vector<vertex> vertices;
        for(int i=0;i<node_num;++i){
            vertices.emplace_back(i);
        }
        t1=std::chrono::high_resolution_clock::now();
        PLLConstructV(Neighbor);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"PLL construction time: "<<runT<<" s."<<endl;

        LabelV.postProcess(Label,threadnum);
        if(ifDynamic){
            PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
        }

    }

    IndexSizePLL();
}
//PLL construction function
void Graph::PLLConstruct(){
    cout<<"Hash table-based Construction."<<endl;
    Label.assign(node_num,unordered_map<vertex,int>());

    if(ifDynamic){
        PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
    }

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/10000)*2000;
    stepShow = max(stepShow,1000);
    cout<<"Step for show: "<<stepShow<<endl;

    Timer tt;
    tt.start();
    for(int i=node_num-1;i>=0;i--){//start from the least important vertex

        ID=vNodeOrder[i];

        vector<pair<int,int>> vp;
        DijksPrune1(ID,vp);
        if(cnt%stepShow==0){
            tt.stop();
            cout<<"Node "<<cnt<<": "<<ID<<" ; vp.size: "<<vp.size()<<"; time: "<<tt.GetRuntime()<<" s."<<endl;
            tt.start();
            //cout<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
        }
        cnt+=1;
    }

}

void Graph::PLLConstructV(vector<vector<pair<vertex,int>>>& Neighbor){
    cout<<"Vector-based construction."<<endl;

    LabelV.resize(node_num);

    if(ifDynamic){
        PPRV.resize(node_num);
    }


    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/10000)*2000;
    stepShow = max(stepShow,1000);
    cout<<"Step for show: "<<stepShow<<endl;

    Timer tt;
    tt.start();
    for(int i=node_num-1;i>=0;i--){//start from the least important vertex

        ID=vNodeOrder[i];

        DijksPrune1V(ID,Neighbor);
        if(cnt%stepShow==0){
            tt.stop();
            cout<<"Node "<<cnt<<": "<<ID<<" ; time: "<<tt.GetRuntime()<<" s."<<endl;
            tt.start();
            //cout<<"ID "<<ID<<" vp.size "<<vp.size()<<endl;
        }
        cnt+=1;
    }

    LabelV.sort(threadnum);
}
//pruned Dijkstra from nodeID
void Graph::DijksPrune1(int nodeID, vector<pair<int,int>>& vp){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    if(ifDynamic){
        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            closed[topNodeID]=true;

            int TempDis; vector<int> SupNode;
            PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
            if(TempDis<=topNodeDis){
                //for index update
            if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
                for(int k=0;k<SupNode.size();k++){
                    int supn=SupNode[k];
                    PruningPointSet[topNodeID][supn].insert(nodeID);
                    PruningPointSet[nodeID][supn].insert(topNodeID);
                }
            }
                continue;
            }

            Label[topNodeID].insert({nodeID, topNodeDis});
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
    }else{
        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            closed[topNodeID]=true;

            int TempDis; vector<int> SupNode;
            PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
            if(TempDis<=topNodeDis){
                continue;
            }

            Label[topNodeID].insert({nodeID, topNodeDis});
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
    }

}
//vector-based implementation
void Graph::DijksPrune1V(int nodeID, vector<vector<pair<vertex,int>>>& Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    unordered_map<vertex,int> Lh;
    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
        Lh.insert({it->first,it->second});
    }

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    if(ifDynamic){
        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            closed[topNodeID]=true;

            ull TempDis; vector<int> SupNode;
//        PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
            TempDis = PLLDisQuery1V(topNodeID,Lh,SupNode);

            if(TempDis<=topNodeDis){
                //for index update
                if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
                    for(int k=0;k<SupNode.size();k++){
                        int supn=SupNode[k];

                        PPRV.add(nodeID,supn,topNodeID);
                    }
                }
                continue;
            }

            LabelV.add(topNodeID,nodeID,topNodeDis);
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
    }else{
        while(!pqueue.empty()){
            pqueue.extract_min(topNodeID, topNodeDis);
            closed[topNodeID]=true;

            ull TempDis; vector<int> SupNode;
//        PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
            TempDis = PLLDisQuery1V(topNodeID,Lh,SupNode);

            if(TempDis<=topNodeDis){
                continue;
            }

            LabelV.add(topNodeID,nodeID,topNodeDis);
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
    }

}
//query processing with current PLL index
int Graph::PLLDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;

    int hub, dis1, dis2;
    if(ifDynamic){
        for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(Label[ID2].find(hub)!=Label[ID2].end()){
                dis2=Label[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                    SupNode.clear();
                    SupNode.push_back(hub);
                }
                else if(dis1+dis2==d){
                    SupNode.push_back(hub);
                }
            }
        }
    }else{
        for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(Label[ID2].find(hub)!=Label[ID2].end()){
                dis2=Label[ID2][hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                }
            }
        }
    }


    return d;
}
//vector-based implementation
int Graph::PLLDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode){
    int d=INF;

    unsigned int hub;
    unsigned long long dis1, dis2;
    if(ifDynamic){
        for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(Lh.find(hub)!=Lh.end()){
                dis2=Lh[hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                    SupNode.clear();
                    SupNode.push_back(hub);
                }
                else if(dis1+dis2==d){
                    SupNode.push_back(hub);
                }
            }
        }
    }else{
        for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();it++){
            hub=(*it).first;
            dis1=(*it).second;
            if(Lh.find(hub)!=Lh.end()){
                dis2=Lh[hub];
                if(dis1+dis2<d){
                    d=dis1+dis2;
                }
            }
        }
    }


    return d;
}

void Graph::writePLL(string filename, string filenameP){
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open "<<filename<<endl;
        exit(1);
    }
    cout<<"Writing index to disk..."<<endl;
    OF<<node_num<<endl;
    for(int nodeID=0;nodeID<node_num;nodeID++){
        OF<<nodeID<<" "<<NodeOrder[nodeID]<<" "<<Label[nodeID].size();
        for(auto it=Label[nodeID].begin();it!=Label[nodeID].end();it++){
            OF<<" "<<(*it).first<<" "<<(*it).second;
        }
        OF<<endl;
    }
    OF.close();
//    ofstream OF1(filenameP);
//    for(int nodeID=0;nodeID<node_num;nodeID++){
//        OF1<<nodeID<<" "<<PruningPointNew[nodeID].size();///PruningPointNew
//        for(auto itp=PruningPointNew[nodeID].begin();itp!=PruningPointNew[nodeID].end();itp++){
//            OF1<<" "<<(*itp).first<<" "<<(*itp).second.size();
////			for(int k=0;k<(*itp).second.size();k++){
////				OF1<<" "<<(*itp).second[k];
////			}
//            for(auto k=itp->second.begin();k!=itp->second.end();k++){
//                int temp1=*k;
//                OF1<<" "<< temp1;
//            }
//        }
//        OF1<<endl;
//    }
//    OF1.close();
    cout<<"Finish index writing!"<<endl;
}
void Graph::readPLL(string filename, string filenameP){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open "<<filename<<endl;
        exit(1);
    }
    cout<<"Reading PLL index..."<<endl;
    int num;
    IF>>num;

    Label.assign(num,unordered_map<vertex,int>());
    NodeOrder.assign(num,0);
    vNodeOrder.assign(num,0);

    int nodeID, order, pnum, hubid, dis;
    for(int i=0;i<num;i++){
        IF>>nodeID>>order>>pnum;
        NodeOrder[nodeID]=order;
        vNodeOrder[order]=nodeID;
        for(int j=0;j<pnum;j++){
            IF>>hubid>>dis;
            Label[nodeID].insert(make_pair(hubid,dis));
        }
    }
    IF.close();
    unordered_map<int,vector<int>> unorderm; unorderm.clear();
//    PruningPointNew.assign(num,unorderm);

//    ifstream IF1(filenameP);
//
//    int pairnum, pairsize;
//    int c,u;
//    for(int i=0;i<num;i++){
//        IF1>>nodeID>>pairnum;
//        for(int p=0;p<pairnum;p++){
//            IF1>>c>>pairsize;
//            vector<int> vec;
//            for(int k=0;k<pairsize;k++){
//                IF1>>u;
//                vec.push_back(u);
//            }
//            PruningPointNew[nodeID][c]=vec;
//        }
//    }
//    IF1.close();
    cout<<"PLL index finish reading!"<<endl;
}

void Graph::IndexSizePLL(){
    unsigned long long m=0,m1=0,m2=0;

    //core index
    for(int k=0;k<Label.size();k++){
        m1+=Label[k].size()*2*sizeof(int);
    }

    for(int i=0;i<PruningPointSet.size();i++){
        for(auto it=PruningPointSet[i].begin();it!=PruningPointSet[i].end();it++){
            m2+=(1+(*it).second.size())*sizeof(int);
        }
    }

    m=m1+m2;
    cout<<"PLL label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"Pruning point size: "<<(double)m2/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}
//function of PLL query
int Graph::QueryPLL(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;

    int hub, dis1, dis2;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                //cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
            }
        }
    }
    return d;
}



/// PSL
// PSL construction entry
void Graph::PSLConstruction(int strategy, bool ifDynamic){
    string indexfile=graph_path+".PSL";
    string PPRfile=graph_path+"PSL.PPR";

    VertexOrdering(strategy);

    if(ifDynamic){
        cout<<"For dynamic graph."<<endl;
    }else{
        cout<<"For static graph."<<endl;
    }


    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;
    bool ifHash=true;
    ifHash= false;

    if(ifHash){
        t1=std::chrono::high_resolution_clock::now();
        PSLConstruct();//hash table-based implementation
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"PSL construction time: "<<runT<<" s."<<endl;

    }else{
        vector<vertex> vertices;
        for(int i=0;i<node_num;++i){
            vertices.emplace_back(i);
        }
        t1=std::chrono::high_resolution_clock::now();
        PSLConstructV(Neighbor);//vector-based implementation
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"PSL construction time: "<<runT<<" s."<<endl;

        LabelV.postProcess(Label,threadnum);
        PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
    }

    IndexSizePLL();
}
// PSL construction function
void Graph::PSLConstruct(){
    cout<<"Hash table-based implementation."<<endl;
    Label.assign(node_num, unordered_map<vertex,int>());
    Dhop.assign(node_num, unordered_map<vertex,int>());

    DvertexNew.assign(node_num, true);

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }

    for(int i=0;i<node_num;i++){
        Label[i].insert(make_pair(i,0));
        Dhop[i].insert(make_pair(i,0));
        //Dvectex.push_back(i);
    }

    if(ifDynamic){
        PruningPointSet.assign(node_num,unordered_map<vertex, unordered_set<vertex>>());
        PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
    }

    //LabelStep.push_back(Dhop);
    Timer tt;
    bool flag=true;
    int step=0;
    cout<<"step "<<step<<" finish!"<<endl;
    step+=1;
    while(flag){
//		LabelStep.push_back(Dhop);
        tt.start();
        flag=DhopLableRefreshStep(step+1);
        tt.stop();
        cout<<"step "<<step<<" finish! "<<tt.GetRuntime()<<" s."<< endl;
        step+=1;
    }
    cout<<"Index finish construction"<<endl;
}
// vector-based label
void Graph::PSLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Vector-based implementation."<<endl;
    unordered_map<int,int> map0; map0.clear();
//    Label.assign(node_num, map0);
//    Dhop.assign(node_num, map0);
    LabelV.resize(node_num);
    PPRV.resize(node_num);
    DhopV.assign(node_num,vector<pair<vertex,int>>());

    DvertexNew.assign(node_num, true);

    #pragma omp parallel for num_threads(threadnum)
    for(int i=0;i<node_num;i++){
//        Label[i].insert(make_pair(i,0));
        LabelV.add(i,i,0);
//        Dhop[i].insert(make_pair(i,0));
        DhopV[i].emplace_back(i,0);
        //Dvectex.push_back(i);
    }

    //LabelStep.push_back(Dhop);
    Timer tt;
    bool flag=true;
    int step=0;
    cout<<"step "<<step<<" finish!"<<endl;
    step+=1;
    int stepShow=1;
    if(edge_num/node_num<3){
        stepShow=100;
    }
    while(flag){
//		LabelStep.push_back(Dhop);
        tt.start();
        flag=DhopLableRefreshStepV(step+1, Neighbor);
        tt.stop();
        if(step%stepShow==0){
            cout<<"step "<<step<<" finish! "<<tt.GetRuntime()<<" s."<< endl;
        }

        step+=1;
    }
    LabelV.sort(threadnum);
    cout<<"Index finish construction"<<endl;
}

bool Graph::DhopLableRefreshStep(int step){
    bool flag=false;
    vector<unordered_map<vertex,int>> newDhop;
    newDhop.assign(node_num,unordered_map<vertex,int>());
    vector<int> newDvec;

    vector<vector<int>> ProcessID;
    vector<int> vvv; ProcessID.assign(threadnum,vvv);
    threadDistribute(ProcessID);

    /// boost-based implementation
    boost::thread_group thread;
    for(int i=0;i<ProcessID.size();i++){
        vector<int> p=ProcessID[i];
        thread.add_thread(new boost::thread(&Graph::labelMultiThread2New, this, boost::ref(newDhop), p,step));
    }
    thread.join_all();
//    for(int i=0;i<ProcessID.size();i++){
//        vector<int> p=ProcessID[i];
//        labelMultiThread2New(newDhop, p,step);
//    }
    /// openmp-based implementation
//    #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//    for(int i=0;i<ProcessID.size();i++){
//        vector<int> p=ProcessID[i];
//        labelMultiThread2New(newDhop, p,step);
//    }

    int zerolabel=0;
    Dhop.assign(newDhop.begin(), newDhop.end());
    set<int> Dset;
    DvertexNew.assign(node_num,false);
    for(int nodeID=0;nodeID<node_num;nodeID++){
        if(Dhop[nodeID].size()>0){
            flag=true;
            for(auto it=Dhop[nodeID].begin();it!=Dhop[nodeID].end();it++){
                if(Label[nodeID].find((*it).first)!=Label[nodeID].end()){//if found
                    Label[nodeID][(*it).first]=(*it).second;
                }else{
                    Label[nodeID].insert(*it);
                }
            }
            for(int i=0;i<Neighbor[nodeID].size();i++){//Neighbors
                if(Dset.find(Neighbor[nodeID][i].first)==Dset.end()){
                    Dset.insert(Neighbor[nodeID][i].first);
                    //Dvectex.push_back(Neighbors[nodeID][i].first);
                    DvertexNew[Neighbor[nodeID][i].first]=true;
                }
            }
        }else if(Dhop[nodeID].size()==0)
            zerolabel+=1;
    }

    //cout<<"zero "<<zerolabel<<" Vertex to change "<<Dvectex.size()<<endl;
    return flag;
}
// vector-based implementation
bool Graph::DhopLableRefreshStepV(int step, vector<vector<pair<vertex,int>>> &Neighbor){
    bool flag=false;

    vector<vector<pair<vertex,int>>> newDhopV;
    newDhopV.assign(node_num,vector<pair<vertex,int>>());

    /// boost-based implementation
    vector<int> ProcessID;
    for(int r=node_num-1;r>=0;r--){
        int ID=vNodeOrder[r];
        ProcessID.push_back(ID);
    }
    vector<vector<int>> ProcessIDs;
    ProcessIDs.assign(threadnum,vector<int>() );
    for(int i=0;i<ProcessID.size();++i){
        int p=i%threadnum;
        ProcessIDs[p].push_back(ProcessID[i]);
    }
    boost::thread_group thread;
    for(int i=0;i<ProcessIDs.size();i++){
        thread.add_thread(new boost::thread(&Graph::labelMultiThread2NewV2, this, boost::ref(newDhopV), boost::ref(ProcessIDs[i]), boost::ref(Neighbor)));
    }
    thread.join_all();

    /// openmp-based implementation
//    vector<int> ProcessID;
//    for(int r=0;r<node_num;r++){
//        int ID=vNodeOrder[r];
//        if(DvertexNew[ID]){
//            ProcessID.push_back(ID);
//        }
//    }
//    #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//    for (int i = 0; i < ProcessID.size(); ++i) {
//        int ID=ProcessID[i];
//        labelMultiThread2NewV(newDhopV,ID,Neighbor);
//    }


    DhopV.assign(newDhopV.begin(), newDhopV.end());
    DvertexNew.assign(node_num,false);

    int hub;
    for(int nodeID=0;nodeID<node_num;nodeID++){
        if(DhopV[nodeID].size()>0){
            flag=true;
            unordered_map<int,int> imap; imap.clear();
            for(int i=0;i<LabelV.Labels[nodeID].size();++i){
                hub=LabelV.Labels[nodeID][i].first;
                imap.insert({hub,i});
            }
            for(auto it=DhopV[nodeID].begin();it!=DhopV[nodeID].end();++it){
                hub=it->first;
                if(imap.find(hub) != imap.end()){//if found
                    LabelV.Labels[nodeID][imap[hub]].second=it->second;
                }else{
                    LabelV.add(nodeID,it->first,it->second);
                }

            }
            for(int i=0;i<Neighbor[nodeID].size();i++){//Neighbors
                int neiV = Neighbor[nodeID][i].first;
                if(neiV<0 || neiV >= node_num){
                    cout<<"Wrong "<<neiV<<endl;
                }

                if(!DvertexNew[neiV]){
                    DvertexNew[neiV]=true;
                }
            }
        }
    }

    //cout<<"zero "<<zerolabel<<" Vertex to change "<<Dvectex.size()<<endl;
    return flag;
}

void Graph::threadDistribute(vector<vector<int>>& processID){
    int ID;
    int cnt=0;
    int threadOrder;
    int a;
    for(int r=0;r<node_num;r++){
        ID=vNodeOrder[r];
        if(DvertexNew[ID]){
            a=cnt%(2*threadnum);
            if(a>=threadnum){
                threadOrder=2*threadnum-1-a;
            }else{
                threadOrder=a;
            }
            //cout<<"a "<<a<<" threadOrder "<<threadOrder<<endl;
            processID[threadOrder].push_back(ID);
            cnt+=1;
        }
    }
}

void Graph::labelMultiThread2New(vector<unordered_map<vertex,int>>& newDhop, vector<int>& p,int step){
//	sm->wait();
    unordered_map<vertex,int> Dhop0;
    for(int i=0;i<p.size();i++){
        int nodeID=p[i];
         Dhop0.clear();
        int neighID, neighW;

        if(ifDynamic){
            for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
                neighID=Neighbor[nodeID][Index].first;
                //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
                if(Dhop[neighID].size()>0){
                    neighW=Neighbor[nodeID][Index].second;

                    int hub, dis, d;
                    for(auto it=Dhop[neighID].begin();it!=Dhop[neighID].end();it++){
                        hub=(*it).first; dis=(*it).second; d=neighW+dis;
                        if(NodeOrder[hub]>NodeOrder[nodeID]){
                            int TempDis; vector<int> SupNode;
                            ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
//                            TempDis = ShortestDisQuery2(nodeID, hub);
                            if(TempDis>d){
                                if(Dhop0.find(hub)!=Dhop0.end()){
                                    if(Dhop0[hub]>d) Dhop0[hub]=d;
                                }else{
                                    Dhop0.insert(make_pair(hub,d));
                                }
                            }
                            else{
                                for(int k=0;k<SupNode.size();k++){
                                    int supn=SupNode[k];

                                    if(supn!=nodeID && supn!=hub){
                                        vSm[nodeID]->wait();
                                        PruningPointSet[nodeID][supn].insert(hub);
                                        PruningPointSet2[nodeID][hub]=supn;
                                        vSm[nodeID]->notify();

                                        vSm[hub]->wait();
                                        PruningPointSet[hub][supn].insert(nodeID);
                                        PruningPointSet2[hub][nodeID]=supn;
                                        vSm[hub]->notify();
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }else{
            for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
                neighID=Neighbor[nodeID][Index].first;
                //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
                if(Dhop[neighID].size()>0){
                    neighW=Neighbor[nodeID][Index].second;

                    int hub, dis, d;
                    for(auto it=Dhop[neighID].begin();it!=Dhop[neighID].end();it++){
                        hub=(*it).first; dis=(*it).second; d=neighW+dis;
                        if(NodeOrder[hub]>NodeOrder[nodeID]){
                            int TempDis; vector<int> SupNode;
//                        ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
                            TempDis = ShortestDisQuery2(nodeID, hub);
                            if(TempDis>d){
                                if(Dhop0.find(hub)!=Dhop0.end()){
                                    if(Dhop0[hub]>d) Dhop0[hub]=d;
                                }else{
                                    Dhop0.insert(make_pair(hub,d));
                                }
                            }

                        }
                    }
                }
            }
        }

        newDhop[nodeID]=Dhop0;
    }

    //cout<<"one thread finish running!"<<endl;

//	sm->notify();
}
//vector-based implementation
void Graph::labelMultiThread2NewV(vector<vector<pair<vertex,int>>>& newDhop, int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
//	sm->wait();
    unordered_map<vertex,int> Lh;
    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
        Lh.insert({it->first,it->second});
    }

    unordered_map<int,int> Dhop0; Dhop0.clear();
    int neighID, neighW;

    if(ifDynamic){
        for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
            neighID=Neighbor[nodeID][Index].first;
            //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
            if(DhopV[neighID].size()>0){//if neighbor is active
                neighW=Neighbor[nodeID][Index].second;

                int hub, dis, d;
                for(auto it=DhopV[neighID].begin();it!=DhopV[neighID].end();++it){
                    hub=(*it).first; dis=(*it).second; d=neighW+dis;
                    if(NodeOrder[hub]>NodeOrder[nodeID]){
                        int TempDis; vector<int> SupNode;
//                        ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
//                    TempDis = ShortestDisQuery2(nodeID, hub);
                        TempDis = PLLDisQuery1V(hub, Lh, SupNode);
                        if(TempDis>d){
                            if(Dhop0.find(hub)!=Dhop0.end()){
                                if(Dhop0[hub]>d) Dhop0[hub]=d;
                            }else{
                                Dhop0.insert(make_pair(hub,d));
                            }
                        }
                        else{
                            for(int k=0;k<SupNode.size();k++){
                                int supn=SupNode[k];

                                if(supn!=nodeID && supn!=hub){
                                    PPRV.add(nodeID,supn,hub);
                                }
                            }
                        }

                    }
                }
            }
        }
    }else{
        for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
            neighID=Neighbor[nodeID][Index].first;
            //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
            if(DhopV[neighID].size()>0){//if neighbor is active
                neighW=Neighbor[nodeID][Index].second;

                int hub, dis, d;
                for(auto it=DhopV[neighID].begin();it!=DhopV[neighID].end();++it){
                    hub=(*it).first; dis=(*it).second; d=neighW+dis;
                    if(NodeOrder[hub]>NodeOrder[nodeID]){
                        int TempDis; vector<int> SupNode;
//                        ShortestDisQuery1(nodeID, hub,SupNode,TempDis); //cout<<"nodeID "<<nodeID<<" TempDis "<<TempDis<<" d "<<d<<endl;
//                    TempDis = ShortestDisQuery2(nodeID, hub);
                        TempDis = PLLDisQuery1V(hub, Lh, SupNode);
                        if(TempDis>d){
                            if(Dhop0.find(hub)!=Dhop0.end()){
                                if(Dhop0[hub]>d) Dhop0[hub]=d;
                            }else{
                                Dhop0.insert(make_pair(hub,d));
                            }
                        }

                    }
                }
            }
        }
    }

    for(auto it=Dhop0.begin();it!=Dhop0.end();++it){
        newDhop[nodeID].emplace_back(it->first,it->second);
    }

}

void Graph::labelMultiThread2NewV2(vector<vector<pair<vertex,int>>>& newDhop, vector<int>& p, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<p.size();++i){
        int ID=p[i];
        labelMultiThread2NewV(newDhop,ID,Neighbor);
    }
}


//query by current labels
int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return d;
}

//query by current labels
int Graph::ShortestDisQuery1V(int ID1,unordered_map<vertex,int>& Lh,vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;

    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();++it){
        hub=(*it).first;
        dis1=(*it).second;
        if(Lh.find(hub)!=Lh.end()){//if found
            dis2=Lh[hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return d;
}

//query by current labels
int Graph::ShortestDisQuery2(int ID1,int ID2){
    int d=INF;
    int hub, dis1, dis2;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
            }
        }
    }

    return d;
}

//function of querying by current label, query peak
int Graph::ShortestDisQueryPeak(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;
    int finalHub=-1;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                finalHub=hub;
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return finalHub;
}

int Graph::ShortestDisQueryPeakV(int ID1,int ID2, unordered_map<vertex,int>& Lh, vector<int>& SupNode, int& d){
    d=INF;
    int hub, dis1, dis2;
    int finalHub=-1;

    for(auto it=LabelV.Labels[ID1].begin();it!=LabelV.Labels[ID1].end();++it){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Lh.find(hub)!=Lh.end()){
            dis2=Lh[hub];
            if(dis1+dis2<d){
                finalHub=hub;
                d=dis1+dis2;
                SupNode.clear();
                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
                SupNode.push_back(hub);
            }
        }
    }

    return finalHub;
}

//// PCL algorithm
// PSL construction entry
void Graph::PCLConstruction(int strategy, bool ifDynamic){
    string indexfile=graph_path+".PCL";
    string PPRfile=graph_path+"PCL.PPR";

    VertexOrdering(strategy);

    if(ifDynamic){
        cout<<"For dynamic graph."<<endl;

    }else{
        cout<<"For static graph."<<endl;
    }

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;
    bool ifHash=true;
    ifHash= false;
    if(ifHash){
//        cout<<"Hash table-based implementation."<<endl;
        t1=std::chrono::high_resolution_clock::now();
        runT = BPCLConstruct(Neighbor);//hash table-based implementation
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        cout<<"PCL construction time: "<<runT<<" s. Overall time: "<<time_span.count()<< endl;

    }else{
//        cout<<"Vector-based implementation."<<endl;
        t1=std::chrono::high_resolution_clock::now();
        runT = BPCLConstructV(Neighbor);//vector-based implementation
//        runT = PCLConstructV(Neighbor);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        cout<<"PCL construction time: "<<runT<<" s. Overall time: "<<time_span.count()<< endl;
    }

    IndexSizePLL();
}
//function for index construction of batch PLL
void Graph::PCLConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;

    Label.assign(node_num,unordered_map<vertex,int>());

    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

//    PruningPointVector.clear();
//    PruningPointVector.assign(node_num,vector<pair<int,int>>());


    int ID;
    vector<int> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        vertices.emplace_back(ID);
    }

    vector<vector<int>> processID;
    ThreadDistribute(vertices, processID);

    cout<<"Batch size: "<<processID[0].size()<<endl;

    tt.start();
    if(ifParallel){//use multiple thread
//        if(false){
        boost::thread_group thread;
        for(int i=0;i<processID.size();i++){
            thread.add_thread(new boost::thread(&Graph::PCLDijks, this, processID[i], boost::ref(Neighbor)));

        }
        thread.join_all();
    }else{//use single thread
        for(int i=0;i<processID.size();i++){
            PCLDijks(processID[i], Neighbor);
        }
    }

    tt.stop();
    cout<<"Time used for label construction: "<<tt.GetRuntime()<<" s."<<endl;

    tt.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel,processID, Neighbor);
//    PruningPointBuild(false);
    tt.stop();
    cout<<"Time used for pruning point construction: "<<tt.GetRuntime()<<" s."<<endl;
}
//function for index construction of batch PLL, vector-based implementation
double Graph::PCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Vector-based implementation of PCL."<<endl;
    bool ifParallel = true;
    double t1,t2;

    LabelV.resize(node_num);


    int ID;
    vector<int> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        vertices.emplace_back(ID);
    }

    vector<vector<int>> processID;
    ThreadDistribute(vertices, processID);

    cout<<"Batch size: "<<processID[0].size()<<endl;

    tt.start();
    #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
    for(int i=0;i<vertices.size();++i){
        int ID=vertices[i];
        PCLDijkV(ID,Neighbor);
    }
    LabelV.sort(threadnum);
    tt.stop();
    t1=tt.GetRuntime();
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;

//    cout<<"Begin to construct pruning points..."<<endl;
//    tt.start();
//    PruningPointBuild(ifParallel,processID, Neighbor);
//    tt.stop();
//    cout<<"Time used for pruning point construction: "<<tt.GetRuntime()<<" s."<<endl;

/// openmp-based based implementation, post-processing
    cout<<"Begin to construct pruning points... OpenMP-based implementation with post-processing."<<endl;
    PPRV.resize(node_num);
    tt.start();
#pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
    for(int id=0;id<node_num;++id){
        PPRConstructionV(id, Neighbor);
    }
    tt.stop();
    t2=tt.GetRuntime();
    cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;

    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
    PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);

    return t1+t2;
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based
void Graph::PCLDijks(vector<int> & IDs, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<IDs.size();++i){
        PCLDijk(IDs[i], Neighbor);
    }
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based
void Graph::PCLDijk(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path
    vector<bool> pruned(node_num,false);//record the PPR pruning information

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=12516, hid=11449;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        if(pruned[pre[topNodeID]]){
            pruned[topNodeID]= true;
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        //Query of dispeak
        if(NodeOrder[nodeID] <= NodeOrder[topNodeID] && NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]){//if valley path
//            vSm[nodeID]->wait();
            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
//            vSm[nodeID]->notify();
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;

            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario

                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
//            if(NNodeID==lid){
//                cout<<topNodeID<< " "<<NNodeID<<" "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
//            }
        }
    }
}
//function of original Dijkstra of PCL, 2023-06-16, pure search-based, vector-based implementation
void Graph::PCLDijkV(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path
    vector<bool> pruned(node_num,false);//record the PPR pruning information

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=12516, hid=11449;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        if(pruned[pre[topNodeID]]){
            pruned[topNodeID]= true;
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        //Query of dispeak
        if(NodeOrder[nodeID] <= NodeOrder[topNodeID] && NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]){//if valley path
//            vSm[nodeID]->wait();
//            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
            LabelV.add_lockfree(nodeID,topNodeID,topNodeDis);
//            vSm[nodeID]->notify();
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;

            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario

                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
//            if(NNodeID==lid){
//                cout<<topNodeID<< " "<<NNodeID<<" "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
//            }
        }
    }
}
// hash table-based implementation
double Graph::BPCLConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    bool ifParallel = true;
    double t1,t2;
    cout<<"Hash table-based implementation."<<endl;
    
    Label.assign(node_num,unordered_map<vertex,int>());
    PruningPointSet2.clear();
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());
//    PruningPointSetOrder.clear();
//    PruningPointSetOrder.assign(node_num,unordered_map<int,set<OrderComp>>());
    PruningPointSet.clear();
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*10000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = batchsize;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<int>> batches;
    vector<int> bNodes;//vertices for current batch
    vector<vector<int>> ProcessID;


    vector<int> vertices;

    Timer tt;
    double time=0;
    int a=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];
        
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
            }
        }
    }

    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
        hID = bNodes[0];
        unordered_set<int> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();

        /// boost-based implementation
        if(ifParallel){//use multiple thread
//        if(false){
            if(batchSize > threadnum){
                ProcessID.assign(threadnum,vector<int>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%threadnum;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }else if(batchSize == threadnum){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijk(ID, setNodes, hID, Neighbor);
            }
        }

        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    t1=tt2.GetRuntime();
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;

    ///boost-based implementation
    cout<<"Begin to construct pruning points... Boost-based implementation."<<endl;
    vector<vector<int>> processID;
    processID.assign(threadnum, vector<int>());
    ThreadDistribute(vertices, processID);

    tt2.start();
    PruningPointBuild(ifParallel, processID, Neighbor);
    tt2.stop();
    t2=tt2.GetRuntime();
    cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;

    return t1+t2;
}
// vector-based implementation
double Graph::BPCLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    double t1,t2;

    cout<<"Vector-based implementation of BPCL."<<endl;
    bool ifParallel = true;

    Label.assign(node_num,unordered_map<vertex,int>());
    LabelV.resize(node_num);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*10000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = batchsize;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<int>> batches;
    vector<int> bNodes;//vertices for current batch
    vector<vector<int>> ProcessID;


    vector<int> vertices;

    Timer tt;
    double time=0;
    int a=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        vertices.emplace_back(ID);
//        vertices.emplace_back(i);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
            }
        }
    }

    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
        hID = bNodes[0];
        unordered_set<int> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<pair<int,int>> vp;
        // process each batch
        tt.start();


        /// boost-based implementation
        if(ifParallel){//use multiple thread
//        if(false){
            if(batchSize > threadnum){
                /// boost-based implementation
                ProcessID.assign(threadnum,vector<int>());
                for(int j=0;j<bNodes.size();++j){
                    a = j%threadnum;
                    ProcessID[a].emplace_back(bNodes[j]);
                }
                boost::thread_group thread;
                for(int i=0;i<ProcessID.size();i++){
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijk2V, this, boost::ref(ProcessID[i]), boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
                /// openmp-based implementation
//                #pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//                for(int j=0;j<bNodes.size();++j){
//                    int id=bNodes[j];
//                    BatchPCLDijkV(id, setNodes, hID, Neighbor);
//                }

            }else if(batchSize == threadnum){
                boost::thread_group thread;
                for(int i=0;i<bNodes.size();i++){
                    ID = bNodes[i];
                    thread.add_thread(new boost::thread(&Graph::BatchPCLDijkV, this, ID, boost::ref(setNodes), hID, boost::ref(Neighbor)));
                }
                thread.join_all();
            }

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
                BatchPCLDijk(ID, setNodes, hID, Neighbor);
            }
        }

        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    t1=tt2.GetRuntime();
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;
    LabelV.sort(threadnum);

    if(ifDynamic){
        /// boost-based based implementation, post-processing
        cout<<"Begin to construct pruning points... boost-based implementation with post-processing."<<endl;
        PPRV.resize(node_num);
        tt2.start();
        ThreadDistribute(vertices,ProcessID);
        boost::thread_group thread;
        for(int i=0;i<ProcessID.size();i++){
            thread.add_thread(new boost::thread(&Graph::PPRConstructionV2, this, boost::ref(ProcessID[i]),  boost::ref(Neighbor)));
        }
        thread.join_all();
        tt2.stop();
        t2=tt2.GetRuntime();
        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;
        cout<<"Core index construction time (labeling + PPR): "<<t1+t2<<" s."<<endl;


        PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;
    }


    return t1+t2;
}
//function of pruning Dijkstra of batch PLL, 2023-06-07, optimal correct version 1
void Graph::BatchPCLDijk(int nodeID, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=200291, hid=200093;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[hID]) {
            continue;
        }
        else{//if NodeOrder[topNodeID] <= NodeOrder[hID]
            //Query of dispeak
            if(NodeOrder[nodeID] >= NodeOrder[topNodeID]){
                vSm[nodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->wait();
                }
                ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);
                vSm[nodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->notify();
                }


                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {

                        vSm[topNodeID]->wait();
                        if (Label[topNodeID].find(nodeID) == Label[topNodeID].end()) {//if not found
                            Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
                        }
                        vSm[topNodeID]->notify();
                    }

                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){
                        continue;
                    }


                }
            }
            else{//if NodeOrder[nodeID] < NodeOrder[topNodeID]
                vSm[topNodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[nodeID]->wait();
                }
                ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);
                vSm[topNodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[nodeID]->notify();
                }

                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {//if the topNodeID is the highest vertex in the path, reversely insert it.

                        bool update=false;
                        vSm[nodeID]->wait();
                        if (Label[nodeID].find(topNodeID) == Label[nodeID].end()) {//if not found
                            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
                            update=true;
                        }else if(Label[nodeID][topNodeID] > topNodeDis){
                            cout<<"Larger label!! "<<Label[nodeID][topNodeID]<<" "<<topNodeDis<<endl;
                            Label[nodeID][topNodeID] = topNodeDis;
                            update=true;
                        }
                        vSm[nodeID]->notify();

                    }
                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){//>
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){//if the dispeak is sourced from higher vertex
                        continue;
                    }

                }
            }
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    if(ifDebug){
                        if(nodeID == hid && NNodeID == lid){
                            cout<<"here! "<<nodeID<<" "<<topNodeID<<"("<<peak[topNodeID]<<","<<NodeOrder[peak[topNodeID]]<<"), "<<NNodeID<<"("<<peak[NNodeID]<<","<<NodeOrder[peak[NNodeID]]<<"): "<<distance[NNodeID]<<" "<<NWeigh<<" "<<endl;
                        }
                    }

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario
                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
        }
    }
}
void Graph::BatchPCLDijkV(int nodeID, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

//    unordered_map<vertex,uint> Lh;
//    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
//        Lh.insert({it->first,it->second});
//    }

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> peak(node_num,-1);//record the peak hub from nodeID to topNodeID
    vector<int> pre(node_num,-1);//record the predecessor vertex in the shortest path

    pre[nodeID]=nodeID;
    peak[nodeID]=nodeID;
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=181451, hid=190161;
    lid=200291, hid=200093;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        if(NodeOrder[peak[topNodeID]] < NodeOrder[peak[pre[topNodeID]]]) {
            peak[topNodeID] = peak[pre[topNodeID]];
        }
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[hID]) {
            continue;
        }
        else{//if NodeOrder[topNodeID] <= NodeOrder[hID]
            //Query of dispeak
            if(NodeOrder[nodeID] >= NodeOrder[topNodeID]){
                vSm[nodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->wait();
                }
                ShortestDisQuery1V(topNodeID, Label[nodeID], SupNode,TempDis);
                vSm[nodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[topNodeID]->notify();
                }


                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {

                        vSm[topNodeID]->wait();
                        if (Label[topNodeID].find(nodeID) == Label[topNodeID].end()) {//if not found
                            Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
                            LabelV.add(topNodeID,nodeID,topNodeDis);
                        }

                        vSm[topNodeID]->notify();
                    }

                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){
                        continue;
                    }


                }
            }
            else{//if NodeOrder[nodeID] < NodeOrder[topNodeID]
                vSm[topNodeID]->wait();
                if(topNodeID != nodeID){
                    vSm[nodeID]->wait();
                }
                ShortestDisQuery1V(topNodeID, Label[nodeID],SupNode,TempDis);
                vSm[topNodeID]->notify();
                if(topNodeID != nodeID){
                    vSm[nodeID]->notify();
                }

                if(TempDis>topNodeDis) {//if the peak path is longer than this vally path
                    if((NodeOrder[peak[topNodeID]] <= NodeOrder[topNodeID]) || (NodeOrder[peak[topNodeID]] <= NodeOrder[nodeID])) {//if the topNodeID is the highest vertex in the path, reversely insert it.

                        bool update=false;
                        vSm[nodeID]->wait();
                        if (Label[nodeID].find(topNodeID) == Label[nodeID].end()) {//if not found
                            Label[nodeID].insert(make_pair(topNodeID, topNodeDis));
                            LabelV.add(nodeID,topNodeID,topNodeDis);
                            update=true;
                        }else if(Label[nodeID][topNodeID] > topNodeDis){
                            cout<<"Larger label!! "<<Label[nodeID][topNodeID]<<" "<<topNodeDis<<endl;
                            Label[nodeID][topNodeID] = topNodeDis;
                            update=true;
                            exit(1);
                        }
                        vSm[nodeID]->notify();

                    }
                }
                else{//if TempDis<=topNodeDis
                    int highestHub = -1;
                    for (int k = 0; k < SupNode.size(); k++) {
                        int supn = SupNode[k];
                        if(highestHub == -1){
                            highestHub = supn;
                        }else if(NodeOrder[supn]>NodeOrder[highestHub]){//>
                            highestHub = supn;
                        }

                        ifPrune = true;
                        //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
                    }

                    if(NodeOrder[highestHub] > NodeOrder[hID]){//if the dispeak is sourced from higher vertex
                        continue;
                    }

                }
            }
        }

        //Label[topNodeID].insert(nodeID, topNodeDis);
        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    if(ifDebug){
                        if(nodeID == hid && NNodeID == lid){
                            cout<<"here! "<<nodeID<<" "<<topNodeID<<"("<<peak[topNodeID]<<","<<NodeOrder[peak[topNodeID]]<<"), "<<NNodeID<<"("<<peak[NNodeID]<<","<<NodeOrder[peak[NNodeID]]<<"): "<<distance[NNodeID]<<" "<<NWeigh<<" "<<endl;
                        }
                    }

                    pre[NNodeID]=topNodeID;
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);

                    if(peak[NNodeID] == -1){
                        peak[NNodeID] = NNodeID;
                    }

                }else if(distance[NNodeID]==NWeigh){//deal with equal distance scenario
                    if(NodeOrder[peak[topNodeID]] > NodeOrder[peak[pre[NNodeID]]]){//record the highest vertex among all shortest paths
//                        cout<<"change pre. L("<<nodeID<<","<<NNodeID<<"): "<<topNodeID<<"("<<NodeOrder[topNodeID]<<") "<<pre[NNodeID]<<"("<<NodeOrder[pre[NNodeID]]<<")"<<endl;
                        pre[NNodeID]=topNodeID;
                    }

                }
            }
        }
    }
}
//function of pruning Dijkstra of batch PLL, 2023-06-07, optimal correct version 1
void Graph::BatchPCLDijk2(vector<int>& p, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(auto it=p.begin();it!=p.end();++it){
        int nodeID = *it;
        BatchPCLDijk(nodeID, setNodes, hID, Neighbor);
    }
}
//vector-based implementation
void Graph::BatchPCLDijk2V(vector<int>& p, unordered_set<int>& setNodes, int hID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(auto it=p.begin();it!=p.end();++it){
        int nodeID = *it;
        BatchPCLDijkV(nodeID, setNodes, hID, Neighbor);
    }
}

//function of vertex allocation
void Graph::ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID){
    processID.assign(threadnum, vector<int>());
    int pid=0;
    for(int i=0;i<vertices.size();++i){
        pid=i%threadnum;
        processID[pid].emplace_back(vertices[i]);
    }
}
//function of building the pruning points
void Graph::PruningPointBuild(bool ifParallel, vector<vector<int>> & processID, vector<vector<pair<vertex,int>>> &Neighbor){

    if(ifParallel){
//    if(false){
        cout<<"Batch size: "<<processID[0].size()<<endl;
        boost::thread_group thread;
        for(auto j=0;j<processID.size();++j){
            thread.add_thread(new boost::thread(&Graph::PPRConstruction2, this, boost::ref(processID[j]), boost::ref(Neighbor) ));
        }
        thread.join_all();
    }else{
        for(auto j=0;j<processID.size();++j){
            PPRConstruction2(processID[j], Neighbor);
        }
    }
}
void Graph::PPRConstruction2(vector<int> & p, vector<vector<pair<vertex,int>>> &Neighbor){
//    cout<<"1"<<endl;
    for(int i=0;i<p.size();++i){
        PPRConstruction(p[i], Neighbor);
    }
}
//function of building the pruning point records for a given vertex
void Graph::PPRConstruction(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    pqueue.update(nodeID,0);
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        int TempDis; vector<int> SupNode;
        //if NodeOrder[topNodeID] <= NodeOrder[nodeID]
        ShortestDisQueryPeak(topNodeID, nodeID, SupNode,TempDis);//DisQueryPeak
//        ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);

        if(TempDis<=topNodeDis){//if dispeak <= disvally
            for(int k=0;k<SupNode.size();k++){
                int supn=SupNode[k];

                vSm[nodeID]->wait();
//                PruningPointVector[nodeID].emplace_back(supn,topNodeID);
//                PruningPointSetOrder[nodeID][supn].insert(topNodeID);
                PruningPointSet[nodeID][supn].insert(topNodeID);
                PruningPointSet2[nodeID][topNodeID] = supn;
                vSm[nodeID]->notify();

                vSm[topNodeID]->wait();
//                PruningPointVector[topNodeID].emplace_back(supn,nodeID);
//                PruningPointSetOrder[topNodeID][supn].insert(nodeID);
                PruningPointSet[topNodeID][supn].insert(nodeID);
                PruningPointSet2[topNodeID][nodeID] = supn;
                vSm[topNodeID]->notify();

                //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
            }
            continue;

        }

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
}
//vector-based implementation
void Graph::PPRConstructionV(int nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    pqueue.update(nodeID,0);
    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    unordered_map<vertex,int> Lh;
    for(int i=0;i<LabelV.Labels[nodeID].size();++i){
        Lh.insert({LabelV.Labels[nodeID][i].first,LabelV.Labels[nodeID][i].second});
    }

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        int TempDis; vector<int> SupNode;
        //if NodeOrder[topNodeID] <= NodeOrder[nodeID]
        ShortestDisQueryPeakV(topNodeID, nodeID, Lh, SupNode,TempDis);//DisQueryPeak
//        ShortestDisQueryPeak(topNodeID, nodeID, SupNode,TempDis);//DisQueryPeak

        if(TempDis<=topNodeDis){//if dispeak <= disvally
            for(int k=0;k<SupNode.size();k++){
                int supn=SupNode[k];
                unordered_map<int,vector<int>> map1;

                PPRV.add_lockfree(nodeID,supn,topNodeID);

                //cout<<"pruning point "<<topNodeID<<" "<<nodeID<<" "<<supn<<endl;
            }
            continue;

        }

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
}
void Graph::PPRConstructionV2(vector<vertex>& ProcessID, vector<vector<pair<vertex,int>>> &Neighbor){
    for(int i=0;i<ProcessID.size();++i){
        int id=ProcessID[i];
        PPRConstructionV(id,Neighbor);
    }
}

/// GLL algorithm
// GLL construction entry
void Graph::GLLConstruction(int strategy, bool ifDynamic){
    string indexfile=graph_path+".GLL";
    string PPRfile=graph_path+"GLL.PPR";

    VertexOrdering(strategy);

    if(ifDynamic){
        cout<<"For dynamic graph."<<endl;

    }else{
        cout<<"For static graph."<<endl;
    }

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;
    bool ifHash=true;
    ifHash= false;
    if(ifHash){
        t1=std::chrono::high_resolution_clock::now();
        runT = GLLConstruct(Neighbor);//hash table-based implementation
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        cout<<"GLL construction time: "<<runT<<" s. Overall time: "<<time_span.count()<< endl;

    }else{
        t1=std::chrono::high_resolution_clock::now();
        runT = GLLConstructV(Neighbor);//vector-based implementation
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        cout<<"GLL construction time: "<<runT<<" s. Overall time: "<<time_span.count()<< endl;
    }

    IndexSizePLL();
}
//GLL algorithm
double Graph::GLLConstruct(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Hash table-based implementation."<<endl;
    bool ifParallel = true;
    double t1,t2;

    Label.assign(node_num,unordered_map<vertex,int>());
    PruningPointSet.assign(node_num,unordered_map<vertex,unordered_set<vertex>>());
    PruningPointSet2.assign(node_num,unordered_map<vertex,vertex>());

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*1000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = threadnum;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<int>> batches;
    vector<int> bNodes;

    vector<int> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];
        
        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
//                bNodes.emplace_back(ID);
            }
        }
    }
    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
//        hID = bNodes[0];
        unordered_set<int> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<vector<pair<vertex,int>>> localTable;
        localTable.assign(bNodes.size(),vector<pair<vertex,int>>());

        // process each batch
        tt.start();
        if(ifParallel){//use multiple thread
//        if(false){


            boost::thread_group thread;
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                thread.add_thread(new boost::thread(&Graph::DijksPrune2, this, ID, boost::ref(vp)));
                thread.add_thread(new boost::thread(&Graph::DijksPrune4, this, ID,  boost::ref(Neighbor)));
            }
            thread.join_all();


        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                DijksPrune2(ID);
                DijksPrune4(ID, Neighbor);
            }

        }

        BatchClean(bNodes, setNodes, ifParallel);
//        BatchClean(bNodes, setNodes, false);

        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }
    tt2.stop();
    t1=tt2.GetRuntime();
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;

    //original version
    vector<vector<int>> processID;
    processID.assign(threadnum, vector<int>());
    ThreadDistribute(vertices, processID);
    tt2.start();
    cout<<"Begin to construct pruning points..."<<endl;
    PruningPointBuild(ifParallel, processID, Neighbor);
    tt2.stop();
    t2=tt2.GetRuntime();
    cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;

    return t1+t2;
}

double Graph::GLLConstructV(vector<vector<pair<vertex,int>>> &Neighbor){
    cout<<"Vector-based implementation."<<endl;
    bool ifParallel = true;
    double t1=0,t2=0;

    Label.assign(node_num,unordered_map<vertex,int>());
    LabelV.resize(node_num);

    int ID;
    int cnt=0;
    int stepShow = ceil(node_num/100000)*1000;
    stepShow = max(stepShow,1000);

    cout<<"Step for show: "<<stepShow<<endl;

    int batchSize = 512;
    batchSize = threadnum;
    if(batchSize>node_num){
        batchSize=node_num;
    }
//    batchSize = 3;
    cout<<"Batch size: "<<batchSize<<endl;

    vector<vector<int>> batches;
    vector<int> bNodes;

    vector<int> vertices;

    Timer tt;
    double time=0;

    for(int i=node_num-1;i>=0;i--){
        ID=vNodeOrder[i];

        vertices.emplace_back(ID);
        if(bNodes.size()<batchSize){
            bNodes.emplace_back(ID);
            if(bNodes.size()==batchSize || i==0){
                batches.push_back(bNodes);
                bNodes.clear();
//                bNodes.emplace_back(ID);
            }
        }
    }
    if(!bNodes.empty()){
        batches.push_back(bNodes);
        bNodes.clear();
    }

    int batchNum = batches.size();
    int batchShow = 1;
    if (batchShow < batchNum / 10){
        batchShow = batchNum / 10;
    }

    cout<<"Total batch number: "<<batchNum<<endl;
    cout<<"Batch gap for show: "<<batchShow<<endl;
    int lid, hid;
    lid=200291, hid=200093;

    Timer tt2;
    tt2.start();
    int hID = vNodeOrder[node_num-1];
    for(int b1=0;b1<batches.size();++b1){
        bNodes=batches[b1];
//        hID = bNodes[0];
        unordered_set<int> setNodes; setNodes.clear();
        setNodes.insert(bNodes.begin(),bNodes.end());

        vector<vector<pair<vertex,int>>> localTable;
        localTable.assign(bNodes.size(),vector<pair<vertex,int>>());

        // process each batch
        tt.start();
        if(ifParallel){//use multiple thread
//        if(false){

            boost::thread_group thread;
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                thread.add_thread(new boost::thread(&Graph::DijksPrune2, this, ID, boost::ref(vp)));
                thread.add_thread(new boost::thread(&Graph::DijksPrune4V, this, ID, boost::ref(Neighbor)));
            }
            thread.join_all();

        }else{//use single thread
            for(int i=0;i<bNodes.size();i++){
                ID = bNodes[i];
//                DijksPrune2(ID);
                DijksPrune4V(ID,  Neighbor);
//                CommitLabel(ID,localTable[i]);
            }

        }

        BatchClean(bNodes, setNodes, ifParallel);
//        BatchClean(bNodes, setNodes, false);


        tt.stop();
        time+=tt.GetRuntime();
        if(b1%batchShow==0){
//        if(b1%batchShow==0 && b1<451){
            cout<<"Batch "<<b1<<" finish! Order("<<NodeOrder[bNodes[0]]<<","<<NodeOrder[bNodes[bNodes.size()-1]]<<"). "<<time<<" s."<<endl;
            time = 0;

        }

    }

    tt2.stop();
    t1=tt2.GetRuntime();
    cout<<"Time used for label construction: "<<t1<<" s."<<endl;
    LabelV.sort(threadnum);
    //original version

//    ThreadDistribute(vertices, processID);
//    tt2.start();
//    cout<<"Begin to construct pruning points..."<<endl;
//    PruningPointBuild(ifParallel, processID, Neighbor);
//    tt2.stop();
//    t2=tt2.GetRuntime();
//    cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;

    if(ifDynamic){
        /// boost-based based implementation, post-processing
        cout<<"Begin to construct pruning points... boost-based implementation with post-processing."<<endl;
        PPRV.resize(node_num);
        vector<vector<int>> ProcessID;
        ProcessID.assign(threadnum, vector<int>());
        tt2.start();
        ThreadDistribute(vertices,ProcessID);
        boost::thread_group thread;
        for(int i=0;i<ProcessID.size();i++){
            thread.add_thread(new boost::thread(&Graph::PPRConstructionV2, this, boost::ref(ProcessID[i]),  boost::ref(Neighbor)));
        }
        thread.join_all();
        tt2.stop();
        t2=tt2.GetRuntime();
        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;
        cout<<"Core index construction time (labeling + PPR): "<<t1+t2<<" s."<<endl;


        PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;
//        /// openmp-based based implementation, post-processing
//        cout<<"Begin to construct pruning points... OpenMP-based implementation with post-processing."<<endl;
//        PPRV.resize(node_num);
//        tt2.start();
//#pragma omp parallel for num_threads(threadnum) schedule (dynamic, threadnum)
//        for(int id=0;id<node_num;++id){
//            PPRConstructionV(id, Neighbor);
//        }
//        tt2.stop();
//        t2=tt2.GetRuntime();
//        cout<<"Time used for pruning point construction: "<<t2<<" s."<<endl;
//
//        PPRV.postProcess(PruningPointSet,PruningPointSet2,vertices,threadnum);
    }

    return t1+t2;
}
//function of pruning Dijkstra of LCC
void Graph::DijksPrune4(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }

        vSm[nodeID]->wait();
        if(topNodeID != nodeID){
            vSm[topNodeID]->wait();
        }
        ShortestDisQuery1(topNodeID, nodeID, SupNode,TempDis);
        vSm[nodeID]->notify();
        if(topNodeID != nodeID){
            vSm[topNodeID]->notify();
        }

        if(TempDis<=topNodeDis) {//if the peak path is shorter than this vally path
            continue;
        }

        vSm[topNodeID]->wait();
        Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
        vSm[topNodeID]->notify();


        //Label[topNodeID].insert(nodeID, topNodeDis);
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
}
//vector-based implementation
void Graph::DijksPrune4V(vertex nodeID, vector<vector<pair<vertex,int>>> &Neighbor){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

//    unordered_map<vertex,int> Lh;
//    for(auto it=LabelV.Labels[nodeID].begin();it!=LabelV.Labels[nodeID].end();++it){
//        Lh.insert({it->first,it->second});
//    }

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;
    bool ifPrune;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;
        ifPrune = false;
        int TempDis; vector<int> SupNode;

        if(NodeOrder[topNodeID] > NodeOrder[nodeID]){
            continue;
        }


        vSm[nodeID]->wait();
        if(topNodeID != nodeID){
            vSm[topNodeID]->wait();
        }

        ShortestDisQuery1V(topNodeID, Label[nodeID], SupNode,TempDis);

        vSm[nodeID]->notify();
        if(topNodeID != nodeID){
            vSm[topNodeID]->notify();
        }


        if(TempDis<=topNodeDis) {//if the peak path is shorter than this vally path
            continue;
        }

        vSm[topNodeID]->wait();
        Label[topNodeID].insert(make_pair(nodeID, topNodeDis));
        LabelV.add(topNodeID,nodeID,topNodeDis);
        vSm[topNodeID]->notify();

//        vp.emplace_back(make_pair(topNodeID, topNodeDis));


        //Label[topNodeID].insert(nodeID, topNodeDis);
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
}


//function of cleaning redundant labels generated in current batch
void Graph::BatchClean(vector<int> &bNodes, unordered_set<int>& setNodes, bool ifParallel){
    int ID;
//    Timer tt;
//    tt.start();
    if(ifParallel){
        boost::thread_group thread;
        for(int i=0;i<bNodes.size();i++){
            ID = bNodes[i];
            thread.add_thread(new boost::thread(&Graph::LabelClean, this, ID, boost::ref(setNodes)));
        }
        thread.join_all();
    }else{
        for(auto it=bNodes.begin();it!=bNodes.end();++it){
            ID = *it;
            LabelClean(ID, setNodes);
        }
    }
//    tt.stop();
//    cout<<"Time used for label cleaning: "<<tt.GetRuntime()<<" s."<<endl;
}

//function of clean redundant labels of ID1
void Graph::LabelClean(int ID1, unordered_set<int>& setNodes){
    pair<int,vertex> peakPair;
    int dispeak, Dis;
    vertex ID2, peakhub;
    vector<pair<int,int>> eraseLabels;

    vSm[ID1]->wait();
    for(auto it2=Label[ID1].begin();it2!=Label[ID1].end();++it2){
        ID2=it2->first, Dis=it2->second;
        if(ID1 == ID2) continue;
        if(setNodes.find(ID2) != setNodes.end()){//if found
            vSm[ID2]->wait();
            peakPair=DisQueryPeak2(ID1,ID2);
            vSm[ID2]->notify();
            dispeak=peakPair.first, peakhub=peakPair.second;
            if(dispeak<=Dis){
                eraseLabels.emplace_back(ID1,ID2);
            }
        }else{
            peakPair=DisQueryPeak2(ID1,ID2);
            dispeak=peakPair.first, peakhub=peakPair.second;
            if(dispeak<=Dis){
                eraseLabels.emplace_back(ID1,ID2);
            }
        }
    }
    for(auto it2=eraseLabels.begin();it2!=eraseLabels.end();++it2){
        ID1=it2->first, ID2=it2->second;
        Label[ID1].erase(ID2);
    }
    vSm[ID1]->notify();
}


/*-------------Index-Maintenance Algorithms--------------*/
/// PLL
void Graph::IndexMaintenancePLL( int updateType, int batchNum){
    cout<<"Update Batch: "<<batchNum<<endl;
    vector<pair<pair<int,int>,int>> testdata;
    string updateFile=graph_path+".update";

    ReadUpdate(updateFile, testdata);

    int ID1,ID2, oldW,newW;
    double runT1=0;
    Timer tt1;
    bool ifDebug=false;
//    ifDebug=true;

    switch (updateType) {
        case 0:{
            cout<<"Update type: Decrease"<<endl;
            Graph g1=*this;
            runT1 = 0;
            for(int u=0;u<batchNum;u++){
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*0.5;
                if(newW < 1) {
                    cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    exit(1);
                }
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                }

                tt1.start();
                g1.DecreasePSL(ID1,ID2,oldW,newW,g1.Neighbor);
                tt1.stop();

                runT1 += tt1.GetRuntime();
                if(ifDebug){
                    g1.CorrectnessCheck(100);
                }

//                testdata[u].second = newW;
            }
            cout<<"Average Decrease update Time: "<<runT1/batchNum<<" s.\n"<<endl;
//            break;
        }
        case 1:{
            cout<<"Update type: Increase"<<endl;
            runT1 = 0;
            for(int u=0;u<batchNum;u++){
//            for(int u=406;u<updateBatch;u++){//core 30
                ID1=testdata[u].first.first;
                ID2=testdata[u].first.second;
                oldW=testdata[u].second;
                newW=oldW*1.5;
//                ID1=232897; ID2=232898; oldW=3113; newW=4650;
                if(ifDebug){
                    cout<<"Batch "<<u<<": "<<ID1<<"("<<NodeOrder[ID1]<<") "<<ID2<<"("<<NodeOrder[ID2]<<") "<<oldW<<" "<<newW<<endl;
                }

                tt1.start();
                IncreasePSL(ID1,ID2,oldW,newW,Neighbor);
                tt1.stop();

                runT1 += tt1.GetRuntime();
//                cout<<"After update: "<<g.Query(261222,146486)<<" "<<g.Dijkstra(261222,146486,g.Neighbor)<<endl;
//                if(g.AdjaCoreMap[261222].find(146486) != g.AdjaCoreMap[261222].end())
//                    cout<< g.AdjaCoreMap[261222][146486]<<" "<<g.AdjaCoreMap[146486][261222]<<endl;
                /*if(g.Tree[g.rank[197979]].disInf.find(144762) != g.Tree[g.rank[197979]].disInf.end()){
                    cout<<g.Tree[g.rank[197979]].disInf[144762]<<" ";
                }
                for(auto it=g.Tree[g.rank[197979]].vert.begin();it!=g.Tree[g.rank[197979]].vert.end();++it){
                    if(it->first == 144762){
                        cout<<it->second.first<<" "<<it->second.second;
                        break;
                    }
                }
                cout << endl;*/

//                cout<<"After: "<<g.Query(id1,id2)<<" "<<g.DijkstraCore(id1,id2)<<" "<<g.Dijkstra(id1,id2,g.Neighbor)<<endl;

                if(u==43){//core 15
//                    g.WriteCoreGraph(graphfile+"Core153");
//                    g.WriteLabels(graphfile+"Core153");
                }

//                g.DijkstraPath(63055,233384);
//                g.QueryCoreDebug(142488,143850);
                if(ifDebug){
                    CorrectnessCheck(100);
                }

            }
            cout<<"Average Increase update Time: "<<runT1/batchNum<<" s.\n"<<endl;
            break;
        }
        default:
            break;
    }
}

//function for maintaining the edge weight decrease updates of 2-hop labeling
void Graph::DecreasePSL(vertex a, vertex b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors){
    //update the edges on original graph
    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    //check the dis(a,b)
    int Dab=QueryPLL(a,b);//query by PSL label

    int LID,HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    if(Dab>newW){//the index change is triggered
        vector<vector<pair<int,int>>> Change;
        vector<pair<int,int>> vec;
        Change.assign(node_num,vec);
        set<int> WaitPro;

        Label[LID][HID]=newW;
        Change[LID].push_back(make_pair(HID,newW));
        WaitPro.insert(LID);

        //check the label of a,b
        int hubid, hubdis;

        for(auto it1=Label[LID].begin();it1!=Label[LID].end();++it1){
            hubid=(*it1).first; hubdis=(*it1).second;
            if(NodeOrder[hubid]>NodeOrder[HID] && newW+hubdis<QueryPLL(HID,hubid)){
                Label[HID][hubid]=newW+hubdis;
                Change[HID].push_back(make_pair(hubid, newW+hubdis));
                WaitPro.insert(HID);
            }
        }

        for(auto it2=Label[HID].begin();it2!=Label[HID].end();it2++){
            hubid=(*it2).first; hubdis=(*it2).second;
            if(newW+hubdis<QueryPLL(LID, hubid)){
                Label[LID][hubid]=newW+hubdis;
                Change[LID].push_back(make_pair(hubid, newW+hubdis));
                WaitPro.insert(LID);
            }
        }

        //check the label of their neighbors step by step
        while(WaitPro.size()>0){
            set<int> WaitProTem;
            vector<vector<pair<int,int>>> ChangeTem;
            vector<pair<int,int>> vec;
            ChangeTem.assign(node_num,vec);

            for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
                int curID=*it;
                vector<pair<int,int>> curChange=Change[curID];
                int neiID, neiDis, hID, hDis;
                for(int i=0;i<Neighbors[curID].size();i++){
                    neiID=Neighbors[curID][i].first;
                    neiDis=Neighbors[curID][i].second;

                    for(int j=0;j<curChange.size();j++){//for each decreased label, we check all its neighbors
                        hID=curChange[j].first; hDis=curChange[j].second;
                        if(NodeOrder[hID]>NodeOrder[neiID] && QueryPLL(neiID, hID)>neiDis+hDis){
                            Label[neiID][hID]=neiDis+hDis;
                            WaitProTem.insert(neiID);
                            ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
                        }
                    }
                }
            }

            WaitPro=WaitProTem;
            Change=ChangeTem;
        }
    }
}
//new version: set version with NoSupportedPair, clean label, 2023-05-14, iterative update, queue version, correct
void Graph::IncreasePSL(vertex a, vertex b, int oldW, int newW, vector<vector<pair<vertex,int>>> &Neighbors){//
    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
//            cout<<"Core "<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
//            cout<<"Core "<<b<<" "<<a<<" "<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbors[b][i].second=newW;
            break;
        }
    }

    int LID, HID;
    if(NodeOrder[a]>NodeOrder[b]){
        LID=b; HID=a;
    }else{
        LID=a; HID=b;
    }

    ///for debug
    bool ifDebug = false;
//    ifDebug = true;

    int lid, hid;
    lid=203827, hid=207273;//core 29-1

    int lid2,hid2;
    lid2=147082, hid2=196720;

//    cout<<"LID: "<<LID<<" ("<<NodeOrder[LID]<<") ; HID: "<<HID<<" ("<<NodeOrder[HID]<<")"<<endl;
    int dislower,disvally,dispeak;
    vertex peakhub,vallyID;
    pair<int,vertex> peakPair;//distance, hubID
    pair<int,vertex> vallyPair;//distance, vallyID
    set<tuple<vertex,vertex,vertex>> outdatedPruning;//<nodeID,supportNode,prunedID>
    map<pair<vertex,vertex>,vertex> newPruningPoints;//<<nodeID,prunedID>,supportNode>
    outdatedPruning.clear();
    newPruningPoints.clear();
    //activate or not
    dislower = DisQueryLower1(LID,HID,Neighbors);
    dispeak = DisQueryPeak(LID,HID);
//    cout<<"dis "<<dis<<" "<<oldW<<endl;

    if(dispeak<=oldW){
        if(ifDebug){
            cout<<"Not triggered! "<<LID<<" "<<HID<<": "<<dislower<<" "<<dispeak<<" "<<oldW<<" "<<newW<<endl;
        }
        return;
    }
    if((Label[LID].find(HID)!=Label[LID].end()) && (dislower>oldW)){//index update triggered
        if(ifDebug){
            cout<<"Triggered! "<<Label[LID][HID]<<endl;
        }


        set<tuple<vertex,vertex,vertex>> outdatedPruning;//<nodeID,supportNode,prunedID>
        map<pair<vertex,vertex>,vertex> newPruningPoints;//<<nodeID,prunedID>,supportNode>
        outdatedPruning.clear();
        newPruningPoints.clear();
        NoSupportedPair.clear();

        queue<pair<int,pair<int,int>>> WaitPro;
        queue<pair<int,pair<int,int>>> WaitProP;
        vector<pair<vertex,vertex>> AL1; AL1.clear();
        vector<pair<vertex,vertex>> AL2; AL2.clear();
        vector<pair<vertex,vertex>> AL2Check; AL2Check.clear();
        ChangedLabels.assign(node_num,unordered_set<int>());

        int curID, hubID, hubDis;
        int dis, cnt;



        /// weight change source 1, Label(a,b)
        vallyPair = DisQueryVally2(LID,HID,Neighbors);
        disvally=vallyPair.first, vallyID=vallyPair.second;
        peakPair = DisQueryPeak2(LID,HID);
        dispeak=peakPair.first, peakhub=peakPair.second;

        if(disvally >= dispeak){
            cout<<"! dispeak <= disvally for L("<<LID<<","<<HID<<"): "<<disvally<<" "<<dispeak<<endl;
        }


        if(Label[LID][HID] < disvally){
            AL1.emplace_back(LID,HID);
//            WaitPro.push(make_pair(LID, make_pair(HID,Label[LID][HID])));
//            Label[LID][HID]=disvally;//correct to the new value
//            ChangedLabels[LID].insert(HID);
        }

//        cout<<"Begin iterative propagation..."<<endl;
        int iterations=1;

        while(!AL1.empty()){
            if(iterations > 1){
                cout<<"Iteration "<<iterations<<endl;
            }
            iterations++;

            //update AL1: AL1->AL1
            CoarseUpdate(LID, HID, oldW, WaitPro, WaitProP, AL1, AL2, AL2Check, Neighbors, ifDebug, lid, hid);
//            cout<<"After AL1: "<<AL1.size()<<" "<<AL2Check.size()<<" "<<WaitPro.size()<<" "<<WaitProP.size()<< endl;


            //update AL2: AL1->AL2, AL2->AL2
            RefineUpdate( WaitPro, WaitProP, AL1, AL2, AL2Check, outdatedPruning, newPruningPoints, Neighbors, ifDebug, lid, hid);
//            cout<<"After AL2: "<<AL1.size()<<" "<<AL2Check.size()<<" "<<WaitPro.size()<<" "<<WaitProP.size()<<endl;
        }
//        cout<<"Done."<<endl;

//        cout<<"Remove and add pruning point."<<endl;
        PPRClean(Neighbors, newPruningPoints, outdatedPruning, ifDebug, lid, hid);

    }
    else{
        if(ifDebug){
            cout<<"Not triggered! "<<LID<<" "<<HID<<": "<<dislower<<" "<<dispeak<<" "<<oldW<<" "<<newW<<endl;
        }
    }

}

void Graph::PPRClean(vector<vector<pair<vertex,int>>> &Neighbors, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid){
    int dislower,disvally,dispeak;
    vertex peakhub,vallyID;
    pair<int,vertex> peakPair;//distance, hubID
    pair<int,vertex> vallyPair;//distance, vallyID
    int id1, hub, id2;
    ///remove old pruning point
    for(auto it2=outdatedPruning.begin();it2!=outdatedPruning.end();++it2){
        id1 = get<0>(*it2); id2 = get<2>(*it2); hub = get<1>(*it2);

        disvally=DisQueryVally(id1, id2,Neighbors);
        peakPair=DisQueryPeak2(id1, id2);
        dispeak=peakPair.first; peakhub=peakPair.second;

        if((disvally < dispeak)){
            if(ifDebug){
                if(lid == id1 && hid == id2){
                    cout<<"Remove pruning point: ("<<id1<<","<<id2<<"): hub "<<hub<<", "<<dispeak<<" "<<disvally<<endl;
                }
            }
            if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                PruningPointSet2[id1].erase(id2);
            }
            if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                PruningPointSet2[id2].erase(id1);
            }
            if(PruningPointSet[id1].find(hub) != PruningPointSet[id1].end()){//if found
                if(PruningPointSet[id1][hub].find(id2) != PruningPointSet[id1][hub].end() ){//if found
                    PruningPointSet[id1][hub].erase(id2);
                    if(PruningPointSet[id1][hub].empty()){
                        PruningPointSet[id1].erase(hub);
                    }
                }
            }
            if(PruningPointSet[id2].find(hub) != PruningPointSet[id2].end()){//if found
                if(PruningPointSet[id2][hub].find(id1) != PruningPointSet[id2][hub].end() ){//if found
                    PruningPointSet[id2][hub].erase(id1);
                    if(PruningPointSet[id2][hub].empty()){
                        PruningPointSet[id2].erase(hub);
                    }
                }
            }
        }else {//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for outdatedPruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
            PruningPointSet[id1][peakhub].insert(id2);
            PruningPointSet[id2][peakhub].insert(id1);
            PruningPointSet2[id1][id2]=peakhub;
            PruningPointSet2[id2][id1]=peakhub;
            if(Label[id1].find(id2) != Label[id1].end()){//if found
//                    cout<<"Outdated Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                Label[id1].erase(id2);
            }
        }
    }
    outdatedPruning.clear();
    ///add new pruning point
    for(auto it=newPruningPoints.begin();it!=newPruningPoints.end();++it){
        id1 = it->first.first; id2 = it->first.second; hub = it->second;
        vallyPair=DisQueryVally2(id1, id2,Neighbors);
        disvally=vallyPair.first; vallyID=vallyPair.second;
        peakPair=DisQueryPeak2(id1, id2);
        dispeak=peakPair.first; peakhub=peakPair.second;
        if(disvally < dispeak){//it is normal for this case as we add all peak pairs to newPruningPoint which may not be correct peak pairs
            if((Label[id1].find(id2) == Label[id1].end()) || (Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally)) {//if not found or found but incorrect
                if(ifDebug){
                    int disDijk = Dijkstra(id1,id2);
                    if(disvally == disDijk){
                        if(Label[id1].find(id2) != Label[id1].end()){//if found
                            cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") "<<Label[id1][id2]<<" "<<disvally<<endl;
                        }else{//if not found
                            cout<<"add label by new pruning point. L("<<id1<<","<<id2<<") -1 "<<disvally<<endl;
                        }
                        //                    Label[id1][id2] = disvally;
                    }else if(dispeak == disDijk){
                        cout<<"!! dispeak is correct. "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak;
                        if(Label[id1].find(id2) != Label[id1].end()) {//if found
                            cout<<" !! erase Label("<<id1<<","<<id2<<")";
                            Label[id1].erase(id2);
                        }
                        cout<<endl;
                    }else{

                        if(Label[id1].find(id2) != Label[id1].end()){
                            cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): " <<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                        }else{
                            cout<<"!!!!!!!! Totally wrong !!!!!!!! L("<<id1<<","<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<disDijk<<endl;
                        }
                    }
                }
                if(Label[id1].find(id2) == Label[id1].end()) {//if not found or found but incorrect
                    cout<<"missing label of outdatedPruning. L("<<id1<<" "<<id2<<"): -1 "<<disvally<<" "<<dispeak<<" "<<endl;
                }else if(Label[id1].find(id2) != Label[id1].end() && Label[id1][id2] != disvally){
                    cout<<"incorrect label of outdatedPruning. L("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<" "<<endl;
                }
            }

        }else{//if disvally >= dispeak
//                if(hub != peakhub){
//                    cout<<"hub is inconsistent for new pruning point ("<<id1<<" "<<id2<<") ! "<<hub<<" "<<peakhub<<endl;
//                }
            PruningPointSet[id1][peakhub].insert(id2);
            PruningPointSet[id2][peakhub].insert(id1);
            PruningPointSet2[id1][id2]=peakhub;
            PruningPointSet2[id2][id1]=peakhub;
            if(Label[id1].find(id2) != Label[id1].end()) {//if found
//                    cout<<"New Pruning point erase wrong. "<<id1<<" "<<id2<<" "<<Label[id1][id2]<<" "<<disvally<<" "<<dispeak<<"("<<peakhub<<","<<hub<<") "<<DijkstraCore(id1,id2)<<endl;
                Label[id1].erase(id2);
            }

        }

    }
    /// remove redundant label
    //cout<<"remove redundant labels..."<<endl;
    //        CleanLabel(Label);
    vector<pair<int,int>> eraseLabels; eraseLabels.clear();
    for(int i=node_num-1;i>=0;--i){
        id1 = vNodeOrder[i];
        if(!ChangedLabels[id1].empty()){
            for(auto it=ChangedLabels[id1].begin();it!=ChangedLabels[id1].end();++it){
                id2 = *it;
                vallyPair=DisQueryVally2(id1, id2,Neighbors);
                disvally=vallyPair.first; vallyID=vallyPair.second;
                peakPair=DisQueryPeak2(id1, id2);
                dispeak=peakPair.first; peakhub=peakPair.second;
                if(dispeak <= disvally){
                    PruningPointSet[id1][peakhub].insert(id2);
                    PruningPointSet[id2][peakhub].insert(id1);
                    PruningPointSet2[id1][id2]=peakhub;
                    PruningPointSet2[id2][id1]=peakhub;
                    eraseLabels.emplace_back(id1,id2);
//                    if(Label[id1].find(id2) != Label[id1].end()){//if found
////                        cout<<"erase redundant label ("<<id1<<" "<<id2<<"): "<<Label[id1][id2]<<" "<<dispeak<<endl;
//                        Label[id1].erase(id2);
//                    }
                }
                else{//if disvally < dispeak
                    if(PruningPointSet2[id1].find(id2) != PruningPointSet2[id1].end()){
                        PruningPointSet2[id1].erase(id2);
                    }
                    if(PruningPointSet2[id2].find(id1) != PruningPointSet2[id2].end()){
                        PruningPointSet2[id2].erase(id1);
                    }
                    if(PruningPointSet[id1].find(hub) != PruningPointSet[id1].end()){//if found
                        if(PruningPointSet[id1][hub].find(id2) != PruningPointSet[id1][hub].end() ){//if found
                            PruningPointSet[id1][hub].erase(id2);
                            if(PruningPointSet[id1][hub].empty()){
                                PruningPointSet[id1].erase(hub);
                            }
                        }
                    }
                    if(PruningPointSet[id2].find(hub) != PruningPointSet[id2].end()){//if found
                        if(PruningPointSet[id2][hub].find(id1) != PruningPointSet[id2][hub].end() ){//if found
                            PruningPointSet[id2][hub].erase(id1);
                            if(PruningPointSet[id2][hub].empty()){
                                PruningPointSet[id2].erase(hub);
                            }
                        }
                    }
                    if(Label[id1].find(id2) == Label[id1].end()){//if not found
                        cout<<"missing label! L("<<id1<<","<<id2<<"): -1 "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                    }else if(Label[id1][id2] != disvally){
                        cout<<"incorrect vally label! L("<<id1<<","<<id2<<"): "<<Label[id1][id2]<<" "<<disvally<<endl;
//                        Label[id1][id2]=disvally;
                    }
                }
            }
        }
    }
    if(!eraseLabels.empty()){
        cout<<"eraseLabels is not empty! ";
        bool flag_erase=false;
        for(auto it2=eraseLabels.begin();it2!=eraseLabels.end();++it2){
            id1 = it2->first, id2 = it2->second;
            if(Label[id1].find(id2) != Label[id1].end()){//if found
                Label[id1].erase(id2);
                flag_erase = true;
            }
        }
//            cout<<"label erased."<<endl;
    }
}


//function of propagating the AL1 labels, output the updated labels for AL2 check, queue version, 2023-03-26, new
void Graph::CoarseUpdate(vertex LID, vertex HID, int oldW, queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<vertex,vertex>>& AL1, vector<pair<vertex,vertex>>& AL2, vector<pair<vertex,vertex>>& AL2Check, vector<vector<pair<vertex,int>>> &Neighbor, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID

    int curID,hubID,hubDis;
    int neiID, neiDis;

    int lid2,hid2;
    lid2=147082, hid2=196720;

    /// weight change source 1, Label(a,b)
//    vallyPair = DisQueryVallyVert2(LID,HID,Neighbor,Label);
    vallyPair = DisQueryVally2(LID,HID,Neighbor);
    disvally=vallyPair.first, vallyID=vallyPair.second;
    peakPair = DisQueryPeak2(LID,HID);
    dispeak=peakPair.first, peakhub=peakPair.second;

    if(disvally >= dispeak){
        cout<<"! dispeak <= disvally for L("<<LID<<","<<HID<<"): "<<disvally<<" "<<dispeak<<endl;
    }


    if(Label[LID][HID] < disvally){
//        AL1.emplace_back(LID,HID);
        WaitPro.push(make_pair(LID, make_pair(HID,Label[LID][HID])));
        Label[LID][HID]=disvally;//correct to the new value
        ChangedLabels[LID].insert(HID);
        AL2Check.emplace_back(LID,HID);
    }

    /// weight change source 2, the labels that are directly affected by the e(a,b), using asynchronous propagation
    //Detect the affected labels of HID
    for(auto it=Label[LID].begin();it!=Label[LID].end();it++){
        hubID=(*it).first; hubDis=(*it).second;
        if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hubDis==Label[HID][hubID]){
            AL1.emplace_back(HID,hubID);
            vallyPair= DisQueryVally2(HID,hubID,Neighbor);
            disvally=vallyPair.first, vallyID= vallyPair.second;

            if(Label[HID][hubID] < disvally){
                if(ifDebug){
                    if((lid == HID && hid ==hubID) || (lid2 == HID && hid2 == hubID)){
                        cout<<"Find. "<<HID<<" "<<hubID<<" "<<Label[HID][hubID]<<" "<<disvally<<endl;
                    }
                }
                WaitPro.push(make_pair(HID, make_pair(hubID,Label[HID][hubID])));

                Label[HID][hubID]=disvally;//correct to the new value
                ChangedLabels[HID].insert(hubID);
                //AL1->AL2, check PPR
                AL2Check.emplace_back(HID,hubID);
            }
        }
    }
    //Detect the affected labels of LID
    for(auto it=Label[HID].begin();it!=Label[HID].end();it++){
        hubID=(*it).first; hubDis=(*it).second;
        if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hubDis==Label[LID][hubID]){
            AL1.emplace_back(LID,hubID);
            vallyPair= DisQueryVally2(LID,hubID,Neighbor);
            disvally=vallyPair.first, vallyID= vallyPair.second;

            if(Label[LID][hubID] < disvally){
                if(ifDebug){
                    if((lid == LID && hid ==hubID) || (lid2 == LID && hid2 == hubID)){
                        cout<<"Find. "<<LID<<" "<<hubID<<" "<<Label[LID][hubID]<<" "<<disvally<<endl;
                    }
                }
                WaitPro.push(make_pair(LID, make_pair(hubID,Label[LID][hubID])));

                Label[LID][hubID]=disvally;//correct to the new value
                ChangedLabels[LID].insert(hubID);
                //AL1->AL2, check PPR
                AL2Check.emplace_back(LID,hubID);
            }
        }
    }

    //Update the source labels of AL1
//    for(auto it=AL1.begin();it!=AL1.end();++it){
//        curID = it->first, hubID = it->second;
//        vallyPair= DisQueryVally2(curID,hubID,Neighbor,Label);
//        disvally=vallyPair.first, vallyID= vallyPair.second;
//
//        if((Label[curID].find(hubID) != Label[curID].end()) &&  (Label[curID][hubID] < disvally)){
//            if(ifDebug){
//                if((lid == curID && hid ==hubID) || (lid2 == curID && hid2 == hubID)){
//                    cout<<"Find. "<<curID<<" "<<hubID<<" "<<Label[curID][hubID]<<" "<<disvally<<endl;
//                }
//            }
//            WaitPro.push(make_pair(curID, make_pair(hubID,Label[curID][hubID])));
//
//            Label[curID][hubID]=disvally;//correct to the new value
//            ChangedLabels[curID].insert(hubID);
//            //AL1->AL2, check PPR
//            AL2Check.emplace_back(curID,hubID);
//        }
//    }
    AL1.clear();

    while(!WaitPro.empty()){
        queue<pair<int,pair<int,int>>> WaitProTem;
        while(!WaitPro.empty()){
            auto front = WaitPro.front();
            curID = front.first, hubID = front.second.first, hubDis = front.second.second;
            WaitPro.pop();
//            cout<<"AL("<<curID<<","<<hubID<<"): "<<hubDis<<endl;

            //lower-order vertex
            for (int k = 0; k < Neighbor[curID].size(); k++) {//for each neighbor of curID which has the hub hID
                neiID = Neighbor[curID][k].first;
                neiDis = Neighbor[curID][k].second;
//                if (ifDebug) {
//                        if (neiID == lid && hubID == hid) {
//                            cout << "AL1->AL1 detection. curID " << curID << " L(" << neiID << "," << hubID << ")"
//                                 << ": ";
//                            if (Label[neiID].find(hubID) != Label[neiID].end()) {//if found
//                                cout << Label[neiID][hubID] << " ";
//                            }
//                            cout << DisQueryVallyDebug(neiID, hubID, Neighbor) << " "<<DijkstraCore(neiID,hubID)<< endl;
//                        }
//                }
                //detect and update the affected label of neighbor
                if ((Label[neiID].find(hubID) != Label[neiID].end()) && (neiDis + hubDis == Label[neiID][hubID])) {
//                        AL1.insert(OrderComp3(neiID, hubID));
                    vallyPair = DisQueryVally2(neiID, hubID, Neighbor);
                    disvally = vallyPair.first, vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(neiID, hubID);
                    dispeak = peakPair.first, peakhub = peakPair.second;

                    //update the label
                    if (Label[neiID][hubID] < disvally) {//original setting
                        WaitProTem.push(make_pair(neiID, make_pair(hubID,Label[neiID][hubID])));

                        Label[neiID][hubID] = disvally;
                        ChangedLabels[neiID].insert(hubID);
                        //AL1->AL2, check PPR
//                        AL2Check.emplace_back(neiID,hubID);
                    }
                    //AL1->AL2, check PPR
                    AL2Check.emplace_back(neiID,hubID);
                }
            }

        }

        WaitPro = WaitProTem;
    }

}
//function of propagating the AL2 labels, output the updated labels for AL1 check, queue version, new
void Graph::RefineUpdate(queue<pair<int,pair<int,int>>>& WaitPro, queue<pair<int,pair<int,int>>>& WaitProP, vector<pair<vertex,vertex>>& AL1, vector<pair<vertex,vertex>>& AL2, vector<pair<vertex,vertex>>& AL2Check, set<tuple<vertex,vertex,vertex>>& outdatedPruning, map<pair<vertex,vertex>,vertex>& newPruningPoints, vector<vector<pair<vertex,int>>> &Neighbor, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID

    int curID,hubID,hubDis;

    bool ifUpdate = false;

    //check and get the source AL2, AL1->AL2
    for(int i=0;i<AL2Check.size();++i){//for each increased label
        curID = AL2Check[i].first, hubID = AL2Check[i].second;

        bool temp_bool = PPRCheck(curID,hubID,Neighbor,WaitProP,AL2,newPruningPoints,outdatedPruning,ifDebug,lid,hid);///
        if(temp_bool){
            ifUpdate = true;
        }

//        WaitProP.push(make_pair(curID, make_pair(hubID,-1)));///

//        if(ifUpdate){//if update
//            vallyPair=DisQueryVally2(curID,hubID,Neighbor,Label);
//            disvally=vallyPair.first, vallyID=vallyPair.second;
//            peakPair=DisQueryPeak2(curID,hubID,Label);
//            dispeak=peakPair.first, peakhub=peakPair.second;
//
//            if(disvally < dispeak){
//                int Dijk = DijkstraCore(curID,hubID);
//                if(Dijk != Label[curID][hubID]){
//                    cout<<"Incorrect L("<<curID<<","<<hubID<<"): "<<Label[curID][hubID]<<" "<<disvally<<"("<<vallyID<<") "<< dispeak<<"("<<peakhub<<") "<< Dijk<<endl;
//                }
//            }
//            cout<<"new insertion caused by L("<<curID<<","<<hubID<<") "<<Label[curID][hubID]<<" "<<Dijk<<endl;
//        }

    }

    AL2Check.clear();
    AL1.clear();
    if(ifDebug){
        if(ifUpdate){
            cout<<"!! Have new labels inserted!"<<endl;
        }
        if(Label[lid].find(hid) != Label[lid].end()){
            cout<<lid<<" "<<hid<<": "<<Label[lid][hid]<<endl;
        }
        cout<<"AL2->AL2"<<endl;
    }

    int round=0;
    int u,neiid,neiw,cDis;

    while(!WaitProP.empty()){
        queue<pair<int,pair<int,int>>> WaitProPTem;
        while(!WaitProP.empty()){
//        cout<<"round "<<round<<endl;
//        round++;
            auto front = WaitProP.front();
            curID = front.first, u = front.second.first; //hubID=front.second.second;
            WaitProP.pop();

            //detect the affected neighbors of AL2
            //Lower-order vertex
            for(int l = 0; l < Neighbor[curID].size(); ++l) {
                neiid = Neighbor[curID][l].first;
                neiw = Neighbor[curID][l].second;
                if (NodeOrder[neiid] < NodeOrder[u]) {
                    vallyPair = DisQueryVally2(neiid, u, Neighbor);
                    disvally = vallyPair.first; vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(neiid, u);
                    dispeak = peakPair.first; peakhub = peakPair.second;//dispeak may be incorrect for the time being
                    bool flag = false;
                    if (PruningPointSet2[neiid].find(u) != PruningPointSet2[neiid].end()) {
                        flag = true;
                        hubID = PruningPointSet2[neiid][u];
                    }else{
                        hubID = -1;
                    }


                    if (Label[neiid].find(u) == Label[neiid].end()) {// if not found
                        if (disvally < dispeak) {//if disvally is smaller, add new label
//                        WaitProPTem.insert(neiid);
//                        ChangePTem[neiid].insert(u);
                            WaitProPTem.push(make_pair(neiid, make_pair(u,-1)));

                            AL2.emplace_back(neiid, u);
                            Label[neiid][u] = disvally;
                            ChangedLabels[neiid].insert(u);
                            NoSupportedPair.insert(make_pair(neiid, u));///
//                       AL2Checked.insert(make_pair(neiid,u));
//                       cout<<"AL2->AL2 new insertion. L("<<neiid<<","<<u<<"): "<<Label[neiid][u]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(neiid,u)<<endl;

                            if(flag){
                                outdatedPruning.insert(make_tuple(neiid, hubID, u));//
                            }


                            if (PruningPointSet2[neiid].find(u) != PruningPointSet2[neiid].end()) {///
                                PruningPointSet2[neiid].erase(u);
                            }
                            if (PruningPointSet2[u].find(neiid) != PruningPointSet2[u].end()) {
                                PruningPointSet2[u].erase(neiid);
                            }
                        }
                        else {//if dispeak<=disvally
                            if (peakhub != -1 && peakhub != hubID) {//may add redundant PruningPoints
                                newPruningPoints[make_pair(neiid,u)] = peakhub;// this may not be true as the dispeak may be wrong
                                if (hubID != -1) {
                                    outdatedPruning.insert(make_tuple(neiid, hubID, u));//
                                }
                                PruningPointSet2[neiid][u] = peakhub;
                                PruningPointSet2[u][neiid] = peakhub;
                                PruningPointSet[neiid][peakhub].insert(u);
                                PruningPointSet[u][peakhub].insert(neiid);
                            }
                        }

                    }
                    else if ((Label[neiid].find(u) != Label[neiid].end()) && (Label[neiid][u] != disvally)){//if found, check whether the new label L(curID,u) will affect L(neiid,u)


                        if(Label[neiid][u] > disvally){//for decrease update
                            AL2.emplace_back(neiid, u);
                            WaitProPTem.push(make_pair(neiid, make_pair(u,-1)));///AL2
                            Label[neiid][u]=disvally;
                            ChangedLabels[neiid].insert(u);
                        }else if(Label[neiid][u] < disvally){//for increase update
                            AL1.emplace_back(neiid,u);
                        }

                    }

                }
            }

            //Higher-order vertex
            /*for(int l=0; l<Neighbor[u].size();++l){
                neiid = Neighbor[u][l].first;
                neiw = Neighbor[u][l].second;
                if (NodeOrder[curID] < NodeOrder[neiid]) {
                    vallyPair = DisQueryVally2(curID, neiid, Neighbor);
                    disvally = vallyPair.first; vallyID = vallyPair.second;
                    peakPair = DisQueryPeak2(curID, neiid);
                    dispeak = peakPair.first; peakhub = peakPair.second;//dispeak may be incorrect for the time being
                    bool flag = false;
                    if (PruningPointSet2[curID].find(neiid) != PruningPointSet2[curID].end()) {
                        hubID = PruningPointSet2[curID][neiid];//original peak hub
                        flag = true;
                    } else {
                        hubID = -1;
                    }

                    if (Label[curID].find(neiid) == Label[curID].end()) {// if not found
                        if (disvally < dispeak) {//if disvally is smaller, add new label
                            WaitProPTem.push(make_pair(curID, make_pair(neiid,-1)));

                            AL2.emplace_back(curID,neiid);
                            Label[curID][neiid] = disvally;
                            ChangedLabels[curID].insert(neiid);
                            NoSupportedPair.insert(make_pair(curID,neiid));

                            if (hubID != -1) {
                                outdatedPruning.insert(make_tuple(curID, hubID, neiid));//
                            } else {//hub == -1, if there is no peak hub in PruningPointSet2, but the label from u to neiid may truly pruned
                                outdatedPruning.insert(make_tuple(curID, peakhub, neiid));//
                            }

                            if (PruningPointSet2[neiid].find(curID) != PruningPointSet2[neiid].end()) {///
                                PruningPointSet2[neiid].erase(curID);
                            }
                            if (PruningPointSet2[curID].find(neiid) != PruningPointSet2[curID].end()) {
                                PruningPointSet2[curID].erase(neiid);
                            }
                        }
                        else {//if dispeak<=disvally
                            if (peakhub != -1 && peakhub != hubID) {//may add redundant PruningPoints
                                newPruningPoints[make_pair(curID,neiid)] = peakhub;// this may not be true as the dispeak may be wrong
                                if (hubID != -1) {
                                    outdatedPruning.insert(make_tuple(curID, hubID, neiid));//
                                }
                                PruningPointSet2[neiid][curID] = peakhub;
                                PruningPointSet2[curID][neiid] = peakhub;
                                PruningPointSet[neiid][peakhub].insert(curID);
                                PruningPointSet[curID][peakhub].insert(neiid);
                            }
                        }

                    }
                    else {//if found, check whether the new label L(curID,u) will affect L(neiid,u)

                        if (Label[curID][neiid] != disvally) {
//                            if(ifDebug){
//                                cout<<"AL2->AL2 higher existing label. "<<curID<<"("<<NodeOrder[curID]<<") L("<<neiid<<","<<u<<"): "<<Label[curID][neiid]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<DijkstraCore(curID,neiid)<<endl;
//                            }

                            if(Label[curID][neiid] > disvally){//for decrease update
                                AL2.emplace_back(curID,neiid);
                                WaitProPTem.push(make_pair(curID, make_pair(neiid,-1)));///AL2
                                Label[curID][neiid]=disvally;
                                ChangedLabels[curID].insert(neiid);
                            }else if(Label[curID][neiid] < disvally){//for increase update
                                AL1.emplace_back(curID,neiid);
                                cout<<"!!! Unexpect result: Enter this branch!! L<d_1 for Refine update! L("<<curID<<","<<neiid<<"): "<<Label[curID][neiid]<<" "<<disvally<< endl;

                            }
                        }
                    }

                }
            }*/

        }
        WaitProP=WaitProPTem;
    }

}
//function of checking the pruned label by PPR, pair version, queue version, correct
bool Graph::PPRCheck(vertex curID, vertex hubID, vector<vector<pair<vertex,int>>> &Neighbor, queue<pair<int,pair<int,int>>>& WaitProPTem, vector<pair<vertex,vertex>> &AL2, map<pair<vertex,vertex>,vertex> &newPruningPoints, set<tuple<vertex,vertex,vertex>> &outdatedPruning, bool ifDebug, int lid, int hid){
    int disvally,dispeak,peakhub,vallyID;
    pair<int,int> peakPair;//distance, hubID
    pair<int,int> vallyPair;//distance, vallyID
    bool flag = false;

//    if(ifDebug){
//        if(curID==211099 && hubID==211425){
//            cout<<"PPR source find. L("<<curID<<","<<hubID<<"): "<<Label[curID][hubID]<<endl;
//        }
//    }

    if(PruningPointSet[curID].find(hubID)!=PruningPointSet[curID].end()){
        for(auto snum=PruningPointSet[curID][hubID].begin();snum!=PruningPointSet[curID][hubID].end();++snum){//for each pruned vertex
            int s=*snum;
            if((NodeOrder[s]<NodeOrder[curID]) && (NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end())){//not found
//            if(NodeOrder[s]<NodeOrder[curID]){//not found
//                AL2.insert(OrderComp3(s,curID));///
                vallyPair=DisQueryVally2(s,curID,Neighbor);
                disvally=vallyPair.first; vallyID=vallyPair.second;//disvally may not be correct if there exist neighbor which is also the pruning point of (curID,hID)
                peakPair=DisQueryPeak2(s,curID);//the original dispeak is smaller than disvally
                dispeak=peakPair.first; peakhub=peakPair.second;


                if(disvally<dispeak){//the pruned label should be inserted back
//                    WaitProPTem.insert(s);
//                    ChangePTem[s].insert(curID);
                    WaitProPTem.push(make_pair(s, make_pair(curID,hubID)));

                    if(ifDebug){
//                        int Dijk= DijkstraCore(s,curID);
//                        if(Dijk != disvally){
//                            if(Label[s].find(curID) != Label[s].end()){//if found
//                                cout<<"AL1->AL2. L("<<s<<","<<curID<<") "<<Label[s][curID]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }else{//if not found
//                                cout<<"AL1->AL2. L("<<s<<","<<curID<<") -1 "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }
//                        }
                    }

                    AL2.emplace_back(s,curID);
                    Label[s][curID]=disvally;
                    ChangedLabels[s].insert(curID);
                    flag = true;

                    outdatedPruning.insert(make_tuple(s,hubID,curID));
                    NoSupportedPair.insert(make_pair(s,curID));

                    if(PruningPointSet2[curID].find(s) != PruningPointSet2[curID].end()){///
                        PruningPointSet2[curID].erase(s);
                    }
                    if(PruningPointSet2[s].find(curID) != PruningPointSet2[s].end()){
                        PruningPointSet2[s].erase(curID);
                    }

                }
                else {//if dispeak<=disvally

                    if(peakhub != -1 && peakhub != hubID){
                        newPruningPoints[make_pair(s,curID)] = peakhub;
                        outdatedPruning.insert(make_tuple(s,hubID,curID));
                        PruningPointSet2[curID][s]=peakhub;
                        PruningPointSet2[s][curID]=peakhub;
                        PruningPointSet[curID][peakhub].insert(s);
                        PruningPointSet[s][peakhub].insert(curID);
                    }
                }
            }else if(NodeOrder[curID]<NodeOrder[s] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
//            } else if(NodeOrder[s]>NodeOrder[curID]){
//                AL2.insert(OrderComp3(curID,s));///
                vallyPair=DisQueryVally2(curID,s,Neighbor);
                disvally=vallyPair.first; vallyID=vallyPair.second;
                peakPair=DisQueryPeak2(curID,s);//the old P2 distance dominates P1 distance, but the increase of d(curID,hID)
                dispeak=peakPair.first; peakhub=peakPair.second;


                if(disvally<dispeak){//disvally may be incorrect
//                    WaitProPTem.insert(curID);///
//                    ChangePTem[curID].insert(s);
                    WaitProPTem.push(make_pair(curID, make_pair(s,hubID)));

                    if(ifDebug){
//                        int Dijk= DijkstraCore(s,curID);
//                        if(Dijk != disvally){
//                            if(Label[curID].find(s) != Label[curID].end()){//if found
//                                cout<<"AL1->AL2. L("<<curID<<","<<s<<") "<<Label[curID][s]<<" "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }else{//if not found
//                                cout<<"AL1->AL2. L("<<curID<<","<<s<<") -1 "<<disvally<<"("<<vallyID<<") "<<dispeak<<"("<<peakhub<<") "<<Dijk<<endl;
//                            }
//                        }
                    }

                    AL2.emplace_back(curID,s);
                    Label[curID][s]=disvally;//the label here may be incorrect for the time being
                    ChangedLabels[curID].insert(s);
                    flag = true;

                    outdatedPruning.insert(make_tuple(curID,hubID,s));
                    NoSupportedPair.insert(make_pair(curID,s));

                    if(PruningPointSet2[curID].find(s) != PruningPointSet2[curID].end()){
                        PruningPointSet2[curID].erase(s);
                    }
                    if(PruningPointSet2[s].find(curID) != PruningPointSet2[s].end()){
                        PruningPointSet2[s].erase(curID);
                    }
                }
                else {//if dispeak<=disvally

                    if(peakhub != -1 && peakhub != hubID) {
                        newPruningPoints[make_pair(curID,s)] = peakhub;
                        outdatedPruning.insert(make_tuple(curID, hubID, s));///
                        PruningPointSet2[curID][s]=peakhub;
                        PruningPointSet2[s][curID]=peakhub;
                        PruningPointSet[curID][peakhub].insert(s);
                        PruningPointSet[s][peakhub].insert(curID);
                    }
                }
            }
        }
    }
    return flag;
}

//function of computing distance from ID1(LID) to ID2(HID), via a neighbor of ID1 (which has lower order than HID) on updated graph
int Graph::DisQueryLower1(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors){
    vertex neiID;
    int neiDis,d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if neiID has lower order and ID2 is a hub of neiID
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
            }
        }
    }
    return d;
}
//old version: function of computing d2
int Graph::DisQueryPeak(vertex ID1, vertex ID2){
    int d=INF;

    vertex hub, dis1;
    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
            if(dis1+Label[ID2][hub]<d){
                d=dis1+Label[ID2][hub];
            }
        }
    }
    return d;
}
//old version. function of computing the shortest distance through the label of ID1's neighbors, i.e., d1.
int Graph::DisQueryVally(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors){
    vertex neiID;
    int neiDis,d=INF;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){//if found
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
            }
        }
    }
    return d;
}
//function of computing d2, new version
pair<int,vertex> Graph::DisQueryPeak2(vertex ID1, vertex ID2){
    int d=INF;
    unordered_map<int,int>::iterator it;
    vertex hub, dis1, finalHub=-1;

    for(auto it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){//hub is the hub of ID1 and ID2 and hub != ID2
            if(dis1+Label[ID2][hub]<d){
                d=dis1+Label[ID2][hub];
                finalHub = hub;
            }

        }
    }
    return make_pair(d,finalHub);
}
//function of computing vally distance, return distance and neighbor
pair<int,vertex> Graph::DisQueryVally2(vertex ID1, vertex ID2, vector<vector<pair<vertex,int>>> &Neighbors){
    int neiDis;
    int d=INF;
    vertex neiID,finalNeigh = -1;
    for(int i=0;i<Neighbors[ID1].size();i++){
        neiID=Neighbors[ID1][i].first;
//        if(neiID == 204062){
//            cout<<ID1<<" "<<neiID<<" "<<ID2<<endl;
//            if(PruningPointSet2[neiID].find(ID2)!=PruningPointSet2[neiID].end()){
//                cout<<"Pruned "<<PruningPointSet2[neiID][ID2]<<endl;
//            }
//        }
        neiDis=Neighbors[ID1][i].second;
        if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
            if(neiDis+Label[neiID][ID2]<d){
                d=neiDis+Label[neiID][ID2];
                finalNeigh = neiID;
//                cout<<neiID<<" "<<d<<endl;
            }
        }
    }
    return make_pair(d,finalNeigh);
}

#endif //PLL_HPP