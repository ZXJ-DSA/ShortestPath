//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef H2H_HPP
#define H2H_HPP

#include "head.h"

/*-------------Index Construction Algorithms--------------*/
/// CH algorithm
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
void Graph::MDEOrderGenerate(string orderfile) {
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

    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);

    set<DegComp> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    unsigned int Twidth=0;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
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
            Deg.erase(DegComp(x));
            _DD_[x]=DD[x];
            Deg.insert(DegComp(x));
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
    int ifWrite=false;
//    ifWrite=true;
    if(ifWrite){
        cout<<"Writing vertex order..."<<endl;
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << "\t" << NodeOrder[i] << endl;
        ofile.close();
    }

    cout<<"Time for node ordering: "<<tt.GetRuntime()<<" s. treewidth "<<Twidth<<endl;
}
//Degree-based Ordering
void Graph::DegreeOrderGenerate(string orderfile){
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

    int ifWrite=false;
//    ifWrite=true;
    if(ifWrite){
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << "\t" << NodeOrder[i] << endl;
        ofile.close();
    }

    cout<<"Time for node ordering: "<<tt.GetRuntime()<<" s."<<endl;
}

int Graph::writeShortCutorder(string filename){
    //only write the one with higher order
    //write together with the higher adjacent vertex
    ofstream ofile(filename);
    cout<<"Writing shortcuts!"<<endl;
    for(int i = 0; i < node_num; i++)
    {
        ofile << i << "\t" << NodeOrder[i] << "\t";

        vector<pair<int, int>> shortcut;
        int ID, Wei;
        for(int k=0;k<vvpShortCut[i].size();k++){
            ID=vvpShortCut[i][k].first;
            Wei=vvpShortCut[i][k].second;
            if(NodeOrder[ID]>NodeOrder[i]){
                shortcut.push_back(make_pair(ID, Wei));
            }

        }

        for(int j=0;j<Neighbor[i].size();j++){
            ID=Neighbor[i][j].first;
            Wei=Neighbor[i][j].second;
            if(NodeOrder[ID]>NodeOrder[i]){
                shortcut.push_back(make_pair(ID, Wei));
            }

        }

        if(shortcut.size() != 0)
        {
            ofile << 1;
            ofile << "\t" << shortcut.size();
            vector<pair<int,int>>::iterator iR;
            iR=shortcut.begin();
            while(iR!=shortcut.end()){
                ofile<<"\t"<<(*iR).first<<"\t"<<(*iR).second;
                AdjaShort[i].push_back(make_pair((*iR).first,(*iR).second));
//                AdjaShortR[(*iR).first].push_back(make_pair(i,(*iR).second));
                iR++;
            }
        }
        else{
            ofile << 0;
        }

        ofile << endl;
    }
    return 0;
}

int Graph::ReadShortCut(string filename){
    ifstream inSCH(filename);
    int nodeID, nodeOrder, c, nsc, nrsc, n, d, m;
    cout<<"Reading shortcuts!"<<endl;

    for(int i = 0; i < node_num; i++)
    {
        inSCH >> nodeID >> nodeOrder >> c;
        NodeOrder[nodeID] = nodeOrder;
        if(!c)
            continue;

        inSCH >> nrsc;
        for(int j = 0; j < nrsc; j++)
        {
            inSCH >> n >> d;
            AdjaShort[nodeID].push_back(make_pair(n,d));
//            AdjaShortR[n].push_back(make_pair(nodeID,d));
        }
    }
    inSCH.close();

    cout<<"shortcut finish reading!"<<endl;
    return 0;
}

int Graph::CHContractionSearchMT(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult){
    benchmark::heap<2,int,int> Heap(node_num);
    vector<int> vDistance(node_num, INF);
    vector<int> vPre(node_num,-1);
    int topNodeID, topDistance, neighborNodeID, neighborLength;
    int w;
    vDistance[ID1]=0;
    Heap.update(ID1,0);
    map<int,int> mWDistance;
    vector<pair<int,int>> vpWDistance;
    map<int,int> mDistance;
    int maxWDistance=-1;

    for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
        w=(*ivp).first;
        int d=(*ivp).second+dUV;
        mWDistance[w]=d;
        if(d>maxWDistance)
            maxWDistance=d;
        mDistance[w]=INF;
    }

    if(mDistance.empty()){
        return 0;
    }

    for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
        vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

    int dThreshold=maxWDistance;

    while(!Heap.empty()){
        Heap.extract_min(topNodeID, topDistance);
        if(vbVisited[topNodeID])
            continue;
        if(topDistance>dThreshold)
            break;
        for(auto ivp=vvNode[topNodeID].begin();ivp!=vvNode[topNodeID].end();ivp++){
            neighborNodeID = (*ivp).first;
            neighborLength = (*ivp).second;
            if(vbVisited[neighborNodeID] || neighborNodeID==ID2)
                continue;
            int d = vDistance[topNodeID] + neighborLength;
            if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
            {
                vDistance[neighborNodeID] = d;
                Heap.update(neighborNodeID, d);
                vPre[neighborNodeID]=topNodeID;
                if(mWDistance.find(neighborNodeID) != mWDistance.end())
                    mDistance[neighborNodeID] = d;
            }else if(vDistance[neighborNodeID]>d){
                vDistance[neighborNodeID] = d;
                Heap.update(neighborNodeID, d);
                vPre[neighborNodeID]=topNodeID;
                if(mWDistance.find(neighborNodeID) != mWDistance.end())
                    mDistance[neighborNodeID] = d;
            }
            else if(vDistance[neighborNodeID] == d)
            {
                if(mWDistance.find(neighborNodeID) != mWDistance.end())
                    mDistance[neighborNodeID] = d;
            }
        }
    }

    for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
    {
        if((*imDistance).second > mWDistance[(*imDistance).first])//if Dijkstra distance is higher than the path pass through v, means that shortcut should be added
        {
            int w = (*imDistance).first;
            int distance = mWDistance[(*imDistance).first];

            vvpResult[k].push_back(make_pair(w, distance));
        }
    }

    return 0;
}
//CH shortcut construction
void Graph::CHConstructionMT() {
    vvpShortCut.assign(node_num, vector<pair<vertex ,int>>());
    vector<bool> vbVisited(node_num, false);
    map<int, int> mPosition;
    mPosition[0]=node_num-1;
    bool bUpdated;
    vector<pair<vertex,int>> vU,vW;
    //cout<<"Building static CH "<<endl;

    int v;// current contracting vertex
    int step=node_num/10000;
    step= max(step,1);
    step*=1000;


    for(int i=0;i<vNodeOrder.size()-1;i++){
        if(i%step==0){
            cout<<"Contracting vertex "<<i<<endl;
        }

        v=vNodeOrder[i];
        if(v!=-1){
            vU=vvNode[v];
            //cout<<"contracting node "<<v<<endl;

            //multiple thread
            vector<vector<pair<int,int>>> vvpResult;
            vector<pair<int,int>> vec;
            vec.clear();
            vvpResult.assign(vU.size(), vec);
            //filter the neighbors no need for contraction
            vector<vector<pair<int,int>>> vUeach;
            vUeach.assign(vU.size(), vec);
            int w;
            for(int l=0;l<vU.size();l++){
                int ID1=vU[l].first;
                for(auto ivp=vU.begin();ivp!=vU.end();ivp++){
                    w=(*ivp).first;
                    if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation, only deal with higher-order nodes
                        if(vbVisited[w])
                            continue;
                        vUeach[l].push_back(*ivp);//the higher-order neighbors

                    }
                }
            }

            //better not single thread
            /*if(vU.size()<10){
                //single thread
                for(auto ivpr=vU.begin();ivpr!=vU.end();ivpr++){
                    CHcontractionorder((*ivpr).first, v, vbVisited, (*ivpr).second, vU);
                }
            }else{*/
            boost::thread_group threadf;
            for(int k=0;k<vU.size();k++){
                if(vUeach[k].size()>0)
                    threadf.add_thread(new boost::thread(&Graph::CHContractionSearchMT, this, k, vU[k].first, v, boost::ref(vbVisited), vU[k].second, boost::ref(vUeach[k]), boost::ref(vvpResult)));
            }
            threadf.join_all();
            //merge to new shortcut and write it back
            for(int k=0;k<vU.size();k++)
            {
                int	ID1 = vU[k].first;
                for(int j=0;j<vvpResult[k].size();j++){
                    int w=vvpResult[k][j].first;
                    int distance=vvpResult[k][j].second;
                    vvpShortCut[ID1].push_back(make_pair(w, distance));
                    vvNode[ID1].push_back(make_pair(w, distance));
                    vvNode[w].push_back(make_pair(ID1, distance));
                }
            }
            //}

            //deleting v from G
            for(auto ivp = vvNode[v].begin(); ivp != vvNode[v].end(); ivp++)
            {
                for(auto ivpr = vvNode[(*ivp).first].begin(); ivpr != vvNode[(*ivp).first].end(); ivpr++)
                    if((*ivpr).first == v)
                    {
                        vvNode[(*ivp).first].erase(ivpr);
                        break;
                    }
            }
            vbVisited[v]=true;
        }
    }

    for(int i = 0; i < node_num; i++)
    {

        vector<pair<int, int>> shortcut;
        int ID, Wei;
        for(int k=0;k<vvpShortCut[i].size();k++){
            ID=vvpShortCut[i][k].first;
            Wei=vvpShortCut[i][k].second;
            if(NodeOrder[ID]>NodeOrder[i]){
                shortcut.push_back(make_pair(ID, Wei));
            }

        }

        for(int j=0;j<Neighbor[i].size();j++){
            ID=Neighbor[i][j].first;
            Wei=Neighbor[i][j].second;
            if(NodeOrder[ID]>NodeOrder[i]){
                shortcut.push_back(make_pair(ID, Wei));
            }

        }

        if(shortcut.size() != 0)
        {
            vector<pair<int,int>>::iterator iR;
            iR=shortcut.begin();
            while(iR!=shortcut.end()){
                AdjaShort[i].emplace_back((*iR).first,(*iR).second);
//                AdjaShortR[(*iR).first].push_back(make_pair(i,(*iR).second));
                iR++;
            }
        }

    }
}
//CH shortcut construction
void Graph::CHConstruction(int strategy, bool ifDynamic){
    string indexfile=graph_path+".shortcutCH";
    string orderfile=graph_path+".orderCH";
//    MDEOrderGenerate(orderfile);
    VertexOrdering(strategy);

    AdjaShort.assign(node_num, vector<pair<vertex,int>>());
//    AdjaShortR.assign(node_num, vector<pair<vertex,int>>());
    vvNode.assign(Neighbor.begin(), Neighbor.end());//duplicate of original adjacency list Neighbor

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;

    fstream file;
    file.open(indexfile);
//    if(!file){
    if(true){
        t1=std::chrono::high_resolution_clock::now();
        CHConstructionMT();
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Time for CH construction: "<<runT<<" s."<<endl;
//        writeShortCutorder(indexfile);
//        ReadShortCut(indexfile);
    }
    else
    {
        ReadShortCut(indexfile);
    }
    IndexsizeCH();
}

//CH algorithm with pruning
int Graph::QueryCH(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            for(auto out=AdjaShort[topNodeIDForward].begin();out!=AdjaShort[topNodeIDForward].end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

            for(auto in=AdjaShort[topNodeIDBackward].begin();in!=AdjaShort[topNodeIDBackward].end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//function for computing the index size
void Graph::IndexsizeCH(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<AdjaShort.size();i++){
        m1+=AdjaShort[i].size()*2*sizeof(int);//dis
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

/// CH algorithm without pruning
void Graph::CHConstructionWP(int strategy, bool ifDynamic){

    string indexfile=graph_path+".shortcutCH";
    string orderfile=graph_path+".orderCH";
//    MDEOrderGenerate(orderfile);
    VertexOrdering(strategy);

    int num, nodeID, nodeorder;

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);

    map<int, vector<int>> mi;
    SCconNodesMT.assign(node_num, mi);

    Timer tt;
    tt.start();
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }

    vector<bool> exist; exist.assign(node_num,true);
    //vector<bool> change; change.assign(nodenum,false);

    //cout<<"Begin to contract"<<endl;
    int stepShow=node_num/5;
    for(int nodeorder=0;nodeorder<node_num;nodeorder++){//start from the most important vertex
        int x=vNodeOrder[nodeorder];
        if(x!=-1){//to identify and exclude the isolated vertices
            if(nodeorder%stepShow==0){
                cout<<"Contract node order: "<<nodeorder<<endl;
            }


            exist[x]=false;
            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.emplace_back(*it);
                }
            }
            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteEorder(x,y);
                //change[y]=true;
            }

            if(Neigh.size()<=100){
                //single thread
                for(int i=0;i<Neigh.size();i++){
                    for(int j=i+1;j<Neigh.size();j++){
                        insertEorder(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                        if(Neigh[i].first<Neigh[j].first)
                            SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
                        else if(Neigh[j].first<Neigh[i].first)
                            SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);
                    }
                }
            }else{
                //multiple thread
                int step=Neigh.size()/threadnum;
                boost::thread_group thread;
                for(int i=0;i<threadnum;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==threadnum-1)
                        p.second=Neigh.size();
                    else
                        p.second=(i+1)*step;
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
                }
                thread.join_all();
            }

        }
        else{
            cout<<"x is -1!"<<endl; exit(1);
        }
    }
    tt.stop();
    cout<<"Time for CHWP construction: "<<tt.GetRuntime()<<" s."<<endl;

    IndexsizeCHWP();
}

int	Graph::QueryCHWP(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//function for computing the index size
void Graph::IndexsizeCHWP(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<NeighborCon.size();i++){
        m1+=NeighborCon[i].size()*2*sizeof(int);//dis
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"CH label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

/// H2H algorithm
void Graph::H2HIndexConstruct() {
    string orderfile=graph_path+".order";
    orderfile=graph_path+".vertex_orderMDE2";
    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    tt.start();
    makeTree();
    tt.stop();
    runT2=tt.GetRuntime();
    cout<<"Time for Tree construction: "<<runT2<<" s."<<endl;

    tt.start();
    makeIndex();
    tt.stop();
    runT3=tt.GetRuntime();
    cout<<"Time for Index building: "<<runT3<<" s."<<endl;

    cout<<"Overall index construction time: "<<runT1+runT2+runT3<<" s."<<endl;

    IndexsizeH2H();
}
//function for MDE contraction
void Graph::MDEContraction(string orderfile){
    cout<<"MDE contraction..."<<endl;
    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    ifstream IF(orderfile);
//    if(true){
    if(!IF.is_open()){/// if no order file, use MDE to generate order
        cout<<"Cannot open vertex ordering file "<<orderfile<<endl;
        int Twidth=0;//tree width
        //initialize SCconNodesMT
        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
        }

        _DD_.assign(node_num,0);
        DD.assign(node_num,0);

        set<DegComp> Deg;
        int degree;
        for(int i=0;i<node_num;i++){
            degree=Neighbor[i].size();
            if(degree!=0){
                _DD_[i]=degree;
                DD[i]=degree;
                Deg.insert(DegComp(i));
            }
        }

        vector<bool> exist; exist.assign(node_num,true);
        vector<bool> change; change.assign(node_num,false);

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect); //NeighborCon.clear();
        //SCconNodes.clear();

        //cout<<"Begin to contract"<<endl;
        int count=0;

        while(!Deg.empty()){
            if(count%10000==0)
                cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
            count+=1;
            int x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            vNodeOrder.push_back(x);
            Deg.erase(Deg.begin());
            exist[x]=false;

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }

            if(Neigh.size()>Twidth)
                Twidth=Neigh.size();

            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            //multi threads for n^2 combination
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteEOrderGenerate(x,y);
                change[y]=true;
            }

            int stepf=Neigh.size()/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*stepf;
                if(i==threadnum-1)
                    p.second=Neigh.size();
                else
                    p.second=(i+1)*stepf;
                threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p));
            }
            threadf.join_all();
        }

        NodeOrder.assign(node_num,-1);
        for(int k=0;k<vNodeOrder.size();k++){
            NodeOrder[vNodeOrder[k]]=k;
        }
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << " " << NodeOrder[i] << endl;
        ofile.close();
        cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
        exit(0);
    }
    else{///if there is an order file
        cout<<"Reading vertex ordering... "<<orderfile<<endl;
        NodeOrder.assign(node_num, -1);
        vNodeOrder.assign(node_num, -1);
        int num, nodeID, nodeorder;
        string line;
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        num=stoi(vs[0]);
        for(int i=0;i<num;i++){
            vs.clear();
            getline(IF,line);
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            if(vs.size()!=2){
                cout<<"Wrong syntax for ordering. "<<line<<endl; exit(1);
            }
            nodeID=stoi(vs[0]);nodeorder=stoi(vs[1]);

            NodeOrder[nodeID]=nodeorder;
            if(nodeorder!=-1){
                vNodeOrder[nodeorder]=nodeID;
            }else{
                cout<<"Wrong order! "<<nodeID<<" "<<nodeorder<<endl; exit(1);
            }
        }
        IF.close();
        unordered_set<int> vertices; vertices.clear();
        for(int i=0;i<node_num;++i){
            if(vertices.find(vNodeOrder[i])==vertices.end()){//if not found
                vertices.insert(vNodeOrder[i]);
            }
        }
        if(vertices.size()!=node_num){
            cout<<"Order wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
        }
        vertices.clear();

        for(int i=0;i<2;++i){
            int id=vNodeOrder[node_num-1-i];
            cout<<"Order "<<node_num-1-i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }
        for(int i=0;i<2;++i){
            int id=vNodeOrder[i];
            cout<<"Order "<<i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect);

        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
        }

        vector<bool> exist; exist.assign(node_num,true);
        //vector<bool> change; change.assign(nodenum,false);

        //cout<<"Begin to contract"<<endl;
        for(int nodeorder=0;nodeorder<node_num;nodeorder++){//start from the most important vertex
            int x=vNodeOrder[nodeorder];
            if(x!=-1){//to identify and exclude the isolated vertices
                exist[x]=false;

                vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

                for(auto it=E[x].begin();it!=E[x].end();it++){
                    if(exist[(*it).first]){
                        Neigh.push_back(*it);
                    }
                }
                NeighborCon[x].assign(Neigh.begin(),Neigh.end());

                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteEorder(x,y);
                    //change[y]=true;
                }

                int ID1,ID2;
                if(Neigh.size()<=100){
//                if(true){
                    //single thread
                    for(int i=0;i<Neigh.size();i++){
                        ID1=Neigh[i].first;
                        for(int j=i+1;j<Neigh.size();j++){
                            ID2=Neigh[j].first;
                            insertEorder(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                                }
                                SCconNodesMT[ID1][ID2].push_back(x);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<int>()});
                                }
                                SCconNodesMT[ID2][ID1].push_back(x);
                            }

                        }
                    }
                }else{
                    if(Neigh.size()>threadnum){
                        int step=Neigh.size()/threadnum;
                        boost::thread_group thread;
                        for(int i=0;i<threadnum;i++){
                            pair<int,int> p;
                            p.first=i*step;
                            if(i==threadnum-1)
                                p.second=Neigh.size();
                            else
                                p.second=(i+1)*step;
                            thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }

                }

            }
            else{
                cout<<"Wrong order! "<<x<<" "<<nodeorder<<endl; exit(1);
            }
        }
    }
    NodeOrder_ = NodeOrder;
}
//function for H2H tree construction
void Graph::makeTree(){
    cout<<"Building H2H tree..."<<endl;
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);//record the tree node id whose unique vertex involves this vertex as neighbor

    rank.assign(node_num,0);
    //Tree.clear();
    int len=vNodeOrder.size()-1;
    heightMax=0;

    Node rootn;
    int x=vNodeOrder[len];//from the highest vertex
    //cout<<"len "<<len<<" , ID "<<x<<endl;
    while(x==-1){//to skip those vertices whose ID is -1
        len--;
        x=vNodeOrder[len];
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
    rootn.vert=NeighborCon[x];
    rootn.uniqueVertex=x;
    rootn.pa=-1;
    rootn.height=1;
    rank[x]=0;
    Tree.push_back(rootn);
    len--;

    int nn;
    for(;len>=0;len--){
        int x=vNodeOrder[len];
        Node nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=match(x,NeighborCon[x]);
        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        nod.hdepth=Tree[pa].height+1;

        for(int i=0;i<NeighborCon[x].size();i++){
            nn=NeighborCon[x][i].first;
            VidtoTNid[nn].push_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax) {
            heightMax=nod.height;
        }
        rank[x]=Tree.size();

        Tree.push_back(nod);
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
}
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
    makeRMQ();

    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);
    Tree[0].vAncestor=list;

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}
/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
//            vSm[ID1]->wait();
//            vSm[ID2]->wait();
            insertEMTOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
//            vSm[ID1]->notify();
//            vSm[ID2]->notify();
        }
    }
//    sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
}

void Graph::deleteEorder(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        //DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}

void Graph::insertEorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        //DD[u]++;
        //DD2[u]++;
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        //DD[v]++;
        //DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}

void Graph::insertEMTorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                }
                SCconNodesMT[ID1][ID2].push_back(x);
            }

        }
    }
//    sm->notify();
}
/// Functions for Tree contraction
int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(rank[vert[i].first]>rank[nearest])
            nearest=vert[i].first;
    }
    int p=rank[nearest];
    return p;
}

/// Functions for Tree index contraction
void Graph::makeRMQ(){
    //EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}
void Graph::makeRMQDFS(int p, int height){
    toRMQ[p] = EulerSeq.size();
    EulerSeq.push_back(p);
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Graph::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){//get the original distance information by shortcuts
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();
//    Tree[p].dis.push_back(0);//distance to itself
    Tree[p].vAncestor=list;
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);//the last vertex is the tree node
    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the k-th ancestor is x, the i-th neighbor is also the k-th ancesotr

        for(int j=0;j<list.size();j++){
            int y=list[j];//the j-th ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//if x is the ancestor of y
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//if y is the ancestor of x
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFS(Tree[p].ch[i],list);
    }
    list.pop_back();
}

//function for computing the index size
void Graph::IndexsizeH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m3+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count

    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"H2H label size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"H2H Update information size: "<<(double)(m2+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}


///////////////////////////// Query Processing ////////////////////////////////

//function for efficiency test
void Graph::EffiCheckH2H(string filename,int runtimes){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.push_back(make_pair(ID1, ID2));
    }
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    double runT=0;
    int d1, d2;
    Timer tt;
    clock_t start = clock();

    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
        tt.start();
        d2=QueryH2H(s,t);
        tt.stop();
        runT += tt.GetRuntime();
        results[i]=d2;
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
    }

//    cout<<endl;

    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}

int Graph::QueryH2H(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);
//    cout<<r1<<" "<<r2<<" "<<LCA<<" "<<Tree.size()<<endl;
    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}

int Graph::QueryH2HDebug(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);
    cout<<"r1: "<<ID1<<"("<<r1<<") ; r2: "<<ID2<<"("<<r2<<")"<<endl;
    int hub,dis1,dis2;
    if(LCA==r1){
        cout<<"LCA: r1"<<endl;
        return Tree[r2].dis[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        cout<<"LCA: r2"<<endl;
        return Tree[r1].dis[Tree[r2].pos.back()];
    }
    else{
        cout<<"LCA: "<<Tree[LCA].uniqueVertex<<"("<<NodeOrder[Tree[LCA].uniqueVertex]<<")"<<endl;
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                hub=Tree[r1].vAncestor[Tree[LCA].pos[i]];
                dis1=Tree[r1].dis[Tree[LCA].pos[i]], dis2=Tree[r2].dis[Tree[LCA].pos[i]];
                tmp=dis1+dis2;
            }

        }
        cout<<"hub: "<<hub<<"("<<NodeOrder[hub]<<"); d1: "<<dis1<<"("<<Dijkstra(ID1,hub)<<") ; d2: "<<dis2<<"("<<Dijkstra(hub,ID2)<<") "<<endl;

        return tmp;
    }
}


int Graph::LCAQuery(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

/*-------------Index Maintenance Algorithms--------------*/
/// CHWP algorithm
void Graph::IndexMaintenanceCHWP(int updateType, bool ifBatch, int batchNum, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = graph_path + ".update";
    bool ifDebug=false;
//    ifDebug=true;

    int ID1, ID2, oldW, newW;
    srand (0);
    vector<pair<pair<int,int>,int>> updateData;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    ReadUpdate(file, updateData);
    if(batchNum>updateData.size()/batchSize){
        batchNum=updateData.size()/batchSize;
    }
    cout<<"Update batch number: "<<batchNum<<" ; Batch size: "<<batchSize<<endl;

    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;

            if(ifBatch){//for batch update
                Graph g2=*this;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

                    tt.start();
                    g2.CHdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease batch update Time: "<<runT1/batchNum<<" s; each edge update time: "<<runT1/(batchNum*batchSize)<<" s.\n"<<endl;

            }
            else{//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    CHdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease single-edge update Time: "<<runT1/batchNum<<" s.\n"<<endl;
            }

//            break;
        }
        case 1:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

                    tt.start();
                    CHincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase batch update Time: "<<runT2/batchNum<<" s; each edge update time: "<<runT1/(batchNum*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    CHincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Increase single-edge update Time: "<<runT2/batchNum<<" s.\n"<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::IndexMaintenanceCHWPRealUpdate(string updateFile) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = graph_path + ".update";
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);
    vector<pair<pair<int,int>,int>> updateData;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    string line;
    vector<string> vs;
    int weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==2);
    int batchNum=stoi(vs[0]);
    int batchInterval=stoi(vs[1]);
    int batchSize;
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    long long int aveBatchSize=0;
    int maxBatchSize=0, minBatchSize=INT32_MAX;
    for(int i=0;i<batchNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "));
        batchSize=stoi(vs[0]);
        assert(vs.size()==3*batchSize+1);
        aveBatchSize+=batchSize;
        if(maxBatchSize<batchSize) maxBatchSize=batchSize;
        if(minBatchSize>batchSize) minBatchSize=batchSize;
        for(int j=0;j<batchSize;++j){
            ID1=stoi(vs[3*j+1]), ID2=stoi(vs[3*j+2]), weight=stoi(vs[3*j+3]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();

    Timer tt;
    double runT1=0, runT2 = 0;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    for(int u=0;u<batchNum;u++){
//    for(int u=2;u<4;u++){
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[u].size();++j){
            ID1=batchUpdates[u][j].first.first, ID2=batchUpdates[u][j].first.second, weight=batchUpdates[u][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                if(it->first==ID2){
                    ifFind=true;
                    oldW=it->second;
                    if(oldW>weight){
                        wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }else if(it->second<weight){
                        wBatchInc.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }
                    break;
                }
            }
            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<u<<" . Batch size: "<<wBatchDec.size()+wBatchInc.size()<<"("<<batchUpdates[u].size()<<","<<uEdges.size()<<") ; Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

//        int id1=207116, id2=148496;
//        cout<<"Before update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<QueryCHWP(id1,id2)<<"("<<Dijkstra(id1,id2)<<")"<<endl;
//        for(auto it=NeighborCon[id1].begin();it!=NeighborCon[id1].end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }
//        id1=38900,id2=38901;
//        for(auto it=NeighborCon[id1].begin();it!=NeighborCon[id1].end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }
//        id1=38900,id2=64221;
//        for(auto it=NeighborCon[id1].begin();it!=NeighborCon[id1].end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }

        tt.start();
        if(!wBatchDec.empty()) {
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            CHdecBat(wBatchDec);;
        }

//        if(ifDebug){
//            CorrectnessCheck(100);
//        }
//        cout<<"after decrease update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<QueryCHWP(id1,id2)<<"("<<Dijkstra(id1,id2)<<")"<<endl;

        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
////            vector<pair<pair<int,int>,pair<int,int>>> wBatch;
////            for(int j=0;j<wBatchInc.size();++j){
////                wBatch.clear();
////                wBatch.emplace_back(wBatchInc[j]);
////                cout<<"Update "<<j<<": "<<wBatch[0].first.first<<"("<<NodeOrder[wBatch[0].first.first]<<") "<<wBatch[0].first.second<<"("<<NodeOrder[wBatch[0].first.second]<<") "<<" "<<wBatch[0].second.first<<" "<<wBatch[0].second.second<<endl;
////                H2HincBatMT(wBatch);
////                CorrectnessCheck(100);
////            }
            CHincBatMT(wBatchInc);

            if(u==0){
//                ofstream OF(graph_path+".CH4");
//                if(!OF.is_open()){
//                    cout<<"Fail to open file"<<endl; exit(1);
//                }
//                OF<<node_num<<" "<<edge_num<<endl;
//                for(int i=0;i<Neighbor.size();++i){
//                    for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//                        OF<<i<<" "<<it->first<<" "<<it->second<<endl;
//                    }
//                }
//                OF.close();
//                CHIndexStore(graph_path+".CHIndex4");
//                cout<<"Write done."<<endl;
            }
//            cout<<"after increase update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<QueryCHWP(id1,id2)<<"("<<Dijkstra(id1,id2)<<")"<<endl;

        }
        tt.stop();
        cout<<"update time: "<<tt.GetRuntime()<<" s."<<endl;
        runT1 += tt.GetRuntime();
        if(ifDebug){
            CorrectnessCheck(100);
        }
    }

    cout<<"Average Increase batch update Time: "<<runT1/batchNum<<" s"<<endl;

}

void Graph::CHIndexStore(string filename){
    ofstream OF(filename);
    if(!OF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
    cout<<"Writing labels..."<<endl;
    OF<<NeighborCon.size()<<endl;
    for(int i=0;i<NeighborCon.size();++i){
        OF<<NeighborCon[i].size();
        for(int j=0;j!=NeighborCon[i].size();++j){
            OF<<" "<<NeighborCon[i][j].first<<" "<<NeighborCon[i][j].second.first<<" "<<NeighborCon[i][j].second.second;
        }
        OF<<endl;
    }
    OF.close();
}

void Graph::CHIndexCompare(string filename){
    vector<map<int,pair<int,int>>> CHIndex;
    ifstream IF(filename);
    if(!IF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
//    ofstream OF(filename+".result");
//    if(!OF.is_open()){
//        cout<<"Cannot open file. "<<filename+".result"<<endl; exit(1);
//    }
    cout<<"Index file: "<<filename<<endl;
    cout<<"Read and compare..."<<endl;
    string line;
    vector<string> vs;
    getline(IF,line);
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int vNum=stoi(vs[0]);
    if(vNum!=node_num) {
        cout<<"Wrong node number. "<<vNum<<" "<<node_num<<endl; exit(1);
    }
    CHIndex.assign(vNum,map<int,pair<int,int>>());
    int ID1,ID2,dis,num,cnt;
    for(int i=0;i<vNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        num=stoi(vs[0]);
        ID1=i;
        assert(vs.size()==3*num+1);
        for(int j=0;j<num;++j){
            ID2=stoi(vs[3*j+1]), dis=stoi(vs[3*j+2]), cnt=stoi(vs[3*j+3]);
            CHIndex[ID1].insert({ID2,make_pair(dis,cnt)});
        }
    }
    IF.close();
    for(int i=0;i<NeighborCon.size();++i){
        ID1=i;
        for(auto it=NeighborCon[i].begin();it!=NeighborCon[i].end();++it){
            ID2=it->first, dis=it->second.first, cnt=it->second.second;
            if(CHIndex[ID1].find(ID2)!=CHIndex[ID1].end()){//if found
                if(CHIndex[ID1][ID2].first!=dis){
                    cout<<"Wrong label. "<<ID1<<" "<<ID2<<" "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
//                    exit(1);
                }else if(CHIndex[ID1][ID2].second!=cnt ){
                    cout<<"Wrong cnt. "<<ID1<<" "<<ID2<<" "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
                }
            }else{
                cout<<"Wrong. Not found. "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
            }
        }
    }
}

void Graph::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //maintain the index caused by the weight change
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderComp> OC;
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
    //OC.clear(); OCdis.clear();


    int a,b,newW;//the weight of (a,b) decrease to newW
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first;
        b=wBatch[k].first.second;
        newW=wBatch[k].second.second;

        //modify the information in original graph
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }
                else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[b][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }

        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first>newW){
                        //cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
                        NeighborCon[a][i].second.first=newW;
                        NeighborCon[a][i].second.second=1;

                        OCdis[make_pair(a,b)]=newW;
                        OC.insert(OrderComp(a,b));
                    }else if(NeighborCon[a][i].second.first==newW)
                        NeighborCon[a][i].second.second+=1;
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first>newW){
                        NeighborCon[b][i].second.first=newW;
                        NeighborCon[b][i].second.second=1;

                        OCdis[make_pair(b,a)]=newW;
                        OC.insert(OrderComp(b,a));
                    }else if(NeighborCon[b][i].second.first==newW)
                        NeighborCon[b][i].second.second+=1;
                    break;
                }
            }
        }
    }


    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];
        map<int,int> InM2t; //InM2t.clear();
        vector<pair<int,int>> InMLower; //InMLower.clear();
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
                InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
            else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
                InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
        }
        int inID,inW,inWt;
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(InM2t.find(inID)!=InM2t.end()){
                inW=InM2t[inID];
                inWt=NeighborCon[t][i].second.first;
                if(inWt>inW+wt){
                    NeighborCon[t][i].second.first=inW+wt;
                    NeighborCon[t][i].second.second=1;
//                    if(t==207116 && inID==148496){
//                        cout<<"Flag 1.1 "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<inID<<"("<<NodeOrder[inID]<<") "<<inWt<<" "<<inW<<" "<<wt<<" "<<NeighborCon[t][i].second.second<<endl;
//                    }
                    OCdis[make_pair(t,inID)]=inW+wt;
                    OrderComp oc={t,inID};
                    OC.insert(oc);
                }else if(inWt==inW+wt){
//                    if(t==207116 && inID==148496){
//                        cout<<"Flag 1.2 "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<inID<<"("<<NodeOrder[inID]<<") "<<inWt<<" "<<inW<<" "<<wt<<" "<<NeighborCon[t][i].second.second<<endl;
//                    }
                    NeighborCon[t][i].second.second+=1;
                }
            }
        }

        for(int i=0;i<InMLower.size();i++){
            inID=InMLower[i].first; inW=InMLower[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    inWt=NeighborCon[inID][j].second.first;
                    if(inWt>inW+wt){
                        NeighborCon[inID][j].second.first=inW+wt;
                        NeighborCon[inID][j].second.second=1;
//                        if(inID==207116 && t==148496){
//                            cout<<"Flag 2.1 "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<inID<<"("<<NodeOrder[inID]<<") "<<inWt<<" "<<inW<<" "<<wt<<" "<<NeighborCon[t][i].second.second<<endl;
//                        }
                        OCdis[make_pair(inID,t)]=inW+wt;
                        OrderComp oc={inID,t};
                        OC.insert(oc);
                    }else if(inWt==inW+wt){
//                        if(inID==207116 && t==148496){
//                            cout<<"Flag 2.2 "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<inID<<"("<<NodeOrder[inID]<<") "<<inWt<<" "<<inW<<" "<<wt<<" "<<NeighborCon[t][i].second.second<<endl;
//                        }
                        if(OCdis.find(make_pair(s,inID))==OCdis.end()){//if not found
                            NeighborCon[inID][j].second.second+=1;
                        }

                    }

                    break;
                }
            }
        }
    }//finish change index
}

void Graph::CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderComp> OC; //OC.clear();
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int wb=0;wb<wBatch.size();wb++){
        int a=wBatch[wb].first.first;
        int b=wBatch[wb].first.second;
        int oldW=wBatch[wb].second.first;
        int newW=wBatch[wb].second.second;

        //modify the original graph information
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                if(newW>Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                if(newW>Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[b][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }


        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first==oldW){
                        NeighborCon[a][i].second.second-=1;
                        if(NeighborCon[a][i].second.second<1){
                            OrderComp oc={a,b};
                            OC.insert(oc);
                            OCdis[make_pair(a,b)]=oldW;
                        }
                    }
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first==oldW){
                        NeighborCon[b][i].second.second-=1;
                        if(NeighborCon[b][i].second.second<1){
                            OrderComp oc={b,a};
                            OC.insert(oc);
                            OCdis[make_pair(b,a)]=oldW;
                        }
                    }
                    break;
                }
            }
        }
    }

    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];//distance of s--->t before change
        int inID,inW;
        map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
        //HigherIn.clear(); LowerIn.clear();
        //the shortcuts infected by s-->t
        for(int i=0;i<NeighborCon[s].size();i++){
            inID=NeighborCon[s][i].first;
            inW=NeighborCon[s][i].second.first;
            if(NodeOrder[inID]<NodeOrder[t]){
                LowerIn.emplace_back(inID,inW);
            }else if(NodeOrder[inID]>NodeOrder[t]){
                HigherIn.insert(make_pair(inID,inW));
            }
        }
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(HigherIn.find(inID)!=HigherIn.end()){
                inW=HigherIn[inID];
                if(NeighborCon[t][i].second.first==wt+inW){
                    NeighborCon[t][i].second.second-=1;
                    if(NeighborCon[t][i].second.second<1){
                        OrderComp oc={t,inID};
                        OC.insert(oc);
                        OCdis[make_pair(t,inID)]=wt+inW;
                    }
                }
            }
        }
        for(int i=0;i<LowerIn.size();i++){
            inID=LowerIn[i].first; inW=LowerIn[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    if(NeighborCon[inID][j].second.first==inW+wt){
                        if(OCdis.find(make_pair(s,inID))==OCdis.end()) {//if not found, new
                            NeighborCon[inID][j].second.second -= 1;
                            if (NeighborCon[inID][j].second.second < 1) {
                                OrderComp oc = {inID, t};
                                OC.insert(oc);
                                OCdis[make_pair(inID, t)] = wt + inW;
                            }
                        }
                    }
                    break;
                }
            }
        }

        //get the new weight value of s-->t
        wt=INF; int countwt=0;
        for(int i=0;i<Neighbor[s].size();i++){
            if(Neighbor[s][i].first==t){
                wt=Neighbor[s][i].second;//the weight value in the original graph
                countwt=1;
                break;
            }
        }
        int ssw,wtt,wid;
        vector<int> Wnodes; //Wnodes.clear();
        if(s<t){
            //Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
            Wnodes=SCconNodesMT[s][t];
        }else{
            //Wnodes=SCconNodes[make_pair(t,s)];
            Wnodes=SCconNodesMT[t][s];
        }

        for(int i=0;i<Wnodes.size();i++){
            wid=Wnodes[i];
            for(int j=0;j<NeighborCon[wid].size();j++){
                if(NeighborCon[wid][j].first==s){
                    ssw=NeighborCon[wid][j].second.first;
                }
                if(NeighborCon[wid][j].first==t){
                    wtt=NeighborCon[wid][j].second.first;
                }
            }

            if(ssw+wtt<wt){
                wt=ssw+wtt;
                countwt=1;
            }else if(ssw+wtt==wt){
                countwt+=1;
            }
        }

        //refresh the weight value of s--t in the index
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NeighborCon[s][i].first==t){
//                cout<<"Refresh shortcut: "<<s<<" "<<t<<" "<<NeighborCon[s][i].second.first<<" "<<wt<<endl;
                NeighborCon[s][i].second.first=wt;
                NeighborCon[s][i].second.second=countwt;
                break;
            }
        }
    }
}

/// H2H algorithm
void Graph::IndexMaintenanceH2H(int updateType, bool ifBatch, int batchNum, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = graph_path + ".update";
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);
    if(batchNum>updateData.size()/batchSize){
        batchNum=updateData.size()/batchSize;
    }
    cout<<"Update batch number: "<<batchNum<<" ; Batch size: "<<batchSize<<endl;

    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g2=*this;
            if(ifBatch){//for batch update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

                    tt.start();
                    g2.H2HdecBat(wBatch);
//                    H2HdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease batch update Time: "<<runT1/batchNum<<" s; single edge update time: "<<runT1/(batchNum*batchSize)<<" s.\n"<<endl;
            }
            else{//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    g2.H2HdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease single-edge update Time: "<<runT1/batchNum<<" s.\n"<<endl;
            }

//            break;
        }
        case 1:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                        ++update_i;
                    }

                    tt.start();
                    H2HincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    cout<<"Batch "<<u<<", update time: "<<tt.GetRuntime()<<" s."<<endl;
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase batch update Time: "<<runT2/batchNum<<" s; single edge update time: "<<runT2/(batchNum*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    H2HincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Increase single-edge update Time: "<<runT2/batchNum<<" s.\n"<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::IndexMaintenanceH2HRealUpdate(string updateFile) {
    cout<<"Index update test..."<<endl;
    // read updates
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    string line;
    vector<string> vs;
    int weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==2);
    int batchNum=stoi(vs[0]);
    int batchInterval=stoi(vs[1]);
    int batchSize;
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    long long int aveBatchSize=0;
    int maxBatchSize=0, minBatchSize=INT32_MAX;
    for(int i=0;i<batchNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "));
        batchSize=stoi(vs[0]);
        assert(vs.size()==3*batchSize+1);
        aveBatchSize+=batchSize;
        if(maxBatchSize<batchSize) maxBatchSize=batchSize;
        if(minBatchSize>batchSize) minBatchSize=batchSize;
        for(int j=0;j<batchSize;++j){
            ID1=stoi(vs[3*j+1]), ID2=stoi(vs[3*j+2]), weight=stoi(vs[3*j+3]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();

    Timer tt;
    double runT1=0, runT2 = 0;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    for(int u=0;u<batchNum;u++){
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[u].size();++j){
            ID1=batchUpdates[u][j].first.first, ID2=batchUpdates[u][j].first.second, weight=batchUpdates[u][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                if(it->first==ID2){
                    ifFind=true;
                    oldW=it->second;
                    if(oldW>weight){
                        wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }else if(it->second<weight){
                        wBatchInc.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }
                    break;
                }
            }
            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<u<<" . Batch size: "<<wBatchDec.size()+wBatchInc.size()<<"("<<batchUpdates[u].size()<<","<<uEdges.size()<<") ; Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

//        int id1=144586, id2=144587;
//        cout<<"Before update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") ";
//        for(int ii=0;ii<Tree[rank[id1]].vAncestor.size();++ii){
//            if(Tree[rank[id1]].vAncestor[ii]==id2){
//                cout<<Tree[rank[id1]].dis[ii]<<" "<<Tree[rank[id1]].cnt[ii]<<" "<<Tree[rank[id1]].FN[ii]<<endl;
//            }
//        }
        tt.start();
        if(!wBatchDec.empty()) {
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            H2HdecBat(wBatchDec);
        }

//        if(ifDebug){
//            CorrectnessCheck(100);
//        }
//        cout<<"after decrease update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") ";
//        for(int ii=0;ii<Tree[rank[id1]].vAncestor.size();++ii){
//            if(Tree[rank[id1]].vAncestor[ii]==id2){
//                cout<<Tree[rank[id1]].dis[ii]<<" "<<Tree[rank[id1]].cnt[ii]<<" "<<Tree[rank[id1]].FN[ii]<<endl;
//            }
//        }
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
////            vector<pair<pair<int,int>,pair<int,int>>> wBatch;
////            for(int j=0;j<wBatchInc.size();++j){
////                wBatch.clear();
////                wBatch.emplace_back(wBatchInc[j]);
////                cout<<"Update "<<j<<": "<<wBatch[0].first.first<<"("<<NodeOrder[wBatch[0].first.first]<<") "<<wBatch[0].first.second<<"("<<NodeOrder[wBatch[0].first.second]<<") "<<" "<<wBatch[0].second.first<<" "<<wBatch[0].second.second<<endl;
////                H2HincBatMT(wBatch);
////                CorrectnessCheck(100);
////            }
            H2HincBatMT(wBatchInc);

            if(u==2){
//                ofstream OF(graph_path+".new");
//                if(!OF.is_open()){
//                    cout<<"Fail to open file"<<endl; exit(1);
//                }
//                OF<<node_num<<" "<<edge_num<<endl;
//                for(int i=0;i<Neighbor.size();++i){
//                    for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//                        OF<<i<<" "<<it->first<<" "<<it->second<<endl;
//                    }
//                }
//                OF.close();
//                H2HIndexStore(graph_path+".labelTemp3");
//                cout<<"Write done."<<endl;
            }
//            cout<<"after increase update "<<u<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") ";
//            for(int ii=0;ii<Tree[rank[id1]].vAncestor.size();++ii){
//                if(Tree[rank[id1]].vAncestor[ii]==id2){
//                    cout<<Tree[rank[id1]].dis[ii]<<" "<<Tree[rank[id1]].cnt[ii]<<" "<<Tree[rank[id1]].FN[ii]<<endl;
//                }
//            }
        }
        tt.stop();
        cout<<"update time: "<<tt.GetRuntime()<<" s."<<endl;
        runT1 += tt.GetRuntime();
        if(ifDebug){
            CorrectnessCheck(100);
        }
    }

    cout<<"Average Increase batch update Time: "<<runT1/batchNum<<" s"<<endl;



}


void Graph::H2HIndexStore(string filename){
    ofstream OF(filename);
    if(!OF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
    cout<<"Writing labels..."<<endl;
    OF<<Tree.size()<<endl;
    for(int i=0;i<Tree.size();++i){
        OF<<Tree[i].uniqueVertex<<" "<<Tree[i].dis.size();
        for(int j=0;j!=Tree[i].dis.size();++j){
            OF<<" "<<Tree[i].vAncestor[j]<<" "<<Tree[i].dis[j]<<" "<<Tree[i].cnt[j];
        }
        OF<<endl;
    }
    OF.close();
}

void Graph::H2HIndexCompare(string filename){
    vector<map<int,pair<int,int>>> DisLabels;
    ifstream IF(filename);
    if(!IF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
//    ofstream OF(filename+".result");
//    if(!OF.is_open()){
//        cout<<"Cannot open file. "<<filename+".result"<<endl; exit(1);
//    }
    cout<<"Index file: "<<filename<<endl;
    cout<<"Read and compare..."<<endl;
    string line;
    vector<string> vs;
    getline(IF,line);
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int vNum=stoi(vs[0]);
    if(vNum!=node_num) {
        cout<<"Wrong node number. "<<vNum<<" "<<node_num<<endl; exit(1);
    }
    DisLabels.assign(vNum,map<int,pair<int,int>>());
    int ID1,ID2,dis,num,cnt;
    for(int i=0;i<vNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]), num=stoi(vs[1]);
        assert(vs.size()==3*num+2);
        for(int j=0;j<num;++j){
            ID2=stoi(vs[3*j+2]), dis=stoi(vs[3*j+3]), cnt=stoi(vs[3*j+4]);
            DisLabels[ID1].insert({ID2,make_pair(dis,cnt)});
        }
    }
    IF.close();
    for(int i=0;i<Tree.size();++i){
        ID1=Tree[i].uniqueVertex;
        for(int j=0;j<Tree[i].dis.size();++j){
            ID2=Tree[i].vAncestor[j], dis=Tree[i].dis[j], cnt=Tree[i].cnt[j];
            if(DisLabels[ID1].find(ID2)!=DisLabels[ID1].end()){//if found
                if(DisLabels[ID1][ID2].first!=dis){
                    cout<<"Wrong label. "<<ID1<<" "<<ID2<<" "<<DisLabels[ID1][ID2].first<<"("<<dis<<") "<<DisLabels[ID1][ID2].second<<"("<<cnt<<")"<<endl;
//                    exit(1);
                }else if(DisLabels[ID1][ID2].second!=cnt ){
                    cout<<"Wrong cnt. "<<ID1<<" "<<ID2<<" "<<DisLabels[ID1][ID2].first<<"("<<dis<<") "<<DisLabels[ID1][ID2].second<<"("<<cnt<<")"<<endl;
                }
            }else{
                cout<<"Wrong. Not found. "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
            }
        }
    }
}

void Graph::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear(); not truly define the shortcut priority since set<int> rank the shortcut by vertex ID
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,oldW,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Wrong edge update."<<endl; exit(1);
                }

                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }else{
                    cout<<"Wrong edge update."<<endl; exit(1);
                }

                break;
            }
        }

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProID]].dis[cidH]>Cw){
//                if(ProID==144586 && Cid==144587){
//                    cout<<"Here 0 "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].cnt[cidH]<<" "<< Tree[rank[ProID]].dis[cidH]<<" "<<Cw<<endl;
//                }
                Tree[rank[ProID]].dis[cidH]=Cw;
                Tree[rank[ProID]].FN[cidH]=true;
                Tree[rank[ProID]].cnt[cidH]=1;//new
                ProIDdisCha=true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            }else if(Tree[rank[ProID]].dis[cidH]==Cw){
//                if(ProID==144586 && Cid==144587){
//                    cout<<"Here 01 "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].cnt[cidH]<<" "<< Tree[rank[ProID]].dis[cidH]<<" "<<Cw<<endl;
//                }
//                cout<<"Reduced Cw results in cnt++ "<<ProID<<" "<<Cid<<" "<<Cw<<" "<<Tree[rank[ProID]].cnt[cidH]<<endl;
                Tree[rank[ProID]].FN[cidH]=true;
                Tree[rank[ProID]].cnt[cidH]+=1;//new
            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()){//if not found, avoid repeated count
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
    bool ProIDdisCha=false;
    bool flagFind=false;
//    if(Tree[child].uniqueVertex==144586){
//        cout<<"Find "<<Tree[child].uniqueVertex<<endl;
//        flagFind=true;
//    }
    vector<int> cntNew(line.size(),0);
    vector<bool> flagUpdate(line.size(),false);
    if(Tree[child].DisRe.size()!=0){//if ChildID has changed labels
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){//if dis comes from shortcut
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//if this dis and shortcut from childID to b has changed
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
//                        if(flagFind && Tree[child].vAncestor[i]==144587){
//                            cout<<"Here 1 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[b]].dis[i]<<" "<<Tree[child].cnt[i]<<endl;
//                        }
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];//use cntNew to redress the cnt value since the edge decrease may lead to more ways for dis (i.e., increase the cnt)
                            }
                        }
                    }
//                    if(vbW<Tree[child].dis[i])
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
//                        if(flagFind && Tree[child].vAncestor[i]==144587){
//                            cout<<"Here 2 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[line[i]]].dis[bH]<<" "<<Tree[child].cnt[i]<<endl;
//                        }
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }
                else{//if this dis from childID to b has not changed
                    if(vertexIDChL.find(b)!=vertexIDChL.end()){//if b's label has changed
                        for(int i=0;i<bH;i++){
                            if(flagFind && Tree[child].vAncestor[i]==144587){
                                cout<<"Here 3 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[b]].dis[i]<<" "<<Tree[child].cnt[i]<<endl;
                            }
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                Tree[child].cnt[i]=1;//new
                                ProIDdisCha=true;
                                flagUpdate[i]=true;
                                cntNew[i]=1;
                            }
                            else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                                cntNew[i]++;
                                if(flagUpdate[i]) {
                                    Tree[child].cnt[i]+=1;//new
                                }
                                else if(cntNew[i]>Tree[child].cnt[i]){
                                    Tree[child].cnt[i]=cntNew[i];
                                }
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        if(flagFind && Tree[child].vAncestor[i]==144587){
//                            cout<<"Here 4 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[line[i]]].dis[bH]<<" "<<Tree[child].cnt[i]<<endl;
//                        }
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }
            }
        }
    }
    else{//if childID's labels have not changed
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){//if dis from childID to b comes from shortcut
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
//                        if(flagFind && Tree[child].vAncestor[i]==144587){
//                            cout<<"Here 5 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[b]].dis[i]<<" "<<Tree[child].cnt[i]<<endl;
//                        }
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
//                    if(flagFind && Tree[child].vAncestor[i]==144587){
//                        cout<<"Here 6 "<<Tree[child].uniqueVertex<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[line[i]]].dis[bH]<<" "<<Tree[child].cnt[i]<<endl;
//                    }
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        Tree[child].cnt[i]=1;//new
                        ProIDdisCha=true;
                        flagUpdate[i]=true;
                        cntNew[i]=1;
                    }
                    else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                        cntNew[i]++;
                        if(flagUpdate[i]) {
                            Tree[child].cnt[i]+=1;//new
                        }
                        else if(cntNew[i]>Tree[child].cnt[i]){
                            Tree[child].cnt[i]=cntNew[i];
                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}

void Graph::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear(); the affected shortcut pair
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//the lower-order vertex of the affected shortcut, vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){
            for(int i=0;i<Neighbor[a].size();i++){
                if(Neighbor[a][i].first==b){
                    Neighbor[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbor[b].size();i++){
                if(Neighbor[b][i].first==a){
                    Neighbor[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){//the shortcut needs update, should be increased
                            OCdis[make_pair(lid,hid)]=oldW;//original weight of the affected shortcut
                            SCre[lid].insert(hid);//the affected shortcut pair
                            OC.insert(OrderCompp(lid));//the lower-order vertex of the affected shortcut
                        }
                    }
                    break;
                }
            }
        }
        else{
            cout<<"wrong update. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl; exit(1);
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    /// Shortcut update
    while(!OC.empty()){
        ProID=(*OC.begin()).x;//from the lowest-order vertex
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        // get the ancestors of ProID, each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

//            if(ProID==178597 && Cid==40366){
//                cout<<"here 1. "<<ProID<<" "<<Cid<<" "<<Cw<<endl;
//            }

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.emplace_back(Vert[j].first,Vert[j].second.first);
                }
            }
            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid]].vert[k].second.second -= 1;
                                if (Tree[rank[lid]].vert[k].second.second < 1) {
                                    SCre[lid].insert(Cid);
                                    OC.insert(OrderCompp(lid));
                                    OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }




            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int CwNew=INF; int countwt=0;

            for(int i=0;i<Neighbor[ProID].size();i++){
                if(Neighbor[ProID][i].first==Cid){
                    CwNew=Neighbor[ProID][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<int> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i];
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<CwNew){
                        CwNew=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==CwNew){
                        countwt+=1;
                    }
                }
            }
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=CwNew;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }
//            if(ProID==178597 && Cid==40366){
//                cout<<"here 2. "<<ProID<<" "<<Cid<<" "<<Cw<<endl;
//            }

            //cout<<Cw<<endl;
            if(CwNew>Cw){//Cw may not change
                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){//if the distance label is from shortcut, then the label may be affected.
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
//                            if(ProID==182172 && Tree[rank[ProID]].vAncestor[i]==291263){
//                                cout<<"Find 1. "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].vAncestor[i]<<" "<<Tree[rank[ProID]].cnt[i]<<endl;
//                            }
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;//? may still be the source
                    Tree[rank[ProID]].cnt[cidH]-=1;//cnt may not be accurate after update
//                    if(ProID==182172 && Cid==291263){
//                        cout<<"Find 2. "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].cnt[cidH]<<endl;
//                    }

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
//                            if(ProID==182172 &&  Tree[rank[ProID]].vAncestor[i]==291263){
//                                cout<<"Find 3. "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].vAncestor[i]<<" "<<Tree[rank[ProID]].cnt[i]<<endl;
//                            }
                        }
                    }
                }
            }else{
//                cout<<"Shortcut weight do not increase. "<<ProID<<" "<<Cid<<" "<<Cw<<"("<<CwNew<<") "<<countwt<<endl;
//                exit(1);
            }


        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){//if they are in different branches
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;//identify the roots
        }

    }

    int ProBeginVertexID;
//    cout<<"Root number: "<<ProBeginVertexSet.size()<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){//for each root
        ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<i<<" "<<ProBeginVertexID<<" "<<Tree[rank[ProBeginVertexID]].height<<endl;
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum);
    }
    //return checknum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    bool flagFind=false;
    for(int i=0;i<Tree[children].dis.size();i++){
//        if(childID==243096 && Tree[children].vAncestor[i]==200022){
//            cout<<"here 3. "<<childID<<" "<<Tree[children].vAncestor[i]<<" "<<Tree[children].dis[i]<<" "<<Tree[children].cnt[i]<<endl;
//            flagFind=true;
//        }
        if(Tree[children].cnt[i]<=0){//if the distance label to i-th ancestor should be maintained
//        if(true){
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//for the tree node that contains childID as vert element
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){//if label is from shortcut
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){//for the tree node that contains line[i] as vert element
                PID=VidtoTNid[line[i]][k];
//                if(Tree[PID].height>Tree[children].height){///modified for correctness
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){
//                    if(PID>Tree.size()){
//                        cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
//                    }
//                    if(childH>Tree[PID].dis.size()){
//                        cout<<"childH error! "<<childH<<" "<<Tree[PID].dis.size()<<": "<<children<<"("<<Tree[children].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
//                    }
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){///
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            int hubFinal,dis1Final, dis2Final;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                        hubFinal=b, dis1Final=Dvb, dis2Final=Tree[rank[line[i]]].dis[bH];
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                        hubFinal=b, dis1Final=Dvb, dis2Final=0;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                        hubFinal=b, dis1Final=Dvb, dis2Final=Tree[rank[b]].dis[i];
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) {
                Tree[children].FN[i]=true;
            }
//            if(flagFind){
//                cout<<"Find. "<<childID<<" "<<Tree[children].vAncestor[i]<<": "<<hubFinal<<"("<<NodeOrder[hubFinal]<<") "<<dis1Final<<"("<<Dijkstra(childID,b)<<") "<<dis2Final<<"("<<Dijkstra(b,Tree[children].vAncestor[i])<<") "<<Tree[children].dis[i]<<" "<<dis<<endl;
//                flagFind=false;
//            }
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
        flagFind=false;
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}


#endif //H2H_HPP