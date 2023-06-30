//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef INDEXBASED_HPP
#define INDEXBASED_HPP

#include "head.h"

/*-------------Index-based Algorithms--------------*/
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
    int ifWrite=false;
//    ifWrite=true;
    if(ifWrite){
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
        //if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation
        //if(vbVisited[w])
        //continue;
        int d=(*ivp).second+dUV;
        mWDistance[w]=d;//the distance from w to ID1 via ID2
        if(d>maxWDistance)
            maxWDistance=d;
        mDistance[w]=INF;
        //}
    }

    if(mDistance.empty()){
        return 0;
    }

    for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
        vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

    int dThreshold=maxWDistance;//to restrict the search scope

    while(!Heap.empty()){
        Heap.extract_min(topNodeID, topDistance);
        if(vbVisited[topNodeID])
            continue;
        if(topDistance>dThreshold)
            break;
        for(auto ivp=vvNode[topNodeID].begin();ivp!=vvNode[topNodeID].end();ivp++){
            neighborNodeID = (*ivp).first;
            neighborLength = (*ivp).second;
            if(vbVisited[neighborNodeID] || neighborNodeID==ID2)//skip ID2
                continue;
            int d = vDistance[topNodeID] + neighborLength;
            if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
            {
                vDistance[neighborNodeID] = d;
                Heap.update(neighborNodeID, d);
                if(mWDistance.find(neighborNodeID) != mWDistance.end())
                    mDistance[neighborNodeID] = d;
            }else if(vDistance[neighborNodeID]>d){
                vDistance[neighborNodeID] = d;
                Heap.update(neighborNodeID, d);
                if(mWDistance.find(neighborNodeID) != mWDistance.end())
                    mDistance[neighborNodeID] = d;
            }
            else if(vDistance[neighborNodeID] == d)
            {
                if(mWDistance.find(neighborNodeID) != mWDistance.end())//if found
                    mDistance[neighborNodeID] = d;
            }
        }
    }

    for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
    {
        if((*imDistance).second > mWDistance[(*imDistance).first])
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
    vector<pair<int,int>> nodeinfor;
    vvpShortCut.assign(node_num, nodeinfor);
    vector<bool> vbVisited(node_num, false);
    map<int, int> mPosition;
    mPosition[0]=node_num-1;
    bool bUpdated;
    vector<pair<int,int>> vU,vW;

    //cout<<"Building static CH "<<endl;

    int v;// current contracting vertex
    int step=node_num/10000;
    step*=1000;

    for(int i=0;i<vNodeOrder.size()-1;i++){
        if(i%step==0)
            cout<<"Contracting vertex "<<i<<endl;
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
                        if(w<ID1)
                            SupportNodes[w][ID1].push_back(v);
                        else
                            SupportNodes[ID1][w].push_back(v);
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
}
//CH shortcut construction
void Graph::CHConstruction(){
    string indexfile=graph_path+".shortcutCH";
    string orderfile=graph_path+".orderCH";
    MDEOrderGenerate(orderfile);

    vector<pair<int,int>> vecp;
    AdjaShort.assign(node_num, vecp);
//    AdjaShortR.assign(node_num, vecp);
    vvNode.assign(Neighbor.begin(), Neighbor.end());//duplicate of original adjacency list Neighbor
    map<int, vector<int>> vecm;
    SupportNodes.assign(node_num, vecm);

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;

    fstream file;
    file.open(indexfile);
    if(!file){
        t1=std::chrono::high_resolution_clock::now();
        CHConstructionMT();
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Time for CH construction: "<<runT<<" s."<<endl;
        writeShortCutorder(indexfile);
//        ReadShortCut(indexfile);
    }
    else
    {
        ReadShortCut(indexfile);
    }
}

//CH algorithm with pruning
int Graph::CH(int ID1, int ID2){
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

/// PLL
//PLL construction entry
void Graph::PLLConstruction(int strategy){
    string indexfile=graph_path+".PLL";
    string PPRfile=graph_path+".PPR";

    VertexOrdering(strategy);

//    fstream file;
//    file.open(indexfile);
//    if(!file){
        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration<double> time_span;
        double runT;

        t1=std::chrono::high_resolution_clock::now();
        PLLConstruct();
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"PLL construction time: "<<runT<<" s."<<endl;

//        writePLL(indexfile, PPRfile);
//    }
//    else
//    {
//        readPLL(indexfile, PPRfile);
//    }
}
//PLL construction function
void Graph::PLLConstruct(){
    unordered_map<int,int> m; m.clear();
    Label.assign(node_num,m);

    unordered_map<int,vector<int>> unorderm; unorderm.clear();
//    PruningPointNew.assign(node_num,unorderm);


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
        for(int j=0;j<vp.size();j++){
            Label[vp[j].first].insert(make_pair(ID, vp[j].second));
            //cout<<vp[j].first<<" "<<vp[j].second<<endl;
        }
    }
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
        for(unordered_map<int,int>::iterator it=Label[nodeID].begin();it!=Label[nodeID].end();it++){
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
    unordered_map<int,int> M0; M0.clear();
    Label.assign(num,M0);
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
//query processing with current PLL index
int Graph::PLLDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;
    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
        hub=(*it).first;
        dis1=(*it).second;
        if(Label[ID2].find(hub)!=Label[ID2].end()){
            dis2=Label[ID2][hub];
            if(dis1+dis2<d){
                d=dis1+dis2;
//                SupNode.clear();
//                SupNode.push_back(hub);
            }else if(dis1+dis2==d){
//                SupNode.push_back(hub);
            }
        }
    }

    return d;
}
//pruned Dijkstra from nodeID
void Graph::DijksPrune1(int nodeID, vector<pair<int,int>>& vp){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(nodeID,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);

    distance[nodeID]=0;
    int topNodeID, topNodeDis;
    vector<pair<int,int>>::iterator it;
    int NNodeID,NWeigh;

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        closed[topNodeID]=true;

        int TempDis; vector<int> SupNode;
        PLLDisQuery1(nodeID, topNodeID,SupNode,TempDis);
        if(TempDis<=topNodeDis){
            //for index update
//            if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
//                for(int k=0;k<SupNode.size();k++){
//                    int supn=SupNode[k];
//
//                    PruningPointNew[topNodeID][supn].push_back(nodeID);
//                    PruningPointNew[nodeID][supn].push_back(topNodeID);
//
//                }
//            }
            continue;
        }


        //Label[topNodeID].insert(nodeID, topNodeDis);
        vp.push_back(make_pair(topNodeID,topNodeDis));
        for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
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

//function of PLL query
int Graph::PLL(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;
    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
void Graph::PSLConstruction(int strategy){
    string indexfile=graph_path+".PSL";
    string PPRfile=graph_path+"PSL.PPR";

    VertexOrdering(strategy);

//    fstream file;
//    file.open(indexfile);
//    if(!file){
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span;
    double runT;

    t1=std::chrono::high_resolution_clock::now();
    PSLConstruct();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    runT= time_span.count();
    cout<<"PSL construction time: "<<runT<<" s."<<endl;

//        writePLL(indexfile, PPRfile);
//    }
//    else
//    {
//        readPLL(indexfile, PPRfile);
//    }
}

// PSL construction function
void Graph::PSLConstruct(){
    unordered_map<int,int> map0; map0.clear();
    Label.assign(node_num, map0);
    Dhop.assign(node_num, map0);

    DvertexNew.assign(node_num, true);

    for(int i=0;i<node_num;i++){
        Label[i].insert(make_pair(i,0));
        Dhop[i].insert(make_pair(i,0));
        //Dvectex.push_back(i);
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

bool Graph::DhopLableRefreshStep(int step){
    bool flag=false;
    vector<unordered_map<int,int>> newDhop;
    unordered_map<int,int> m0; m0.clear();
    newDhop.assign(node_num,m0);
    vector<int> newDvec;
//    vector<vector<tri>> PP;
//    vector<tri> prp; prp.clear();
//    PP.assign(threadnum, prp);

    boost::thread_group thread;

    //int stepthr=Dvectex.size()/threadnum; //cout<<" Dvectex size "<<Dvectex.size()<<endl;
    vector<vector<int>> ProcessID;
    vector<int> vvv; ProcessID.assign(threadnum,vvv);
    threadDistribute(ProcessID);

    for(int i=0;i<ProcessID.size();i++){
        vector<int> p=ProcessID[i];
        thread.add_thread(new boost::thread(&Graph::labelMultiThread2New, this, boost::ref(newDhop), p,step));
    }
    thread.join_all();

    //cout<<"lllllllllllllllllllllllllllllll"<<endl;

    int zerolabel=0;
    Dhop.assign(newDhop.begin(), newDhop.end());
    Dvectex.clear(); set<int> Dset;
    DvertexNew.assign(node_num,false);
    for(int nodeID=0;nodeID<node_num;nodeID++){
        if(Dhop[nodeID].size()>0){
            flag=true;
            for(unordered_map<int,int>::iterator it=Dhop[nodeID].begin();it!=Dhop[nodeID].end();it++){
                if(Label[nodeID].find((*it).first)!=Dhop[nodeID].end()){
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

void Graph::labelMultiThread2New(vector<unordered_map<int,int>>& newDhop, vector<int>& p,int step){
//	sm->wait();

    for(int i=0;i<p.size();i++){
        int nodeID=p[i];
        unordered_map<int,int> Dhop0; Dhop0.clear();
        int neighID, neighW;
        for(int Index=0;Index<Neighbor[nodeID].size();Index++){//Neighbors
            neighID=Neighbor[nodeID][Index].first;
            //cout<<"neighID "<<neighID<<" Dhop.size "<<Dhop[neighID].size()<<endl;
            if(Dhop[neighID].size()>0){
                neighW=Neighbor[nodeID][Index].second;

                unordered_map<int,int>::iterator it=Dhop[neighID].begin();
                int hub, dis, d;
                for(;it!=Dhop[neighID].end();it++){
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
//                        else{
//                            for(int k=0;k<SupNode.size();k++){
//                                int supn=SupNode[k];
//
//                                if(supn!=nodeID && supn!=hub){
//                                    vSm[nodeID]->wait();
//                                    PruningPointNew[nodeID][supn].push_back(hub);
////									PruningPointStepNew[nodeID][supn].insert(make_pair(hub,step));
//                                    PruningPointList[nodeID][supn].push_back(hub);
//
//                                    vSm[nodeID]->notify();
//
//                                    vSm[hub]->wait();
//                                    PruningPointNew[hub][supn].push_back(nodeID);
////									PruningPointStepNew[hub][supn].insert(make_pair(nodeID,step));
//                                    PruningPointList[hub][supn].push_back(nodeID);
//
//                                    vSm[hub]->notify();
//                                }
//                            }
//                        }

                    }
                }
            }
        }
        newDhop[nodeID]=Dhop0;
    }

    //cout<<"one thread finish running!"<<endl;

//	sm->notify();
}
//query by current labels
int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
    d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;

    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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
int Graph::ShortestDisQuery2(int ID1,int ID2){
    int d=INF;
    unordered_map<int,int>::iterator it;
    int hub, dis1, dis2;

    for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
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

#endif //INDEXBASED_HPP