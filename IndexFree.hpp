//
// Created by Xinjie ZHOU on 24/05/2023.
//
#ifndef INDEXFREE_HPP
#define INDEXFREE_HPP

#include "head.h"

/*-------------Index-free Algorithms--------------*/
//Dijkstra's algorithm
int Graph::Dijkstra(int ID1, int ID2){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<int> pre(node_num, 0);
    distance[ID1]=0;
    int topNodeID, topNodeDis;

    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(topNodeID==ID2){
            d=distance[ID2];
//            cout<<ID1<<" "<<ID2<<": "<<d<<endl;
            break;
        }
        closed[topNodeID]=true;

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

    return d;
}

//retrieve all shortest paths
int Graph::DijkstraPaths(int ID1, int ID2){
    benchmark::heap<2, int, int> pqueue(node_num);
    pqueue.update(ID1,0);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);
    vector<unordered_set<int>> pres(node_num, unordered_set<int>());
    distance[ID1]=0;
    pres[ID1].insert(ID1);
    int topNodeID, topNodeDis;

    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(topNodeID==ID2){
            d=distance[ID2];
//            cout<<ID1<<" "<<ID2<<": "<<d<<endl;
            break;
        }
        closed[topNodeID]=true;

        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+topNodeDis;
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh);
                    pres[NNodeID].clear();
                    pres[NNodeID].insert(topNodeID);
                }else if(distance[NNodeID]==NWeigh){
                    pres[NNodeID].insert(topNodeID);
                }
            }
        }
    }
    RetrievePaths(ID1,ID2,pres);
    return d;
}

void Graph::DijkstraV(vector<pair<int,int>>& p){
    /// boost-based implementation
    boost::thread_group thread;
    for(int i=0;i<p.size();i++){
        thread.add_thread(new boost::thread(&Graph::Dijkstra, this, p[i].first, p[i].second));
    }
    thread.join_all();
}


//Bidirectional Dijkstra's algorithm
int Graph::BiDijkstra(int ID1, int ID2){
    if(ID1==ID2) return 0;
    //if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;//to avoid the incorrectness caused by the isolated vertex
    benchmark::heap<2, int, int> queueF(node_num), queueB(node_num);
    queueF.update(ID1,0);
    queueB.update(ID2,0);

    vector<bool> closedF(node_num, false), closedB(node_num, false);
    vector<int> distanceF(node_num, INF), distanceB(node_num, INF);

    distanceF[ID1]=0;
    distanceB[ID2]=0;
    int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
    int NNodeIDF,NWeighF, NNodeIDB, NWeighB;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!queueF.empty() || !queueB.empty()){
//        if(queueF.top_key()+queueB.top_key()>=d){//termination condition 1
//            return d;
//        }
        //forward
        queueF.extract_min(topNodeIDF, topNodeDisF);
        closedF[topNodeIDF]=true;
        if(closedB[topNodeIDF]){//termination condition 2
            return d;
        }
        for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
            NNodeIDF=(*itF).first;
            NWeighF=(*itF).second+topNodeDisF;
            if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
                d=NWeighF+distanceB[NNodeIDF];
            }
            if(!closedF[NNodeIDF]){
                if(distanceF[NNodeIDF]>NWeighF){
                    distanceF[NNodeIDF]=NWeighF;
                    queueF.update(NNodeIDF, NWeighF);
                }
            }
        }
        //backward
        queueB.extract_min(topNodeIDB, topNodeDisB);
        closedB[topNodeIDB]=true;
        if(closedF[topNodeIDB]){//termination condition 2
            return d;
        }
        for(auto itB=Neighbor[topNodeIDB].begin();itB!=Neighbor[topNodeIDB].end();itB++){
            NNodeIDB=(*itB).first;
            NWeighB=(*itB).second+topNodeDisB;
            if(closedF[NNodeIDB] && NWeighB+distanceF[NNodeIDB]<d){
                d=NWeighB+distanceF[NNodeIDB];
            }
            if(!closedB[NNodeIDB]){
                if(distanceB[NNodeIDB]>NWeighB){
                    distanceB[NNodeIDB]=NWeighB;
                    queueB.update(NNodeIDB, NWeighB);
                }
            }
        }
    }
    return d;
}
//A* algorithm
int Graph::Astar(int ID1, int ID2){
    if(ID1==ID2) return 0;
    benchmark::heap<2, int, int> pqueue(node_num);
    int heurisDis=EuclideanDis(ID1,ID2);
    pqueue.update(ID1,heurisDis);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);//actual distance from source

    distance[ID1]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(topNodeID==ID2){
            d=distance[ID2];
            break;
        }
        closed[topNodeID]=true;

        for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
            NNodeID=(*it).first;
            NWeigh=(*it).second+distance[topNodeID];
            heurisDis=EuclideanDis2(NNodeID,ID2);
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh+heurisDis);
                }
            }
        }
    }
    return d;
}
//Bidirectional A* algorithm, consistent approach, correct
int Graph::BiAstar(int ID1, int ID2){
    if(ID1==ID2) return 0;
    //if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;//to avoid the incorrectness caused by the isolated vertex
    benchmark::heap<2, int, int> queueF(node_num), queueB(node_num);
    queueF.update(ID1,0);
    queueB.update(ID2,0);

    vector<bool> closedF(node_num, false), closedB(node_num, false);
    vector<int> distanceF(node_num, INF), distanceB(node_num, INF);

    distanceF[ID1]=0;
    distanceB[ID2]=0;
    int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
    int NNodeIDF,NWeighF, NNodeIDB, NWeighB;

    int d=INF;//initialize d to infinite for the unreachable case
    int potentialF, potentialB;

    while(!queueF.empty() || !queueB.empty()){
        if(queueF.top_key()+queueB.top_key()>=d){
            return d;
        }
        //forward
        queueF.extract_min(topNodeIDF, topNodeDisF);
        closedF[topNodeIDF]=true;
        for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
            NNodeIDF=(*itF).first;
            NWeighF=(*itF).second+distanceF[topNodeIDF];
            potentialF=0.5*(EuclideanDis2(ID1,NNodeIDF) - EuclideanDis2(ID2,NNodeIDF));
            if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
                d=NWeighF+distanceB[NNodeIDF];
            }
            if(!closedF[NNodeIDF]){
                if(distanceF[NNodeIDF]>NWeighF){
                    distanceF[NNodeIDF]=NWeighF;
                    queueF.update(NNodeIDF, NWeighF+potentialF);
                }
            }
        }
        //backward
        queueB.extract_min(topNodeIDB, topNodeDisB);
        closedB[topNodeIDB]=true;
        for(auto itB=Neighbor[topNodeIDB].begin();itB!=Neighbor[topNodeIDB].end();itB++){
            NNodeIDB=(*itB).first;
            NWeighB=(*itB).second+distanceB[topNodeIDB];
            potentialB=0.5*(EuclideanDis2(ID2,NNodeIDB) - EuclideanDis2(ID1,NNodeIDB));
            if(closedF[NNodeIDB] && NWeighB+distanceF[NNodeIDB]<d){
                d=NWeighB+distanceF[NNodeIDB];
            }
            if(!closedB[NNodeIDB]){
                if(distanceB[NNodeIDB]>NWeighB){
                    distanceB[NNodeIDB]=NWeighB;
                    queueB.update(NNodeIDB, NWeighB+potentialB);
                }
            }
        }
    }
    return d;
}
//Bellman-Ford algorithm
int Graph::BellmanFord(int ID1, int ID2) {
    vector<int> distance(node_num, INF);
    vector<int> pre(node_num, -1);

    distance[ID1]=0;
    int topNodeID, topNodeDis;
    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    //relaxations
    for(int i=1;i<=node_num-1;++i){
        for(auto it=Edges.begin();it!=Edges.end();++it){
            topNodeID=it->s; NNodeID=it->t;
            if(distance[topNodeID]!=INF){
                NWeigh=distance[topNodeID]+it->weight;
                if(NWeigh < distance[NNodeID]){
                    distance[NNodeID] = NWeigh;
                    pre[NNodeID]=topNodeID;
//                    cout<<topNodeID<<" "<<NNodeID<<" "<<distance[topNodeID]<<" "<<distance[NNodeID]<<endl;

                }
            }
        }

    }
    //judge negative cycle
    for(int i=0;i<node_num;++i){
        int preID=pre[i];
    }
    d=distance[ID2];

    for(auto it=Edges.begin();it!=Edges.end();++it) {
        topNodeID=it->s; NNodeID=it->t;
        if (distance[topNodeID] != INF) {
            NWeigh = distance[topNodeID] + it->weight;
            if (NWeigh < distance[NNodeID]) {
                cout << "The graph contains a negative cycle!" << endl;
                return d;
            }
        }
    }

    return d;
}

//function of computing Euclidean distance
int Graph::EuclideanDis(int s, int t)
{
    double d=0.111319;//distance in equator per degree
    double coe;
    double temp=(Coordinate[s].second+Coordinate[t].second)/2;
    temp=temp*3.1415926/(180*1000000);
//    temp=temp*3.1415926/(180*1000000);
    coe=cos(temp);
    double y1 = Coordinate[s].second * d, y2 = Coordinate[t].second * d;
    double x1 = Coordinate[s].first * d * coe, x2 = Coordinate[t].first * d * coe;
    double xx = x1 - x2, yy = y1 - y2;
    return sqrt(xx * xx + yy * yy);//Euclidean Distance
}
//function of computing Euclidean distance, using conversed function
int Graph::EuclideanDis2(int s, int t) {
    //return abs(start.x - end.x) + abs(start.y - end.y);//Manhattan Distance
    double xx = Coordinate2[s].first - Coordinate2[t].first, yy = Coordinate2[s].second - Coordinate2[t].second;
    return sqrt(xx*xx + yy*yy);//Euclidean Distance
    //return sqrt(pow(start.x - end.x,2) + pow(start.y - end.y, 2));//Euclidean Distance
}
//estimated function
int Graph::EuclideanDis3(int s, int t) {
    int lat=(int)(abs(Coordinate[s].first-Coordinate[t].first)*111319);
    int lon=(int)(abs(Coordinate[s].second-Coordinate[t].second)*83907);
    int min,max;
    min=(lat>lon)?lon:lat;
    max=(lat>lon)?lat:lon;
    int approx=max*1007+min*441;
    if(max<(min<<4))
        approx-=max*40;
    return (approx+512)>>10;

}

//function of retrieving the shortest path
void Graph::RetrievePath(int ID1, int ID2, vector<int> & prece){

    //path retrieval
    vector<int> path;
    path.clear();
    path.push_back(ID2);
    int preID=prece[ID2];
    while(preID!=ID1){
        path.push_back(preID);
        preID=prece[preID];
    }
    path.push_back(ID1);

    pair<int,int> highestVertex(-1,0);//ID, order
    cout<<"path from "<<ID1<<" to "<<ID2<<": "<<endl;
    for(int i=path.size()-1;i>-1;i--){
        cout<<" "<<path[i]<<"("<<NodeOrder[path[i]]<<") ";//<<endl;
        if(NodeOrder[path[i]] > highestVertex.second){
            highestVertex.second = NodeOrder[path[i]];
            highestVertex.first = path[i];
        }
        if(i>0){
            for(int j=0;j<Neighbor[path[i]].size();++j){
                if(Neighbor[path[i]][j].first == path[i-1]){
                    cout<<Neighbor[path[i]][j].second<<endl;
                    break;
                }
            }
        }
    }
    cout<<endl;
    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;
}

//function of retrieving the shortest path
/*void Graph::RetrievePaths(int ID1, int ID2, vector<unordered_set<int>> & prece){
    //path retrieval
    vector<vector<int>> paths;
    vector<int> path;
    vector<bool> ifBranch(node_num,false);
    stack<pair<int,vector<int>>> pathBranch;
    int preID;
    stack<int> st;
    st.push(ID2);
    int tempID;
    while(!st.empty()){
        preID=st.top();
        st.pop();

        if(preID!=ID1){
            path.push_back(preID);
//            if(ifBranch[preID]){
//                pathBranch.push(make_pair(preID,path));
//            }
        }else{
            path.push_back(preID);
            paths.push_back(path);
            path.clear();
            if(!pathBranch.empty() && !st.empty() && prece[pathBranch.top().first].find(st.top())==prece[pathBranch.top().first].end()){
                pathBranch.pop();
            }
            if(pathBranch.empty()){
                path.push_back(ID2);
            }else{
                path=pathBranch.top().second;
            }
            continue;
        }
        if(prece[preID].size()>1){
            ifBranch[preID]=true;
            pathBranch.push(make_pair(preID,path));
            cout<<"Branch vertex: "<<preID<<" "<<prece[preID].size()<<" ; ";
            for(auto it=prece[preID].begin();it!=prece[preID].end();++it){
                cout<<*it<<" ";
            }
            cout<<endl;
        }
        for(auto it=prece[preID].begin();it!=prece[preID].end();++it){
            st.push(*it);
        }
    }

    pair<int,int> highestVertex(-1,0);//ID, order
    cout<<"There are "<<paths.size()<<" shortest paths from "<<ID1<<" to "<<ID2<<": "<<endl;
    for(int j=0;j<paths.size();++j){
        path=move(paths[j]);
        cout<<"Path "<<j<<":\t";
        for(int i=path.size()-1;i>-1;i--){
            cout<<path[i];//<<endl;
//            if(NodeOrder[path[i]] > highestVertex.second){
//                highestVertex.second = NodeOrder[path[i]];
//                highestVertex.first = path[i];
//            }
            if(i>0){
                bool flag=false;
                for(int j=0;j<Neighbor[path[i]].size();++j){
                    if(Neighbor[path[i]][j].first == path[i-1]){
                        cout<<"("<<Neighbor[path[i]][j].second<<")\t";
                        flag=true;
                        break;
                    }
                }
                if(!flag){
                    cout<<"Not find neighbor "<<path[i-1]<<" for "<<path[i]<<endl;
                    cout<<prece[path[i]].size()<<":\t";
                    for(auto it=prece[path[i]].begin();it!=prece[path[i]].end();++it){
                        cout<<*it<<"\t";
                    }
                    cout<<"; ";
                    for(auto it=Neighbor[path[i]].begin();it!=Neighbor[path[i]].end();++it){
                        cout<<it->first<<"("<<it->second<<")\t";
                    }
                    cout<<endl;
                    cout<<prece[path[i-1]].size()<<":\t";
                    for(auto k=prece[path[i-1]].begin();k!=prece[path[i-1]].end();++k){
                        cout<<*k<<"\t";
                    }
                    cout<<"; ";
                    for(auto it=Neighbor[path[i-1]].begin();it!=Neighbor[path[i-1]].end();++it){
                        cout<<it->first<<"("<<it->second<<")\t";
                    }
                    cout<<endl;
                    exit(1);
                }
            }else{
                cout<<endl;
            }
        }
    }
//    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;
}*/
//function of retrieving the shortest path, correct
void Graph::RetrievePaths(int ID1, int ID2, vector<unordered_set<int>> & prece){
    //path retrieval
    vector<vector<int>> paths;
    vector<int> path;
    vector<bool> ifBranch(node_num,false);
    stack<pair<int,vector<int>>> pathBranch;

    int preID;

    stack<int> st;
    st.push(ID2);
    while(!st.empty()){
        preID=st.top();
        st.pop();

        if(preID!=ID1){
            path.push_back(preID);
        }else{//if preID==ID1
            path.push_back(preID);
            paths.push_back(path);
            path.clear();

            if(pathBranch.empty()){
                path.push_back(ID2);
            }else{
                while(!pathBranch.empty() && !st.empty() && pathBranch.top().first!=st.top()){
                    pathBranch.pop();
                    path=pathBranch.top().second;
                }
            }
            continue;
        }
        if(prece[preID].size()>1){
            ifBranch[preID]=true;
        }
        for(auto it=prece[preID].begin();it!=prece[preID].end();++it){
            st.push(*it);
            if(ifBranch[preID]){
                pathBranch.push(make_pair(*it,path));
            }
        }
    }

    pair<int,int> highestVertex(-1,0);//ID, order
    cout<<"There are "<<paths.size()<<" shortest paths from "<<ID1<<" to "<<ID2<<": "<<endl;
    for(int j=0;j<paths.size();++j){
        path=move(paths[j]);
        cout<<"Path "<<j<<":\t";
        for(int i=path.size()-1;i>-1;i--){
            cout<<path[i];//<<endl;
//            if(NodeOrder[path[i]] > highestVertex.second){
//                highestVertex.second = NodeOrder[path[i]];
//                highestVertex.first = path[i];
//            }
            if(i>0){
                bool flag=false;
                for(int j=0;j<Neighbor[path[i]].size();++j){
                    if(Neighbor[path[i]][j].first == path[i-1]){
                        cout<<"("<<Neighbor[path[i]][j].second<<")\t";
                        flag=true;
                        break;
                    }
                }
                if(!flag){
                    cout<<"Not find neighbor "<<path[i-1]<<" for "<<path[i]<<endl;
                    cout<<prece[path[i]].size()<<":\t";
                    for(auto it=prece[path[i]].begin();it!=prece[path[i]].end();++it){
                        cout<<*it<<"\t";
                    }
                    cout<<"; ";
                    for(auto it=Neighbor[path[i]].begin();it!=Neighbor[path[i]].end();++it){
                        cout<<it->first<<"("<<it->second<<")\t";
                    }
                    cout<<endl;
                    cout<<prece[path[i-1]].size()<<":\t";
                    for(auto k=prece[path[i-1]].begin();k!=prece[path[i-1]].end();++k){
                        cout<<*k<<"\t";
                    }
                    cout<<"; ";
                    for(auto it=Neighbor[path[i-1]].begin();it!=Neighbor[path[i-1]].end();++it){
                        cout<<it->first<<"("<<it->second<<")\t";
                    }
                    cout<<endl;
                    exit(1);
                }
            }else{
                cout<<endl;
            }
        }
    }
//    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;
}

#endif //INDEXFREE_HPP