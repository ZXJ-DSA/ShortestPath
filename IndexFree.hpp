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
    vector<pair<int,int>>::iterator it;
    int NNodeID,NWeigh;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!pqueue.empty()){
        pqueue.extract_min(topNodeID, topNodeDis);
        if(topNodeID==ID2){
            d=distance[ID2];
            break;
        }
        closed[topNodeID]=true;

        for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
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
                    cout<<topNodeID<<" "<<NNodeID<<" "<<distance[topNodeID]<<" "<<distance[NNodeID]<<endl;

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


#endif //INDEXFREE_HPP