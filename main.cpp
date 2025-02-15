#include <iostream>

#include "Founda.hpp"
#include "IndexFree.hpp"
#include "H2H.hpp"
#include "PLL.hpp"


int main(int argc, char** argv){

    if( argc <4 || argc>11){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> algorithm, 1: Dijkstra; 2: Bi-Dijkstra; 3: A*; 4: Bi-A*; 5: Bellman-Ford; 6: CH; 7: CHWP; 8:H2H; 9: PLL; 10: PSL; 11: PCL; 12: GLL. e.g. 1\n");
        printf("<arg4> vertex ordering strategy (optional), 0: read from disk; 1: MDE-based; 2:Degree-based. e.g. 0\n");
        printf("<arg5> if dynamic networks (optional), 0: No; 1: Yes. Default: 0\n");
        printf("<arg6> update type (optional), 0: Decrease; 1: Increase. Default: 0\n");
        printf("<arg7> update batch number (optional), e.g. 10\n");
        printf("<arg8> update batch size (optional), e.g. 10\n");
        printf("<arg9> thread number (optional), e.g. 10\n");
        printf("<arg10> batch size of PCL (optional), e.g. 10\n");

        exit(0);
    }
    cout << "Hello, World!" << endl;
    string source_path = "/Users/zhouxj/Documents/1-Research/Datasets";
    string dataset;
    int algo=1;
    int threadnum=15;
    int batchsize=15;
    int batchNum=15;
    int batchSize=10;
    int orderStrategy=2;
    bool ifDynamic=false;
    int updateType=0;//decrease
    updateType=1;//increase;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (source path): " << argv[1] << endl;//source path
        source_path = argv[1];
        if(argc > 2){
            cout << "argv[2] (dataset): " << argv[2] << endl;//dataset
            dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3] (algorithm): " << argv[3] << endl;//algorithm
            algo = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4] (vertex ordering): " << argv[4] << endl;//vertex ordering strategy
            orderStrategy = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5] (if dynamic): " << argv[5] << endl;//if dynamic
            ifDynamic = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (update type): " << argv[6] << endl;//update type
            updateType = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (update batch number): " << argv[7] << endl;//thread number
            batchNum = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8] (update batch size): " << argv[8] << endl;//batch size
            batchSize = stoi(argv[8]);
        }
        if(argc > 9){
            cout << "argv[9] (thread number): " << argv[9] << endl;//thread number
            threadnum = stoi(argv[9]);
        }
        if(argc > 10){
            cout << "argv[10] (batch size): " << argv[10] << endl;//batch size
            batchsize = stoi(argv[10]);
        }

    }


    Timer tt;
    tt.start();
    Graph g;
    g.dataset = dataset;
    g.algo = algo;
    g.threadnum = threadnum;
    if(argc > 5){
        g.batchsize=batchsize;
    }else{
        g.batchsize=g.threadnum;
    }
    g.orderStrategy=orderStrategy;
    cout<<"Dataset: "<<g.dataset<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    cout<<"Batch size: "<<g.batchsize<<endl;
    g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;
    g.ifDynamic = ifDynamic;//for dynamic graph

//    vector<pair<int,unsigned long long int>> weights;
//    for(int i=0;i<10;++i){
//        int wi=rand()%10000;
//        weights.emplace_back(wi, i*10);
//    }
//    sort(weights.begin(),weights.end());
//    for(int i=0;i<weights.size();++i){
//        cout<<weights[i].first<<"("<<weights[i].second<<") ";
//    }
//    cout<<endl;

    g.ReadGraph(g.graph_path);
    g.Preprocess();

//    g.CorrectnessCheck(100);
    g.EfficiencyTest(100);//1000
//    exit(0);

    bool ifBatch=true;

//    g.H2HIndexCompare(g.graph_path+".labelTemp3");
//    g.CHIndexCompare(g.graph_path+".CHIndex4");
//    exit(0);
    if(dataset=="beijing" || dataset=="Guangdong"){//real-life updates
        g.IndexMaintenanceRealUpdate(source_path+"/"+dataset+"/"+dataset+"_20160105.updates");
    }else {
        g.IndexMaintenance(updateType, ifBatch, batchNum, batchSize);//index maintenance
    }

    tt.stop();
    cout<<"\nOverall time: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"-----------------\n"<<endl;

    return 0;
}
