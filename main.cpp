#include <iostream>

#include "Founda.hpp"
#include "IndexFree.hpp"
#include "IndexBased.hpp"

int main(int argc, char** argv){

    if( argc <4 || argc>6){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> algorithm, 1: Dijkstra; 2: Bi-Dijkstra; 3: A*; 4: Bi-A*; 5: Bellman-Ford; 6: CH; 7: PLL; 8: PSL. e.g. 1\n");
        printf("<arg4> thread number (optional), e.g. 10\n");
        printf("<arg6> vertex ordering strategy (optional), 0: read from disk; 1: MDE-based; 2:Degree-based. e.g. 0\n");
        exit(0);
    }
    cout << "Hello, World!" << endl;
    string source_path = "/Users/zhouxj/Documents/1-Research/Datasets";
    Graph g;
    g.dataset = "NY";
    g.algo = 1;
    g.threadnum = 15;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        source_path = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            g.dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3] (Algorithm): " << argv[3] << endl;//algorithm
            g.algo = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//algorithm
            g.threadnum = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//vertex ordering strategy
            g.orderStrategy = stoi(argv[5]);
        }
    }


    Timer tt;
    tt.start();
    cout<<"Dataset: "<<g.dataset<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;

    g.ReadGraph(g.graph_path);
    g.Preprocess();

    g.CorrectnessCheck(100);
//    g.EfficiencyTest(1000);

    tt.stop();
    cout<<"\nOverall time: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"-----------------"<<endl;

    return 0;
}
