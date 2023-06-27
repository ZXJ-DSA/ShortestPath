#include <iostream>

#include "Founda.hpp"
#include "IndexFree.hpp"
#include "IndexBased.hpp"

int main(int argc, char** argv){

    if( argc != 4){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> algorithm, 0: Graph Data Generation; 1: Dijkstra; 2: Bi-Dijkstra; 3: A*; 4: Bi-A*; 5: Bellman-Ford; 6: CH; 7: H2H; 8: PLL. e.g. 1\n");
//        printf("<arg4> thread number(optional), e.g. 10\n");
        exit(0);
    }
    cout << "Hello, World!" << endl;
    string source_path = "/Users/zhouxj/Documents/1-Research/Datasets";
    Graph g;
    g.dataset = "NY";
    g.algo = 1;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        source_path = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            g.dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//algorithm
            g.algo = stoi(argv[3]);
        }
    }


    Timer tt;
    tt.start();
    cout<<"Dataset: "<<g.dataset<<endl;
    g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;

    if(g.algo == 0){
        g.GraphDataGen(source_path + "/" + g.dataset + "/"+ g.dataset);
        return 0;
    }

    g.ReadGraph(g.graph_path);
    g.Preprocess();

    g.CorrectnessCheck(100);
//    g.EfficiencyTest(1000);

    tt.stop();
    cout<<"\nOverall time: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"-----------------"<<endl;

    return 0;
}
