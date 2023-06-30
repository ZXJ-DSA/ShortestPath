#include <iostream>

#include "Preprocess.hpp"

int main(int argc, char** argv){

    if( argc <4 || argc>5){
        printf("usage:\n<arg1> Source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> Dataset, e.g. NY\n");
        printf("<arg3> Task, 1: Graph Data Generation; 2: Query Pairs Generation; 3: Update Pairs Generation; 4: Vertex Ordering Generation. e.g. 1\n");
        printf("<arg4> Vertex ordering strategy (optional), 1: MDE-based ordering; 2: Degree-based ordering. e.g. 2\n");
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
            cout << "argv[3]: " << argv[3] << endl;//task
            g.algo = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//strategy
            g.threadnum = stoi(argv[4]);
        }
    }


    Timer tt;
    tt.start();
    cout<<"Dataset: "<<g.dataset<<endl;
    g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;

    if(g.algo == 1){
        g.GraphDataGen(source_path + "/" + g.dataset + "/"+ g.dataset, 1);//1
    }else{
        g.Preprocessing();
    }


    tt.stop();
    cout<<"\nOverall time: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"-----------------"<<endl;

    return 0;
}
