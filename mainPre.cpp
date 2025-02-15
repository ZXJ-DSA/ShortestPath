#include <iostream>

#include "Preprocess.hpp"

int main(int argc, char** argv){

    if( argc <4 || argc > 8){
        printf("usage:\n<arg1> Source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> Dataset, e.g. NY\n");
        printf("<arg3> Task, 1: Graph Data Generation; 2: Query Pairs Generation; 3: Update Pairs Generation; 4: Vertex Ordering Generation; 5: Road Network Preprocessing; 6: Compute the diameter; 7: Extract a smaller subgraph for scalability test; 8: Generate queries with varying distance; 9: Maximal step number of all shortest paths e.g. 1\n");
        printf("<arg4> Vertex ordering strategy (optional), 1: MDE-based ordering; 2: Degree-based ordering. e.g. 2\n");
        printf("<arg5> Edge generation strategy (optional), 1: random; 2: degree-proportional; 3: degree as numerator; 4: Degree-based ranged random. e.g. 2\n");
        printf("<arg6> Vertex number for graph data generation (optional)\n");
        printf("<arg7> Percentage of scalability test (optional), eg: 20\n");
        exit(0);
    }
    cout << "Hello, World!" << endl;
    string source_path = "/Users/zhouxj/Documents/1-Research/Datasets";
    Graph g;
    g.dataset = "NY";
    g.task = 1;
    int estrategy=2;
    int percent=0;

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
            g.task = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//vertex ordering strategy
            g.strategy = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5]: " << argv[5] << endl;//edge generation strategy
            estrategy = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6]: " << argv[6] << endl;//vertex number
            g.node_num = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7]: " << argv[7] << endl;//percentage
            percent = stoi(argv[7]);
        }
    }


    Timer tt;
    tt.start();
    cout<<"Dataset: "<<g.dataset<<endl;
    if(percent==0){
        g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;
    }
    else{
        g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset+"_"+ to_string(percent);
    }
//    g.graph_path = source_path + "/" + g.dataset;
    cout<<g.graph_path<<endl;
    if(g.task == 1){
        g.GraphDataGen(source_path + "/" + g.dataset + "/"+ g.dataset, estrategy);//1
    }else{
        g.Preprocessing();
    }


    tt.stop();
    cout<<"\nOverall time: "<<tt.GetRuntime()<<" s."<<endl;
    cout<<"-----------------"<<endl;

    return 0;
}
