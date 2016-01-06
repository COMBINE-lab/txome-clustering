#include <GATB_1.hpp>

#include <fstream>
#include <queue>
#include <stack>
#include <map>
#include <vector>
using namespace std;
#define DEBUG(a)  //a
#define INFO(a)   //a
/********************************************************************************/
const char* STR_NODE_TYPE = "-type";
/********************************************************************************/
const char * tr = "";

class DotGeneratorTool : public Tool
{
public:

    // Constructor
    DotGeneratorTool () : Tool ("DotGenerator")
    {
        _parser->push_front (new OptionOneParam (STR_URI_GRAPH,  "graph file", true ));
        _parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "dot file",  false ));
        _parser->push_front (new OptionOneParam (STR_NODE_TYPE,  "node type (0: all,  1:branching)", false, "0" ));
    }



//    bool & operator==(const Graph::Vector<Node>& lhs, const Graph::Vector<Node>& rhs){
//        bool flag = 1;
//        for(Graph::Vector<Node>::iterator it =lhs.begin(); it != lhs.end(); ++it) {
//            if(std::find(rhs.begin(), rhs.end(), *it) == rhs.end()){
//                flag = 0;
//                break;
//            }
//        }
//        if (flag)
//            return true;
//        else
//            return false;
//    }




    template<typename NodeType>
    void process (const char* name)
    {
        string outputFile = getInput()->get(STR_URI_OUTPUT) ?
            getInput()->getStr(STR_URI_OUTPUT) :
            (System::file().getBaseName(getInput()->getStr(STR_URI_GRAPH)) + ".dot");
        ofstream output (outputFile.c_str());
        output << "digraph " << name << "{\n";
        // We load the graph
        Graph graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
        map<Node, u_int64_t> mapping;
        u_int64_t count = 1;
        Graph::Iterator<NodeType> itMap = graph.iterator<NodeType> ();

        for (itMap.first(); !itMap.isDone(); itMap.next())  {
            NodeType current = itMap.item();
            size_t indegree = graph.indegree(current);
            size_t outdegree = graph.outdegree(current);
//            u_int64_t matchCount = 0;
//            Graph::Vector<NodeType> predecessors = graph.predecessors<NodeType>(current);
//            Graph::Vector<NodeType> successors = graph.successors<NodeType>(current);
//
//            for(size_t i = 0; i < predecessors.size(); ++i) {
//                for(size_t j = 0; j < successors.size(); ++j){
//                    if (predecessors[i] == successors[j]){
//                        matchCount += 1;
//                        break;
//                    }
//                }
//            }
//
//            std:cout<<matchCount<<"\t"<<predecessors.size()<<endl;
//
//            if (matchCount == predecessors.size()-1) {
            if (indegree == outdegree == 1 or indegree == 0 or outdegree == 0) {
                Graph::Vector<NodeType> predecessors = graph.predecessors<NodeType>(current);
                Graph::Vector<NodeType> successors = graph.successors<NodeType>(current);

                auto pred = mapping.find(predecessors[0]);
                auto succ = mapping.find(successors[0]);

                if (pred != mapping.end() and current.abundance == predecessors[0].abundance ){
                    mapping[current] = mapping[predecessors[0]];
                }
                else if(succ != mapping.end() and current.abundance == successors[0].abundance){
                    mapping[current] = mapping[successors[0]];
                }
                else{
                    std::vector<NodeType> collapse;
                    collapse.push_back(current);
                    int64_t id = -1;

                    while(graph.indegree(predecessors[0]) == 1 and graph.outdegree(predecessors[0]) == 1 ){
                        collapse.push_back(predecessors[0]);
                        predecessors = graph.predecessors<NodeType>(predecessors[0]);
                        auto pred = mapping.find(predecessors[0]);
                        if (pred != mapping.end()){
                            id = mapping[predecessors[0]];
                            for( auto it = collapse.begin(); it != collapse.end(); ++it) {
                                mapping[*it] = id;
                            }
                            break;
                        }
                    }
                    while(graph.indegree(successors[0]) == 1 and graph.outdegree(successors[0]) == 1 ){
                        if (id != -1){
                            mapping[successors[0]] = id;
                            successors = graph.successors<NodeType>(successors[0]);
                            auto succ = mapping.find(successors[0]);
                            if (succ != mapping.end()){
                                break;
                            }
                        }
                        else{
                            collapse.push_back(successors[0]);
                            successors = graph.successors<NodeType>(successors[0]);
                            auto succ = mapping.find(successors[0]);
                            if (succ != mapping.end()){
                                id = mapping[successors[0]];
                                for( auto it = collapse.begin(); it != collapse.end(); ++it ) {
                                    mapping[*it] = id;
                                    break;
                                }
                            }
                        }
                    }

                    if (id == -1){
                        mapping[current] = count++;
                    }
                    //mapping[current] = count++;
                }
                //std::cout<<graph.toString(itMap.item())<<current.strand<<endl;
            }
            else{
                //std::cout<<graph.toString(current)<<indegree<<outdegree;
                //std::cout<<"avi<<<<<<<<<<<<<<<<<<<<<<\n";
                mapping[current] = count++;
            }
            //std::cout<<graph.toString(current)<<"\t"<<current.strand<<endl;
            //mapping[current] = count++;
        }

        Graph::Iterator<NodeType> it = graph.iterator<NodeType> ();
        for (it.first(); !it.isDone(); it.next())
        {
            NodeType current = it.item();
            //Graph::Vector<NodeType> neighbors = graph.successors<NodeType> (current);
            Graph::Vector<NodeType> neighbors;

            //if (graph.indegree(current) == graph.outdegree(current)){
                neighbors = graph.neighbors<NodeType> (current);
                for (size_t i=0; i<neighbors.size(); i++){
                    if (mapping[current]  != mapping[neighbors[i]]){
                    //for (size_t j=0; j<current.abundance; j++){
                        //string currKmer = graph.toString(current);
                        //string neiKmer = graph.toString(neighbors[i]);

                        //if (currKmer.substr(0, 30).compare(neiKmer.substr(1, 30)))
                        //    output << mapping[neighbors[i]] << " -> " <<  mapping[current] << " ;\n";
                        //else if (currKmer.substr(1, 30).compare(neiKmer.substr(0, 30)))
                            output << mapping[current] << " -> " <<  mapping[neighbors[i]] << " ;\n";
                        //else
                        //    std::cout<<"error";
                    //}
                    }
                }
            //}
            //else if(graph.successors<NodeType> (current).size() > graph.predecessors<NodeType> (current).size()){
            //    neighbors = graph.successors<NodeType> (current);
            //    for (size_t i=0; i<neighbors.size(); i++){
            //    //if (mapping[current]  != mapping[neighbors[i]]){
            //        for (size_t j=0; j<current.abundance; j++){
            //            output << mapping[current] << " -> " <<  mapping[neighbors[i]] << " ;\n";
            //        }
            //    //}
            //    }
            //}
            //else {
            //    neighbors = graph.predecessors<NodeType> (current);
            //    for (size_t i=0; i<neighbors.size(); i++){
            //    //if (mapping[current]  != mapping[neighbors[i]]){
            //        for (size_t j=0; j<current.abundance; j++){
            //            output << mapping[neighbors[i]] << " -> " <<  mapping[current] << " ;\n";
            //        }
            //    //}
            //    }
            //}
        }
        output << "}\n";
        output.close();

    }
    // Actual job done by the tool is here
    void execute ()
    {
        switch (getInput()->getInt(STR_NODE_TYPE))
        {
            case 0: process<Node>          ("all");        break;
            case 1: process<BranchingNode> ("branching");  break;
            default: break;
        }
     }
};
/********************************************************************************/
/*                                                                              */
/*                   Generate dot file from a graph.                            */
/*                                                                              */
/*  This snippet generates a dot file from a graph file. You can then generate  */
/*  a pdf file with "dot -Tpdf graph.dot -o graph.pdf"                          */
/*                                                                              */
/*  NOTE: de Bruijn graphs may be huge and complex, so dot is not the best tool */
/*  to display such graphs. You should use it on small graphs with only a few   */
/*  hundreds of nodes.                                                          */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        DotGeneratorTool().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;

 }
