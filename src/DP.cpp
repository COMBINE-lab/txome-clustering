#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/range/irange.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

namespace rapmap{
size_t dimension;
//taken from http://cpptruths.blogspot.com/2011/10/multi-dimensional-arrays-in-c11.html
//template <class T, uint32_t ROW, size_t COL>
//using Matrix = std::array<std::array<T, COL>, ROW>;

matrix<uint32_t> getDistacneMatrix(std::vector<std::vector<uint32_t>> eqClass){
    //calculates the distance matrix out of eqClass

    matrix<uint32_t> disMatrix(dimension, dimension);
    //assign zeros to all element
    disMatrix = zero_matrix<uint32_t>(dimension, dimension) ;
    //calculates distance Mtrix
    for (int32_t i : boost::irange(0,  static_cast<int32_t>(dimension))){
        for (int32_t j : boost::irange(i+1, static_cast<int32_t>(dimension))){
            if (eqClass[i] == eqClass[j])
                disMatrix(i, j) = disMatrix(i, j-1) + 1;
            else
                disMatrix(i, j) = disMatrix(i, j-1);
        }
    }

    for (int32_t i : boost::irange(0,  static_cast<int32_t>(dimension))){
        for (int32_t j : boost::irange(i+1, static_cast<int32_t>(dimension))){
            if (disMatrix(i, j) != disMatrix(i, j-1))
                disMatrix(i, j) += disMatrix(i, j-1);
            else
                disMatrix(i, j) = disMatrix(i, j-1);
        }
    }

    std::cerr<<"Done Creating Distance Matrix"<<std::endl;
    std::cout<<disMatrix;
    return disMatrix;
}

void get_opt(std::vector<std::vector<uint32_t>> eqClass){
    //not using vector since std::array is faster and uses stack

    matrix<uint32_t> disMatrix(dimension, dimension);
    matrix<uint32_t> scoreMatrix(dimension+1, dimension);
    matrix<uint32_t> opt(dimension+1, dimension);

    //assigning starting cases and zero matrixes
    disMatrix = getDistacneMatrix(eqClass);
    scoreMatrix = zero_matrix<uint32_t>(dimension+1, dimension) ;
    opt = zero_matrix<uint32_t>(dimension+1, dimension) ;

    //assign score for 1 cluster from 0th row of distance matrix
    matrix_row<matrix<uint32_t>> scoreRow(scoreMatrix, 1);
    matrix_row<matrix<uint32_t>> disRow(disMatrix, 0);
    scoreRow = disRow;

    //assign opt direction for l==1
    matrix_row<matrix<uint32_t>> optRow(opt, 1);
    auto iterator = scoreRow.begin();
    for (size_t i = 0; i < optRow.size(); ++i){
        auto max =  std::max_element(iterator, iterator+i+1);
        opt(1, i) = std::distance(iterator, max);
    }

    uint32_t maxScore, maxIndex, score, clusters;
    //start filling matrix from cluster l=[2, dimension)
    for (size_t l : boost::irange(2, static_cast<int32_t>(dimension+1))){
        for (int32_t i : boost::irange(0, static_cast<int32_t>(dimension))){
            maxScore = 0;
            maxIndex = 0;
            for (size_t j : boost::irange(1, i+1)){
                if ((l-1) > (j-1))
                    continue;
                //get score by breaking at location j that will give
                //l-1 cluster from [0, j-1] and one cluster [j, i]
                score = scoreMatrix(l-1, j-1) + disMatrix(j, i);
                if (score > maxScore){
                    maxScore = score;
                    maxIndex = j-1;
                }
                if (score == 0 and maxScore == 0)
                    maxIndex = j-1;
            }
            scoreMatrix(l, i) = maxScore;
            opt(l, i) = maxIndex;
        }
        std::cout<<std::endl<<scoreMatrix(l, dimension-1)<<std::endl;
        std::cerr<<"\r\rDone for "<<l<<"/"<<dimension-1<<std::flush;
        if (scoreMatrix(l, dimension-1) < scoreMatrix(l-1, dimension-1)){
            //clusters = l-1;
            break;
        }
    }
    std::cerr<<std::endl;
    //std::cout<<scoreMatrix<<opt<<std::endl;

    //get the clusters by parsing the matrixes
    uint32_t N, i;
    N = dimension-1;

    //Below commented code if we are running from 2-> l
    matrix_column<matrix<uint32_t>> scoreCol(scoreMatrix, N);
    auto max = std::max_element(scoreCol.begin(), scoreCol.end());
    clusters = std::distance(scoreCol.begin(), max);
    i = opt(clusters, N);


    //Do the same thing until we get interval for all optimal clusters
    while(clusters>1){
        clusters -= 1;
        //std::cout<<"("<<i+1<<","<<N<<")"<<std::endl;
        std::cout<<i+1<<"\t"<<N<<std::endl;
        N = i;
        //get the column N for scoreMatrix
        matrix_column<matrix<uint32_t>> scoreCol(scoreMatrix, N);
        //get the biggest value's index at save it in i
        auto max = std::max_element(scoreCol.begin(), scoreCol.end());
        i = opt(std::distance(scoreCol.begin(), max), N);
    }
    std::cout<<0<<"\t"<<N<<std::endl;
    //std::cout<<"("<<0<<","<<N<<")"<<std::endl;
}

std::vector<std::vector<uint32_t>>
    parseFasta(std::string fileName){

    using namespace boost::algorithm;

    //eqClass contains for each base which contigs does not have -
    std::vector<std::vector<uint32_t>> eqClass;
    uint32_t cIndex = 0;
    uint32_t bIndex = 0;
    std::string sequence;

    std::ifstream fileHandle (fileName);

    if (fileHandle.is_open()){
        cIndex = 0;
        getline(fileHandle, sequence);
        while(getline(fileHandle, sequence)){
            trim_right(sequence);
            while (sequence[0] != '>'){
                //read different line of sequqnce
                bIndex = 0;
                for (auto base:sequence){
                    //for each base in the sequnce line
                    if (cIndex == 0){
                        //inside this for 1st contig only
                        std::vector<uint32_t> tempVector;
                        //check for -
                        if (base != '-'){
                            tempVector.push_back(0);
                            eqClass.push_back(tempVector);
                        }
                        else{
                            eqClass.push_back(tempVector);
                        }
                    }
                    else
                        //all other contigs except 1 goes here
                         if (base != '-'){
                             eqClass[bIndex].push_back(cIndex);
                         }
                    bIndex += 1;
                }
                getline(fileHandle, sequence);
                trim_right(sequence);
                if (sequence.empty())
                    break;
            }
            cIndex += 1;
        }
        std::cerr<<"Done Creating classes\n";
    }
    else {
        std::cerr<<"error reading fasta file\n";
        exit(0);
    }
    //for (auto x: eqClass){
    //    for (auto y : x){
    //        std::cout<<y<<",";
    //    }
    //    std::cout<<"\n";
    //}
    return eqClass;
}
}

int main( int argc, char* argv[]  ) {

    using std::string;

    string fileName = argv[1];
    std::vector<std::vector<uint32_t>> eqClass;
    eqClass = rapmap::parseFasta(fileName);

    rapmap::dimension = eqClass.size();
    rapmap::get_opt(eqClass);
    return 0;
}
