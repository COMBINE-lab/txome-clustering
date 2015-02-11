// This file is used to read the fastq files and returns

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main (int argc, const char* argv[]) {
    fstream readFile;
    string read;
    readFile.open ( "input.fastq" );
    if ( readFile.is_open() ) {
        getline ( readFile, read );
        while ( getline ( readFile, read ) ) {
            cout << read << '\n';
            getline ( readFile, read );
            getline ( readFile, read );
            getline ( readFile, read );
        }
        readFile.close();
    }
    else    {
        cout << "unable to open file";
    }
    return 0;
}
