#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>

using namespace std;
#define RAND_MAX = 100000
#define DOMAIN 10000
#define LENGTH 501

int main(int argc, char **argv)
{
    int opt;
    const char *optstring = "t:n:o:";
    string tplt;
    string oname;
    int num;
    while ((opt = getopt(argc, argv, optstring)) != -1)
    {
        switch (opt)
        {
        case 't':
            tplt = optarg;
            break;
        case 'n':
            num = atoi(optarg);
            break;
        case 'o':
            oname = optarg;
            break;
        }
    }
    // srand(time(0));
    srand(time(NULL));
    ifstream fin(tplt.c_str());
    ofstream fout(oname.c_str(), ios::out);
    string temp;
    int flag = 1;
    for (int i = 0; i < num; i++)
    {
        getline(fin, temp);
        fout << temp;
        fout << endl;
        for (int j = 0; j < LENGTH; j++)
        {
            int randomNum = rand() % DOMAIN;
            if (randomNum < 2866)
            {
                fout << 'A';
            }
            else if (randomNum < 4998)
            {
                fout << 'C';
            }
            else if (randomNum < 7134)
            {
                fout << 'G';
            }
            else
            {
                fout << 'T';
            }
        }
        fout << endl;

        getline(fin, temp);
    }
    fin.close();
    fout.close();
    cout << "Time used: " << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
}