#include <fstream>
#include "tools.h"

using namespace std;

int loadFramesFromFile(vector<Frame>& frames, string path) {
    ifstream inFile;
    inFile.open(path);

    if (!inFile) {
        printf("error opening file\n");
        return 0;
    }
    double x[2];
    double y[2];
    while (inFile >> x[0] >> y[0] >> x[1] >> y[1] ) {
        Frame f;
        f.x[0] = x[0];
        f.x[1] = x[1];
        f.y[0] = y[0];
        f.y[1] = y[1];
        f.marks = 0;
        frames.push_back(f);
    }
    return 0;
};

int loadSupportingAreasFromFile(vector<SupportingArea>& vertical, vector<SupportingArea>& horizontal, string path) {
    ifstream inFile;
    inFile.open(path);

    if (!inFile) {
        printf("error opening file\n");
        return 0;
    }
    double x[2];
    double y[2];
    char c;
    while (inFile >> x[0] >> y[0] >> x[1] >> y[1] >> c) {
        SupportingArea sa;
        sa.x[0] = x[0];
        sa.x[1] = x[1];
        sa.y[0] = y[0];
        sa.y[1] = y[1];
        sa.load = 0;
        switch(c) {
            case 'v':
                sa.dir = orientation::vertical;
                vertical.push_back(sa);
                break;
            case 'h':
                sa.dir = orientation::horizontal;
                horizontal.push_back(sa);
                break;
            default:
                // do nothing
                break;
        }
    }
    return 0;
}

int writeArrayToFile(vector<vector<double>> array, int m, int n, string path) {
    ofstream outFile(path);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            outFile << array[i][j] << ' ';
        }
        outFile << endl;
    }
    outFile.close();
    return 0;
}