#ifndef TOOLS_H
#define TOOLS_H

#include "facade.h"

int loadFramesFromFile(vector<Frame>& frames, string path);
int loadSupportingAreasFromFile(vector<SupportingArea>& vertical, vector<SupportingArea>& horizontal, string path);
int writeArrayToFile(vector<vector<double>>, int, int, string);
#endif // TOOLS_H