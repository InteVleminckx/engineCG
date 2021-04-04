//
// Created by intev on 19/02/2021.
//

#ifndef ENGINE_HULPFILE_H
#define ENGINE_HULPFILE_H
#include "easy_image.h"
#include <cmath>
#include "draw2DLine.h"
#include <vector>
#include <map>
//#include "drieDdrawings.h"
using namespace std;


class Line2D;

int roundToInt(double d);

double gradesToRad(double graden);

bool isOperator(char &iniChara);

bool isOperator3DL(char &iniChara);

pair<pair<double, double>,pair<double, double>> getMinMax(vector<Line2D> &lijnen);

vector<double> moveImageToFit(vector<Line2D> &lijnen, const int &size);

vector<double> moveImageToFitZb(const int &size, const double &Xmin,const double &Xmax, const double &Ymin, const double &Ymax);


#endif //ENGINE_HULPFILE_H
