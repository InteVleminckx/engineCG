//
// Created by intev on 19/02/2021.
//
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998
#include "hulpFile.h"

int roundToInt(double d) {return static_cast<int>(round(d));}

double gradesToRad(double graden){return 2*M_PI*(graden/360);}

bool isOperator(char &iniChara){
    if (iniChara == '+' || iniChara == '-' || iniChara == '(' || iniChara == ')' ){return true;}
    return false;
}

pair<pair<double, double>,pair<double, double>> getMinMax(vector<Line2D> &lijnen){

    double xMax, xMin, yMax, yMin;

    for (auto & allLine : lijnen) {

        if  (allLine.p1.x > xMax){xMax = allLine.p1.x;}

        if  (allLine.p1.x < xMin){xMin = allLine.p1.x;}

        if  (allLine.p2.x > xMax){xMax = allLine.p2.x;}

        if  (allLine.p2.x < xMin){xMin = allLine.p2.x;}

        if  (allLine.p1.y > yMax){yMax = allLine.p1.y;}

        if  (allLine.p1.y < yMin){yMin = allLine.p1.y;}

        if  (allLine.p2.y > yMax){yMax = allLine.p2.y;}

        if  (allLine.p2.y < yMin){yMin = allLine.p2.y;}

    }

    return make_pair(make_pair(xMin, xMax), make_pair(yMin, yMax));

}

vector<double> moveImageToFit(vector<Line2D> &lijnen, const int &size){

    vector<double> props;

    //bepaalt minX, maxX, minY en maxY
    pair<pair<double, double>,pair<double, double>> MinMax = getMinMax(lijnen);

    //range
    double xRange = MinMax.first.second - MinMax.first.first;
    double yRange = MinMax.second.second - MinMax.second.first;

    //breedte en hoogte image
    double imageX = size*(xRange/max(xRange,yRange));
    double imageY = size*(yRange/max(xRange,yRange));

    //scaleFactor
    double scaleF = 0.95*(imageX/xRange);

    //centerpunt
    double DCx = scaleF*((MinMax.first.second+MinMax.first.first)/2);
    double DCy = scaleF*((MinMax.second.second+MinMax.second.first)/2);

    //verschuivings afstand
    double dx = (imageX/2)-DCx; double dy = (imageY/2)-DCy;


    props.push_back(imageX);
    props.push_back(imageY);
    props.push_back(scaleF);
    props.push_back(dx);
    props.push_back(dy);

    return props;

}

bool isOperator3DL(char &iniChara){

    if (iniChara == '+' || iniChara == '-' || iniChara == '(' || iniChara == ')' ||
        iniChara == '^' || iniChara == '&' || iniChara == 92 || iniChara == '/' || iniChara == '|'){return true;}
    return false;


}



vector<double> moveImageToFitZb(const int &size, const double &Xmin,const double &Xmax, const double &Ymin, const double &Ymax){

    vector<double> props;

    //range
    double xRange = Xmax - Xmin;
    double yRange = Ymax - Ymin;

    //breedte en hoogte image
    double imageX = size*(xRange/max(xRange,yRange));
    double imageY = size*(yRange/max(xRange,yRange));

    //scaleFactor
    double scaleF = 0.95*(imageX/xRange);

    //centerpunt
    double DCx = scaleF*((Xmax+Xmin)/2);
    double DCy = scaleF*((Ymax+Ymin)/2);

    //verschuivings afstand
    double dx = (imageX/2)-DCx; double dy = (imageY/2)-DCy;


    props.push_back(imageX);
    props.push_back(imageY);
    props.push_back(scaleF);
    props.push_back(dx);
    props.push_back(dy);

    return props;


}