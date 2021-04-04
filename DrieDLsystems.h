//
// Created by inte on 14.03.21.
//

#ifndef ENGINE_DRIEDLSYSTEMS_H
#define ENGINE_DRIEDLSYSTEMS_H
#include "vector3d.h"
#include <string>
#include <iostream>
#include "l_parser.h"
#include <fstream>
#include <set>
#include <vector>
#include "ini_configuration.h"
#include "easy_image.h"
#include "draw2DLine.h"

using namespace std;

class Face;
class Figure;

class DrieDLsystems {
    Matrix rotatieMatrix;
public:

    //hoek die die richting bepaalt
    double Angle;

    //een string die bestaat uit symbolen van het alfabet
    string Initiator;

    //aantal herhalening
    unsigned int Iterations;

    //een lijst van symbolen
    set<char> Alfabet;

    //stack voor punt en hoek op te slaan
    vector<vector<Vector3D>> stack;

    //voor elk symbool is er 1 replacement rule die het symbool afbeeldt op een string
    //die bestaat uit symbolen die behoren tot het alfabet
//    vector<pair<char, string>> Replacements;

    //voor elk symbool is er 1 replacement rule die het symbool afbeeldt op een string
    //die bestaat uit symbolen die behoren tot het alfabet
    map<char, vector<pair<string, double>>> Replacements;


    //beeldt elk symbool uit het Alfabet uit met 0 of 1
//    vector<pair<char, int>> Draw;

    //beeldt elk symbool uit het Alfabet uit met 0 of 1
    map<char, int> Draw;
    DrieDLsystems(Matrix &rotMatrix);
    void parse3DL(const string& L3DinputFile);

    Figure create3DLsystem();

    string createDrawString(string &initiator);

    pair<vector<Vector3D>,vector<Face>> create3DPointsVector();

    int topStack = -1;

};


#endif //ENGINE_DRIEDLSYSTEMS_H
