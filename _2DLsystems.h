//
// Created by intev on 21/02/2021.
//

#ifndef ENGINE__2DLSYSTEMS_H
#define ENGINE__2DLSYSTEMS_H
#include <string>
#include <iostream>
#include "l_parser.h"
#include <fstream>
#include <set>
#include <vector>
#include "draw2DLine.h"
using namespace std;

class TweeDLSystem{
public:

    //hoek die die richting bepaalt
    double Angle;

    //bepaalt de start richting
    double StartingAngle;

    //een string die bestaat uit symbolen van het alfabet
    string Initiator;

    //aantal herhalening
    unsigned int Iterations;

    //een lijst van symbolen
    set<char> Alfabet;

    //stack voor punt en hoek op te slaan
    vector<pair<Point2D, double>> stack;

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

    void parse2DL(const string& L2DinputFile);

    Lines2D createDrawVector(Color lineColor);

    string createDrawString(const string& Initia);

    int topStack = -1;

};





#endif //ENGINE__2DLSYSTEMS_H
