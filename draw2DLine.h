//
// Created by intev on 19/02/2021.
//

#ifndef ENGINE_DRAW2DLINE_H
#define ENGINE_DRAW2DLINE_H
#include "easy_image.h"
#include <list>
#include "hulpFile.h"
#include <limits>
#include "assert.h"

using namespace std;

class Color{
public:
    double red;
    double green;
    double blue;
};

class Point2D{
public:
    double x;
    double y;


};

class Line2D{
public:

    Point2D p1, p2;
    double z1, z2;

    Color color;

    Line2D(Point2D punt1, Point2D punt2, Color lineColor, double zp1 = 0, double zp2 = 0);



};

class ZBuffer{
public:


    ZBuffer() {}

    vector<vector<double>> zBuffer;


    ZBuffer(const int width, const int height);

};


using Lines2D = list<Line2D>;

img::EasyImage draw2DLines(const Lines2D &lines, const int size, const img::Color& background, bool zBuf);

void draw_zbuf_line(ZBuffer &, img::EasyImage &,  unsigned int x0,  unsigned int y0,
                    double z0,  unsigned int x1,  unsigned int y1,  double z1,
                    const img::Color &color);



#endif //ENGINE_DRAW2DLINE_H
