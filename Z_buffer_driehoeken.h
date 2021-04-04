//
// Created by inte on 27.03.21.
//

#ifndef ENGINE_Z_BUFFER_DRIEHOEKEN_H
#define ENGINE_Z_BUFFER_DRIEHOEKEN_H
#include <iostream>
#include "draw2DLine.h"
#include "drieDdrawings.h"
#include "easy_image.h"

using namespace std;

class Z_driehoek {

    drieDdrawings newFaceAndPoints;
    int imageX;
    int imageY;
    ZBuffer zBuffer;
    double XMin, XMax, YMin, YMax;
    double d;
    double dx;
    double dy;

    void createNewFace(const Face& curFace, Figure &newFigure);

    Point2D convert3Dto2D(const Vector3D &point);

    void drawFaces(img::EasyImage &image, const Point2D &A, const Point2D &B, const Point2D &C, const Color& color, double &_1_over_Zg, const Face &face);

    void drawZ_bufferd_triangle(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, const img::Color &color,
                                double dzdx, double  dzdy, double Xg, double Yg, double _1_over_Zg, img::EasyImage & image);

public:

    int size;

    Z_driehoek() {}

    drieDdrawings triangulate(const drieDdrawings& oldDrawing);

    void doProjection();

    img::EasyImage drawTriangle(Figures3D &figures3D, const img::Color& background);

};


#endif //ENGINE_Z_BUFFER_DRIEHOEKEN_H