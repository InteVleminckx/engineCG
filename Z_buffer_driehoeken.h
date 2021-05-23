//
// Created by inte on 27.03.21.
//

#ifndef ENGINE_Z_BUFFER_DRIEHOEKEN_H
#define ENGINE_Z_BUFFER_DRIEHOEKEN_H
#include <iostream>
#include "draw2DLine.h"
#include "drieDdrawings.h"
#include "easy_image.h"
#include <time.h>
#include <algorithm>
#include <chrono>

using namespace std::chrono;
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

    Lights3D lights3D;

    void createNewFace(const Face& curFace, Figure &newFigure);

    void createLightMatrix(Vector3D &location, Light &light);

    Point2D lightConvert3Dto2D(const Vector3D &point, double &xMax, double  &xMin, double &yMax, double &yMin);

    Point2D convert3Dto2D(const Vector3D &point);

    void createZBufferLights(int shadowMask, Figures3D &figures3D);

    void drawFaces(img::EasyImage &image, const Point2D &A, const Point2D &B, const Point2D &C, const Figure& Figure, double &_1_over_Zg, const Face &face, const Vector3D &n);

    void drawZ_bufferd_triangle(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, const Color &ambientReflection,
                                const Color &diffuseReflection,const Color &specularReflection, double reflectionCoeff,
                                double dzdx, double  dzdy, double Xg, double Yg, double _1_over_Zg, img::EasyImage & image, const Vector3D &n);

    bool interpolatie(Vector3D &point, Light &light);


//    void createAmbientLight(vector<double> &ambientLight, const Color &ambientReflection);
//
//    void createDiffuseInfinityLight(vector<double> &diffuseInfinityLight, const Color &diffuseReflection, const Vector3D &n);
//
//    void createDiffusepointLight(vector<double> &diffusepointLight, const Color &diffuseReflection, const Vector3D &n, double _1_over_Zi,
//                                 unsigned int i, unsigned int y0, pair<string, Light> &light);
//
//    void createDiffuseSpotLight(vector<double> &diffusepointLight, const Color &diffuseReflection, const Vector3D &n, double _1_over_Zi,
//                                unsigned int i, unsigned int y0, pair<string, Light> &light);
//
//    void createSpeculaireLight(vector<double> &speculaireLight, const Color &speculaireReflection, const Vector3D &n,
//                               double _1_over_Zi, unsigned int i, unsigned int y0, pair<string, Light> &light);
//
//    img::Color createLight(vector<double> ambientLight, vector<double> diffuseInfinityLight, vector<double> diffusepointLight, vector<double> diffuseSportLight, vector<double> speculaireLight);

public:

    int size;
    Matrix eyePointMatrix;

    Z_driehoek() {}

    drieDdrawings triangulate(const drieDdrawings& oldDrawing);

    void doProjection();

    img::EasyImage drawTriangle(Figures3D &figures3D, const img::Color& background, Lights3D& lights, int shadowMask = 0);

};


#endif //ENGINE_Z_BUFFER_DRIEHOEKEN_H
