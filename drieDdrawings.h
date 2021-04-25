//
// Created by intev on 5/03/2021.
//

#ifndef ENGINE_DRIEDDRAWINGS_H
#define ENGINE_DRIEDDRAWINGS_H
#include "vector3d.h"
#include <string>
#include <iostream>
#include "l_parser.h"
#include <fstream>
#include <set>
#include <vector>
#include "draw2DLine.h"
#include "ini_configuration.h"
#include "easy_image.h"
#include "DrieDLsystems.h"
using namespace std;


class Face{
public:
    //De indexen refereren naar
    //punten in de ‘points’ vector
    //van de Figure-klasse
    std::vector<int> point_indexes;
    pair<double, double> dzdx_dzdy;
};


class Figure{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    vector<pair<Point2D, double>> points2D;
    Color ambientReflection;
    Color diffuseRelfection;
    Color speculaireReflection;
    double reflectionCoefficient;
};

class Light
{
public:
    //de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de diffuse licht component
    Color specularLight;
};

class InfLight : public Light
{
public:
    //de richting waarin het licht schijnt
    Vector3D ldVector;
};

class PointLight : public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};

typedef std::list<Figure> Figures3D;
typedef std::list<Light> Lights3D;

class drieDdrawings {
public:

    Figures3D drieDfiguren;

    Lights3D lights;

    Matrix scaleFigure;
    Matrix rotateX;
    Matrix rotateY;
    Matrix rotateZ;
    Matrix translatie;
    Matrix rotatieMatrix;
    Matrix MatrixEyepoint;

    Lines2D lines;

    void parse3Ddrawing(const ini::Configuration &configuration, int figureNumber, bool zBufDriehoek, bool belichting);

    void createLights(const ini::Configuration &configuration, int nrLights, bool light);

    void eyePointTrans(const Vector3D &eyepoint);

    void rotateXmatrix(const double &angle);

    void rotateYmatrix(const double &angle);

    void rotateZmatrix(const double &angle);

    void translatieMatrix(Vector3D &centerPoint);

    void createRotatieMatrix(const double &angleX, const double &angleY, const double &angleZ, double &scale, Vector3D &centerPoint);

    void generateFractal(Figure& fig, Figures3D &fractal, const int nr_iterations, double scale);

    void createDrawVector(Figure &figure);

    void createMengerSponge(int nrIterations, bool zBufdriehoek, vector<double> &color);

    Matrix scaleMatrix(double &scale);

    Figure createCube(bool usage);

    Figure createTetrahedron();

    Figure createOctahedron();

    Figure createIcosahedron(bool usage = false);

    Figure createDodecahedron();

    Figure createCone(int n, double height);

    Figure createCylinder(int n, double height);

    Figure createSphere(int n);

    Figure createTorus(int n, int m, double r, double R);

    Figure createBuckyBall();

    Lines2D getLines() {return lines;}

    Point2D doProjection(const Vector3D &point, double d = 1);

};


#endif //ENGINE_DRIEDDRAWINGS_H
