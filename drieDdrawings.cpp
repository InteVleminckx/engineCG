//
// Created by intev on 5/03/2021.
//
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998
#include "drieDdrawings.h"

void drieDdrawings::parse3Ddrawing(const ini::Configuration &configuration, int figureNumber, bool zBufDriehoek) {

    string figurex = "Figure" + to_string(figureNumber);

    string type = configuration[figurex]["type"].as_string_or_die();
    double RotateX = configuration[figurex]["rotateX"].as_double_or_die();
    double RotateY = configuration[figurex]["rotateY"].as_double_or_die();
    double RotateZ = configuration[figurex]["rotateZ"].as_double_or_die();
    double scale = configuration[figurex]["scale"].as_double_or_die();
    vector<double> center = configuration[figurex]["center"].as_double_tuple_or_die();
    vector<double> color = configuration[figurex]["color"].as_double_tuple_or_die();

    Vector3D centerPoint;
    centerPoint.x = center[0]; centerPoint.y = center[1]; centerPoint.z = center[2];

    //maken vande rotatieMatrix
    createRotatieMatrix(RotateX, RotateY, RotateZ, scale, centerPoint);
    Figure newFigure;

    if (type == "LineDrawing"){

        int nrPoints = configuration[figurex]["nrPoints"].as_int_or_die();
        int nrLines = configuration[figurex]["nrLines"].as_int_or_die();

        for (int i = 0; i < nrPoints; ++i) {
            string pointnumber = "point" + to_string(i);
            vector<double> points = configuration[figurex][pointnumber].as_double_tuple_or_die();
            Vector3D vector3D = vector3D.point(points[0], points[1], points[2]);

            //hier al vermenigvuldigen
            vector3D *= rotatieMatrix * MatrixEyepoint;

            newFigure.points.push_back(vector3D);
        }

        for (int i = 0; i < nrLines; ++i) {
            string linenumber = "line" + to_string(i);
            vector<int> line = configuration[figurex][linenumber].as_int_tuple_or_die();
            Face newFace;
            newFace.point_indexes = line;
            newFigure.faces.push_back(newFace);
        }
    }

    else if (type == "Cube"){newFigure = createCube();}

    else if (type == "Cylinder"){
        int n = configuration[figurex]["n"].as_int_or_die();
        double height = configuration[figurex]["height"].as_double_or_die();
        newFigure = createCylinder(n, height);

    }

    else if (type == "Cone"){
        int n = configuration[figurex]["n"].as_int_or_die();
        double height = configuration[figurex]["height"].as_double_or_die();
        newFigure = createCone(n,height);
    }

    else if (type == "Sphere"){
        int n = configuration[figurex]["n"].as_int_or_die();
        newFigure = createSphere(n);

    }

    else if (type == "Torus"){
        int n = configuration[figurex]["n"].as_int_or_die();
        int m = configuration[figurex]["m"].as_int_or_die();
        double r = configuration[figurex]["r"].as_double_or_die();
        double R = configuration[figurex]["R"].as_double_or_die();

        newFigure = createTorus(n,m,r,R);

    }

    else if (type == "Tetrahedron"){newFigure = createTetrahedron();}

    else if (type == "Octahedron"){newFigure = createOctahedron();}

    else if (type == "Icosahedron"){newFigure = createIcosahedron();}

    else if (type == "Dodecahedron"){ newFigure = createDodecahedron();}

    else if (type == "3DLSystem"){
        Matrix volledigeMatrix = rotatieMatrix * MatrixEyepoint;
        DrieDLsystems system(volledigeMatrix);
        system.parse3DL(configuration[figurex]["inputfile"]);
        newFigure = system.create3DLsystem();
    }


    newFigure.color.red = color[0];
    newFigure.color.green = color[1];
    newFigure.color.blue = color[2];

    drieDfiguren.push_back(newFigure);

    if(!zBufDriehoek) createDrawVector(newFigure);


}

void drieDdrawings::eyePointTrans(const Vector3D &eyepoint) {
    double r;
    double phi;
    double theta;


    double xPow = pow(eyepoint.x, 2.0);
    double yPow = pow(eyepoint.y, 2.0);
    double zPow = pow(eyepoint.z, 2.0);
    double rPow = xPow + yPow + zPow;

    r = sqrt(rPow);
    theta = atan2(eyepoint.y, eyepoint.x);
    phi = acos(eyepoint.z/r);


    MatrixEyepoint(1,1) = -sin(theta);
    MatrixEyepoint(1,2) = -cos(theta)*cos(phi);
    MatrixEyepoint(1,3) = cos(theta)*sin(phi);
    MatrixEyepoint(2,1) = cos(theta);
    MatrixEyepoint(2,2) = -sin(theta)*cos(phi);
    MatrixEyepoint(2,3) = sin(theta)*sin(phi);
    MatrixEyepoint(3,2) = sin(phi);
    MatrixEyepoint(3,3) = cos(phi);
    MatrixEyepoint(4,3) = -r;

}

void drieDdrawings::scaleMatrix(double &scale) {
    scaleFigure(1,1) = scale;
    scaleFigure(2,2) = scale;
    scaleFigure(3,3) = scale;
}

void drieDdrawings::rotateXmatrix(const double &angle) {
    double radAngle = gradesToRad(angle);
    rotateX(2,2) = cos(radAngle);
    rotateX(2,3) = sin(radAngle);
    rotateX(3,2) = -sin(radAngle);
    rotateX(3,3) = cos(radAngle);

}

void drieDdrawings::rotateYmatrix(const double &angle) {

    double radAngle = gradesToRad(angle);
    rotateY(1,1) = cos(radAngle);
    rotateY(3,1) = sin(radAngle);
    rotateY(1,3) = -sin(radAngle);
    rotateY(3,3) = cos(radAngle);

}

void drieDdrawings::rotateZmatrix(const double &angle) {
    double radAngle = gradesToRad(angle);
    rotateZ(1,1) = cos(radAngle);
    rotateZ(1,2) = sin(radAngle);
    rotateZ(2,1) = -sin(radAngle);
    rotateZ(2,2) = cos(radAngle);

}

void drieDdrawings::translatieMatrix(Vector3D &centerPoint) {

    translatie(4,1) = centerPoint.x;
    translatie(4,2) = centerPoint.y;
    translatie(4,3) = centerPoint.z;
}

void drieDdrawings::createRotatieMatrix(const double &angleX, const double &angleY, const double &angleZ, double &scale, Vector3D &centerPoint) {

    scaleMatrix(scale);
    rotateXmatrix(angleX);
    rotateYmatrix(angleY);
    rotateZmatrix(angleZ);
    translatieMatrix(centerPoint);

//    scaleFigure.print(cout);cout << endl;
//    translatie.print(cout);cout << endl;
//    rotateX.print(cout);cout << endl;
//    rotateY.print(cout);cout << endl;
//    cout << endl;
//    MatrixEyepoint.print(cout); cout << endl;


    Matrix emptyMatrix;
    rotatieMatrix = emptyMatrix;
    rotatieMatrix = scaleFigure * rotateX * rotateY * rotateZ * translatie;

}

Point2D drieDdrawings::doProjection(const Vector3D &point, double d) {

    Point2D point2D{};

    point2D.x = (d*point.x)/(-point.z);
    point2D.y = (d*point.y)/(-point.z);

    return point2D;

}

void drieDdrawings::createDrawVector(Figure &figure) {

//    Lines2D lines;

//    for (auto &figure : drieDfiguren) {

       vector<Vector3D> points = figure.points;
       vector<Face> faces = figure.faces;
       vector<Point2D> puntenLijn;
       vector<double> zPoints;

        if (faces[0].point_indexes.size() == 2){
            for (int i = 0; i < faces.size(); ++i) {
                vector<int> punten = faces[i].point_indexes;
                for (int j = 0; j < punten.size(); ++j) {
                    Point2D punt = doProjection(points[punten[j]]);
                    zPoints.push_back(points[punten[j]].z);
                    puntenLijn.push_back(punt);
                }

                Line2D lijn = Line2D(puntenLijn[0], puntenLijn[1], figure.color, zPoints[zPoints.size()-2], zPoints[zPoints.size()-1]);
                lines.push_back(lijn);
                puntenLijn.clear();
            }
        }

        else{
            for (int i = 0; i < faces.size(); ++i) {
                Face vlak = faces[i];
                for (int j = 0; j < vlak.point_indexes.size(); ++j) {
                    vector<int> punten = vlak.point_indexes;

                    if (j < punten.size()-1){
                        Point2D punt1 = doProjection(points[punten[j]-1]);
                        Point2D punt2 = doProjection(points[punten[j+1]-1]);
                        double zPoint1 = points[punten[j]-1].z;
                        double zPoint2 = points[punten[j+1]-1].z;
                        Line2D line2D = Line2D(punt1, punt2, figure.color, zPoint1,zPoint2);
                        lines.push_back(line2D);
                    }
                    else {
                        Point2D punt1 = doProjection(points[punten[j]-1]);
                        Point2D punt2 = doProjection(points[punten[0]-1]);
                        double zPoint1 = points[punten[j]-1].z;
                        double zPoint2 = points[punten[0]-1].z;
                        Line2D line2D = Line2D(punt1, punt2, figure.color, zPoint1,zPoint2);
                        lines.push_back(line2D);
                    }
                }
            }
        }
//    }

//    return lines;
}

Figure drieDdrawings::createCube() {
    Figure cube;

    vector<vector<double>> punten;
    vector<double> rij1{1,-1,1,-1,1,-1,1,-1};
    vector<double> rij2{-1,1,1,-1,1,-1,-1,1};
    vector<double> rij3{-1,-1,1,1,-1,-1,1,1};
    punten.push_back(rij1);
    punten.push_back(rij2);
    punten.push_back(rij3);

    vector<vector<int>> vlakken;
    vector<int>rij1v{1,5,2,6,7,1};
    vector<int>rij2v{5,2,6,1,3,6};
    vector<int>rij3v{3,8,4,7,8,2};
    vector<int>rij4v{7,3,8,4,4,5};
    vlakken.push_back(rij1v);
    vlakken.push_back(rij2v);
    vlakken.push_back(rij3v);
    vlakken.push_back(rij4v);

    for (int i = 0; i < punten[0].size(); ++i) {
        Vector3D vector3D = Vector3D::point(punten[0][i], punten[1][i], punten[2][i]);
        vector3D *= rotatieMatrix * MatrixEyepoint;
        cube.points.push_back(vector3D);
    }

    for (int i = 0; i < vlakken[0].size(); ++i) {
        Face newFace;
        newFace.point_indexes = {vlakken[0][i],vlakken[1][i],vlakken[2][i],vlakken[3][i]};
        cube.faces.push_back(newFace);
    }
    return cube;
}

Figure drieDdrawings::createTetrahedron() {
    Figure tetrahedron;

    vector<vector<double>> punten;
    vector<double> rij1{1,-1,1,-1};
    vector<double> rij2{-1,1,1,-1};
    vector<double> rij3{-1,-1,1,1};
    punten.push_back(rij1);
    punten.push_back(rij2);
    punten.push_back(rij3);

    vector<vector<int>> vlakken;
    vector<int>rij1v{1,2,1,1};
    vector<int>rij2v{2,4,4,3};
    vector<int>rij3v{3,3,2,4};
    vlakken.push_back(rij1v);
    vlakken.push_back(rij2v);
    vlakken.push_back(rij3v);

    for (int i = 0; i < punten[0].size(); ++i) {
        Vector3D vector3D = vector3D.point(punten[0][i], punten[1][i], punten[2][i]);
        vector3D *= rotatieMatrix * MatrixEyepoint;
        tetrahedron.points.push_back(vector3D);
    }

    for (int i = 0; i < vlakken[0].size(); ++i) {
        Face newFace;
        newFace.point_indexes = {vlakken[0][i],vlakken[1][i],vlakken[2][i]};
        tetrahedron.faces.push_back(newFace);
    }
    return tetrahedron;
}

Figure drieDdrawings::createOctahedron() {
    Figure octahedron;

    vector<vector<double>> punten;
    vector<double> rij1{1,0,-1,0,0,0};
    vector<double> rij2{0,1,0,-1,0,0};
    vector<double> rij3{0,0,0,0,-1,1};
    punten.push_back(rij1);
    punten.push_back(rij2);
    punten.push_back(rij3);

    vector<vector<int>> vlakken;
    vector<int>rij1v{1,2,3,4,2,3,4,1};
    vector<int>rij2v{2,3,4,1,1,2,3,4};
    vector<int>rij3v{6,6,6,6,5,5,5,5};
    vlakken.push_back(rij1v);
    vlakken.push_back(rij2v);
    vlakken.push_back(rij3v);

    for (int i = 0; i < punten[0].size(); ++i) {
        Vector3D vector3D = vector3D.point(punten[0][i], punten[1][i], punten[2][i]);
        vector3D *= rotatieMatrix * MatrixEyepoint;
        octahedron.points.push_back(vector3D);
    }

    for (int i = 0; i < vlakken[0].size(); ++i) {
        Face newFace;
        newFace.point_indexes = {vlakken[0][i],vlakken[1][i],vlakken[2][i]};
        octahedron.faces.push_back(newFace);
    }
    return octahedron;
}

Figure drieDdrawings::createIcosahedron(bool usage) {
    Figure icosahedron;

    vector<double> x;
    vector<double> y;
    vector<double> z;

    for (int i = 1; i < 13; ++i) {
        if (i == 1){
            x.push_back(0);
            y.push_back(0);
            z.push_back(sqrt(5)/2);
        }

        else if (2 <= i && i <= 6){
            x.push_back(cos(((i-2)*2*PI)/5));
            y.push_back(sin(((i-2)*2*PI)/5));
            z.push_back(0.5);
        }

        else if (7 <= i && i <= 11){
            x.push_back(cos((PI/5)+((i-7)*2*PI)/5));
            y.push_back(sin((PI/5)+((i-7)*2*PI)/5));
            z.push_back(-0.5);
        }

        else if (i == 12){
            x.push_back(0);
            y.push_back(0);
            z.push_back(-(sqrt(5)/2));
        }
    }


    vector<vector<double>> punten;

    punten.push_back(x);
    punten.push_back(y);
    punten.push_back(z);

    vector<vector<int>> vlakken;
    vector<int>rij1v{1,1,1,1,1,2,3,3,4,4,5,5,6,6,2,12,12,12,12,12};
    vector<int>rij2v{2,3,4,5,6,7,7,8,8,9,9,10,10,11,11,8,9,10,11,7};
    vector<int>rij3v{3,4,5,6,2,3,8,4,9,5,10,6,11,2,7,7,8,9,10,11};
    vlakken.push_back(rij1v);
    vlakken.push_back(rij2v);
    vlakken.push_back(rij3v);

    for (int i = 0; i < punten[0].size(); ++i) {
        Vector3D vector3D = vector3D.point(punten[0][i], punten[1][i], punten[2][i]);
        if (!usage){
            vector3D *= rotatieMatrix * MatrixEyepoint;
        }
        icosahedron.points.push_back(vector3D);
    }

    for (int i = 0; i < vlakken[0].size(); ++i) {
        Face newFace;
        newFace.point_indexes = {vlakken[0][i],vlakken[1][i],vlakken[2][i]};
        icosahedron.faces.push_back(newFace);
    }
    return icosahedron;
}

Figure drieDdrawings::createDodecahedron() {
    Figure dodecahedron;
    Figure icosahedron = createIcosahedron(true);

    vector<double> x;
    vector<double> y;
    vector<double> z;

    for (int i = 0; i < icosahedron.faces.size(); ++i) {
        int point1 = icosahedron.faces[i].point_indexes[0]-1;
        int point2 = icosahedron.faces[i].point_indexes[1]-1;
        int point3 = icosahedron.faces[i].point_indexes[2]-1;

        double xPoint1 = icosahedron.points[point1].x;
        double yPoint1 = icosahedron.points[point1].y;
        double zPoint1 = icosahedron.points[point1].z;

        double xPoint2 = icosahedron.points[point2].x;
        double yPoint2 = icosahedron.points[point2].y;
        double zPoint2 = icosahedron.points[point2].z;

        double xPoint3 = icosahedron.points[point3].x;
        double yPoint3 = icosahedron.points[point3].y;
        double zPoint3 = icosahedron.points[point3].z;

        x.push_back((xPoint1 + xPoint2 + xPoint3)/3);
        y.push_back((yPoint1 + yPoint2 + yPoint3)/3);
        z.push_back((zPoint1 + zPoint2 + zPoint3)/3);

    }

    vector<vector<double>> punten;

    punten.push_back(x);
    punten.push_back(y);
    punten.push_back(z);

    vector<vector<int>> vlakken;
    vector<int>rij1v{1,1,2 ,3 ,4 ,5 ,20,20,19,18,17,16};
    vector<int>rij2v{2,6,8 ,10,12,14,19,15,13,11,9 ,7 };
    vector<int>rij3v{3,7,9 ,11,13,15,18,14,12,10,8 ,6 };
    vector<int>rij4v{4,8,10,12,14,6 ,17,13,11,9 ,7 ,15};
    vector<int>rij5v{5,2,3 ,4 ,5 ,1 ,16,19,18,17,16,20};
    vlakken.push_back(rij1v);
    vlakken.push_back(rij2v);
    vlakken.push_back(rij3v);
    vlakken.push_back(rij4v);
    vlakken.push_back(rij5v);


    for (int i = 0; i < punten[0].size(); ++i) {
        Vector3D vector3D = vector3D.point(punten[0][i], punten[1][i], punten[2][i]);
        vector3D *= rotatieMatrix * MatrixEyepoint;
        dodecahedron.points.push_back(vector3D);
    }

    for (int i = 0; i < vlakken[0].size(); ++i) {
        Face newFace;
        newFace.point_indexes = {vlakken[0][i],vlakken[1][i],vlakken[2][i],vlakken[3][i],vlakken[4][i]};
        dodecahedron.faces.push_back(newFace);
    }

    return dodecahedron;
}

Figure drieDdrawings::createCone(int n, double height) {

    Figure cone;


    for (int i = 0; i < n; ++i) {

        Vector3D punti = punti.point(cos((2*i*PI)/n), sin((2*i*PI)/n), 0);
        punti *= rotatieMatrix * MatrixEyepoint;
        cone.points.push_back(punti);
    }

    Vector3D bovenpunt = bovenpunt.point(0,0,height);
    bovenpunt *= rotatieMatrix * MatrixEyepoint;
    cone.points.push_back(bovenpunt);


    for (int i = 1; i < n+1; ++i) {
        Face vlak;
        vlak.point_indexes.push_back(i);

        if ((i+1)%n == 0){
            vlak.point_indexes.push_back((i+1)%(n+1));
        }
        else{
            vlak.point_indexes.push_back((i+1)%n);
        }

        vlak.point_indexes.push_back(n+1);
        cone.faces.push_back(vlak);
    }

    Face vlak;


    for (int i = -1; i < n-1; ++i) {
        vlak.point_indexes.push_back(n-(1+i));
    }

    cone.faces.push_back(vlak);

    return cone;
}

Figure drieDdrawings::createCylinder(int n, double height) {

    Figure cylinder;

    for (int i = 0; i < n; ++i) {

        Vector3D punti = punti.point(cos((2*i*PI)/n), sin((2*i*PI)/n), 0);
        punti *= rotatieMatrix * MatrixEyepoint;
        cylinder.points.push_back(punti);
    }

    for (int i = 0; i < n; ++i) {

        Vector3D punti = punti.point(cos((2*i*PI)/n), sin((2*i*PI)/n), height);
        punti *= rotatieMatrix * MatrixEyepoint;
        cylinder.points.push_back(punti);
    }


    for (int i = 1; i < n+1; ++i) {
        Face vlak;

        vlak.point_indexes.push_back(i);
        if ((i+1)%n == 0){
            vlak.point_indexes.push_back((i+1)%(n+1));
            vlak.point_indexes.push_back(((i+1)%(n+1))+n);
        }
        else{
            vlak.point_indexes.push_back((i+1)%n);
            vlak.point_indexes.push_back(((i+1)%n)+n);
        }

        vlak.point_indexes.push_back(i+n);
        cylinder.faces.push_back(vlak);
    }

    Face vlak;

    for (int i = -1; i < n-1; ++i) {
        vlak.point_indexes.push_back(n-(1+i));
    }
    cylinder.faces.push_back(vlak);
    vlak.point_indexes.clear();

    for (int i = -1; i < n-1; ++i) {
        vlak.point_indexes.push_back((n-(1+i))+n);
    }
    cylinder.faces.push_back(vlak);

    return cylinder;
}

Figure drieDdrawings::createSphere(int n) {
    Figure sphere;
    Figure icosahedron = createIcosahedron(true);

    for (int i = 0; i < n; ++i) {

        Figure tempFigure;

        for (int j = 0; j < icosahedron.faces.size(); ++j) {

            int Ap = icosahedron.faces[j].point_indexes[0];
            int Bp = icosahedron.faces[j].point_indexes[1];
            int Cp = icosahedron.faces[j].point_indexes[2];

            //ABC => ADE , BFD , CEF, DFE
            //ADE  = A, (A+B)/2, (A+C)/2
            //BFD  = B, (B+C)/2, (A+B)/2
            //CEF  = C, (A+C)/2, (B+C)/2
            //DFE  = (A+B)/2, (B+C)/2, (A+C)/2

            Vector3D A = icosahedron.points[Ap-1];
            Vector3D B = icosahedron.points[Bp-1];
            Vector3D C = icosahedron.points[Cp-1];
            Vector3D D = (A+B)/2;
            Vector3D E = (A+C)/2;
            Vector3D F = (B+C)/2;

            tempFigure.points.push_back(A); //1+j*6
            tempFigure.points.push_back(B); //2+j*6
            tempFigure.points.push_back(C); //3+j*6
            tempFigure.points.push_back(D); //4+j*6
            tempFigure.points.push_back(E); //5+j*6
            tempFigure.points.push_back(F); //6+j*6

            Face ADE;
            Face BFD;
            Face CEF;
            Face DFE;

            ADE.point_indexes.push_back(1+(j*6));
            ADE.point_indexes.push_back(4+(j*6));
            ADE.point_indexes.push_back(5+(j*6));

            BFD.point_indexes.push_back(2+(j*6));
            BFD.point_indexes.push_back(6+(j*6));
            BFD.point_indexes.push_back(4+(j*6));

            CEF.point_indexes.push_back(3+(j*6));
            CEF.point_indexes.push_back(5+(j*6));
            CEF.point_indexes.push_back(6+(j*6));

            DFE.point_indexes.push_back(4+(j*6));
            DFE.point_indexes.push_back(6+(j*6));
            DFE.point_indexes.push_back(5+(j*6));

            tempFigure.faces.push_back(ADE);
            tempFigure.faces.push_back(BFD);
            tempFigure.faces.push_back(CEF);
            tempFigure.faces.push_back(DFE);
        }
        icosahedron = tempFigure;
    }


    for (int i = 0; i < icosahedron.points.size(); ++i) {
        Vector3D punt = icosahedron.points[i];
        double xP = pow(punt.x, 2.0);
        double yP = pow(punt.y, 2.0);
        double zP = pow(punt.z, 2.0);
        double r = sqrt(xP + yP + zP);
        sphere.points.push_back(Vector3D().point((punt.x)/r,(punt.y)/r,(punt.z)/r));
        sphere.points[i] *= rotatieMatrix * MatrixEyepoint;
    }

    sphere.faces = icosahedron.faces;

    return sphere;
}

Figure drieDdrawings::createTorus(int n, int m, double r, double R) {

    Figure torus;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {

            double u = (2*i*PI)/n;
            double v = (2*j*PI)/n;
            double Xuv = (R + (r * cos(v)))*cos(u);
            double Yuv = (R + (r * cos(v)))*sin(u);
            double Zuv = r * sin(v);
            Vector3D punt1 = Vector3D::point(Xuv, Yuv, Zuv);
            punt1 *= rotatieMatrix * MatrixEyepoint;

            u = (2*(i+1)*PI)/n;
            v = (2*j*PI)/n;
            Xuv = (R + (r * cos(v)))*cos(u);
            Yuv = (R + (r * cos(v)))*sin(u);
            Zuv = r * sin(v);
            Vector3D punt2 = Vector3D::point(Xuv, Yuv, Zuv);
            punt2 *= rotatieMatrix * MatrixEyepoint;

            u = (2*(i+1)*PI)/n;
            v = (2*(j+1)*PI)/n;
            Xuv = (R + (r * cos(v)))*cos(u);
            Yuv = (R + (r * cos(v)))*sin(u);
            Zuv = r * sin(v);
            Vector3D punt3 = Vector3D::point(Xuv, Yuv, Zuv);
            punt3 *= rotatieMatrix * MatrixEyepoint;

            u = (2*i*PI)/n;
            v = (2*(j+1)*PI)/n;
            Xuv = (R + (r * cos(v)))*cos(u);
            Yuv = (R + (r * cos(v)))*sin(u);
            Zuv = r * sin(v);
            Vector3D punt4 = Vector3D::point(Xuv, Yuv, Zuv);
            punt4 *= rotatieMatrix * MatrixEyepoint;

            torus.points.push_back(punt1);
            torus.points.push_back(punt2);
            torus.points.push_back(punt3);
            torus.points.push_back(punt4);
        }
    }


    for (int i = 0; i < (torus.points.size())/4; ++i) {

        Face vlak;

        vlak.point_indexes.push_back((i*4)+1);
        vlak.point_indexes.push_back((i*4)+2);
        vlak.point_indexes.push_back((i*4)+3);
        vlak.point_indexes.push_back((i*4)+4);

        torus.faces.push_back(vlak);
    }

    return torus;
}

