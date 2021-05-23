//
// Created by intev on 5/03/2021.
//
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998
#include "drieDdrawings.h"

void drieDdrawings::parse3Ddrawing(const ini::Configuration &configuration, int figureNumber, bool zBufDriehoek, bool belichting) {

    string figurex = "Figure" + to_string(figureNumber);

    string type = configuration[figurex]["type"].as_string_or_die();
    double RotateX = configuration[figurex]["rotateX"].as_double_or_die();
    double RotateY = configuration[figurex]["rotateY"].as_double_or_die();
    double RotateZ = configuration[figurex]["rotateZ"].as_double_or_die();
    double scale = configuration[figurex]["scale"].as_double_or_die();
    vector<double> center = configuration[figurex]["center"].as_double_tuple_or_die();
    vector<double> ambientReflection;
    vector<double> diffuseReflection;
    vector<double> speculaireReflection;

    if (!belichting) ambientReflection = configuration[figurex]["color"].as_double_tuple_or_die();
    else
    {
        ambientReflection = configuration[figurex]["ambientReflection"].as_double_tuple_or_die();

        if (configuration[figurex]["diffuseReflection"].exists()){
            diffuseReflection = configuration[figurex]["diffuseReflection"].as_double_tuple_or_die();
        }

        if (configuration[figurex]["specularReflection"].exists()){
            speculaireReflection = configuration[figurex]["specularReflection"].as_double_tuple_or_die();
        }

    }



    string fractal = type.substr(0, 7);
    int nrIterations;
    double fractalScale;
    Figures3D fractalen;

    if (fractal == "Fractal"){
        nrIterations = configuration[figurex]["nrIterations"].as_int_or_die();
        fractalScale = configuration[figurex]["fractalScale"].as_double_or_die();
    }


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
            Vector3D vector3D = Vector3D::point(points[0], points[1], points[2]);

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

    else if (type == "Cube" || type == "FractalCube"){
        newFigure = createCube(false);
        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);
    }

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

    else if (type == "Tetrahedron" || type == "FractalTetrahedron"){
        newFigure = createTetrahedron();

        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);

    }

    else if (type == "Octahedron" || type == "FractalOctahedron"){
        newFigure = createOctahedron();
        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);
    }

    else if (type == "Icosahedron" || type == "FractalIcosahedron"){
        newFigure = createIcosahedron();
        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);
    }

    else if (type == "Dodecahedron" || type == "FractalDodecahedron"){
        newFigure = createDodecahedron();
        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);
    }

    else if ( type == "BuckyBall" || type == "FractalBuckyBall"){

        newFigure = createBuckyBall();
        if (fractal == "Fractal") generateFractal(newFigure, fractalen, nrIterations, fractalScale);
    }

    else if (type == "MengerSponge"){

        nrIterations = configuration[figurex]["nrIterations"].as_int_or_die();

        double refl = 0;

        if (configuration[figurex]["reflectionCoefficient"].exists())
        {
            refl = configuration[figurex]["reflectionCoefficient"];
        }

        createMengerSponge(nrIterations, zBufDriehoek, ambientReflection, diffuseReflection, speculaireReflection, refl);
        return;

    }

    else if (type == "3DLSystem"){
        Matrix volledigeMatrix = rotatieMatrix * MatrixEyepoint;
        DrieDLsystems system(volledigeMatrix);
        system.parse3DL(configuration[figurex]["inputfile"]);
        newFigure = system.create3DLsystem();
    }




    if(fractal != "Fractal"){
        newFigure.ambientReflection.red = ambientReflection[0];
        newFigure.ambientReflection.green = ambientReflection[1];
        newFigure.ambientReflection.blue = ambientReflection[2];

        if (configuration[figurex]["diffuseReflection"].exists()) {
            newFigure.diffuseRelfection.red = diffuseReflection[0];
            newFigure.diffuseRelfection.green = diffuseReflection[1];
            newFigure.diffuseRelfection.blue = diffuseReflection[2];
        }

        if (configuration[figurex]["specularReflection"].exists()) {
            newFigure.speculaireReflection.red = speculaireReflection[0];
            newFigure.speculaireReflection.green = speculaireReflection[1];
            newFigure.speculaireReflection.blue = speculaireReflection[2];
            newFigure.reflectionCoefficient = configuration[figurex]["reflectionCoefficient"];
        }

        drieDfiguren.push_back(newFigure);

        if(!zBufDriehoek && !newFigure.points.empty()) createDrawVector(newFigure);
    }

    else{
        for (auto &fractaal : fractalen) {


            fractaal.ambientReflection.red = ambientReflection[0];
            fractaal.ambientReflection.green = ambientReflection[1];
            fractaal.ambientReflection.blue = ambientReflection[2];

            if (configuration[figurex]["diffuseReflection"].exists()) {
                fractaal.diffuseRelfection.red = diffuseReflection[0];
                fractaal.diffuseRelfection.green = diffuseReflection[1];
                fractaal.diffuseRelfection.blue = diffuseReflection[2];
            }

            if (configuration[figurex]["specularReflection"].exists()) {
                fractaal.speculaireReflection.red = speculaireReflection[0];
                fractaal.speculaireReflection.green = speculaireReflection[1];
                fractaal.speculaireReflection.blue = speculaireReflection[2];
                fractaal.reflectionCoefficient = configuration[figurex]["reflectionCoefficient"];
            }
            drieDfiguren.push_back(fractaal);

            if(!zBufDriehoek) createDrawVector(fractaal);

        }
    }


}

void drieDdrawings::createLights(const ini::Configuration &configuration, int nrLights, bool light) {


    if (light){
        for (int i = 0; i < nrLights; ++i) {

            string light = "Light" + to_string(i);

            bool infinity, ambient, diffuse, specular, location, spotAngle, direction;

            bool infinityExist = configuration[light]["infinity"].exists();
            ambient = configuration[light]["ambientLight"].exists();
            diffuse = configuration[light]["diffuseLight"].exists();
            specular = configuration[light]["specularLight"].exists();
            location = configuration[light]["location"].exists();
            direction = configuration[light]["direction"].exists();
            spotAngle = configuration[light] ["spotAngle"].exists();

            if (infinityExist) infinity = configuration[light]["infinity"].as_bool_or_die();
            else infinity = false;

            string lightType;

            //geen diffuselicht dus ook al geen specular dus het is ambientlight
            if (!diffuse) lightType = "Ambient";

            /*
             * we hebben diffuse maar geen specular
             * we kunnen dit nog opsplitsen in:
             * infinity
             * finity point
             * finity spot
             */
            else if (diffuse && !specular)
            {
                //infinity
                if (infinity) lightType = "infDiff";

                else
                {
                    //spot
                    if (spotAngle) lightType = "fDiffSpot";

                    //point
                    else lightType = "fDiffPoint";
                }
            }

            /*
             * specularelight bevat ook 2 types
             * infinity diff
             * finity diff point
             */
            else if (specular)
            {
                //infinity
                if (infinity) lightType = "infSpec";

                //finity point
                else lightType = "fSpecPoint";
            }

            if (lightType == "Ambient") {
                Light newLight;
                newLight.type = lightType;

                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                lights.push_back(make_pair(lightType, newLight));

            }

            else if (lightType == "infDiff") {
                Light newLight;

                //type
                newLight.type = lightType;

                //ambientlight
                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                //diffuselight
                lightColors = configuration[light]["diffuseLight"].as_double_tuple_or_die();
                newLight.diffuseLight.red = lightColors[0];
                newLight.diffuseLight.green = lightColors[1];
                newLight.diffuseLight.blue = lightColors[2];

                //direction
                vector<double> directionLight = configuration[light]["direction"].as_double_tuple_or_die();
                newLight.ldVector = Vector3D::vector(directionLight[0],directionLight[1],directionLight[2]);

                lights.push_back(make_pair(lightType, newLight));
            }

            else if (lightType == "fDiffPoint") {
                Light newLight;

                //type
                newLight.type = lightType;

                //ambientlight
                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                //diffuselight
                lightColors = configuration[light]["diffuseLight"].as_double_tuple_or_die();
                newLight.diffuseLight.red = lightColors[0];
                newLight.diffuseLight.green = lightColors[1];
                newLight.diffuseLight.blue = lightColors[2];

                //location
                vector<double> locationLight = configuration[light]["location"].as_double_tuple_or_die();
                newLight.location = Vector3D::point(locationLight[0],locationLight[1],locationLight[2]);
                lights.push_back(make_pair(lightType, newLight));
            }

            else if (lightType == "fDiffSpot") {
                Light newLight;

                //type
                newLight.type = lightType;

                //ambientlight
                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                //diffuselight
                lightColors = configuration[light]["diffuseLight"].as_double_tuple_or_die();
                newLight.diffuseLight.red = lightColors[0];
                newLight.diffuseLight.green = lightColors[1];
                newLight.diffuseLight.blue = lightColors[2];

                //location
                vector<double> locationLight = configuration[light]["location"].as_double_tuple_or_die();
                newLight.location = Vector3D::point(locationLight[0],locationLight[1],locationLight[2]);

                //angle
                double angle = configuration[light]["spotAngle"].as_double_or_die();
                newLight.spotAngle = angle;

                lights.push_back(make_pair(lightType, newLight));
            }

            else if (lightType == "infSpec") {
                Light newLight;

                //type
                newLight.type = lightType;

                //ambientlight
                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                //diffuselight
                lightColors = configuration[light]["diffuseLight"].as_double_tuple_or_die();
                newLight.diffuseLight.red = lightColors[0];
                newLight.diffuseLight.green = lightColors[1];
                newLight.diffuseLight.blue = lightColors[2];

                //speculareLight
                lightColors = configuration[light]["specularLight"].as_double_tuple_or_die();
                newLight.specularLight.red = lightColors[0];
                newLight.specularLight.green = lightColors[1];
                newLight.specularLight.blue = lightColors[2];

                //direction
                vector<double> directionLight = configuration[light]["direction"].as_double_tuple_or_die();
                newLight.ldVector = Vector3D::vector(directionLight[0],directionLight[1],directionLight[2]);

                lights.push_back(make_pair(lightType, newLight));
            }

            else if (lightType == "fSpecPoint") {
                Light newLight;

                //type
                newLight.type = lightType;

                //ambientlight
                vector<double> lightColors = configuration[light]["ambientLight"].as_double_tuple_or_die();
                newLight.ambientLight.red = lightColors[0];
                newLight.ambientLight.green = lightColors[1];
                newLight.ambientLight.blue = lightColors[2];

                //diffuselight
                lightColors = configuration[light]["diffuseLight"].as_double_tuple_or_die();
                newLight.diffuseLight.red = lightColors[0];
                newLight.diffuseLight.green = lightColors[1];
                newLight.diffuseLight.blue = lightColors[2];

                //speculareLight
                lightColors = configuration[light]["specularLight"].as_double_tuple_or_die();
                newLight.specularLight.red = lightColors[0];
                newLight.specularLight.green = lightColors[1];
                newLight.specularLight.blue = lightColors[2];

                //location
                vector<double> locationLight = configuration[light]["location"].as_double_tuple_or_die();
                newLight.location = Vector3D::point(locationLight[0],locationLight[1],locationLight[2]);
                lights.push_back(make_pair(lightType, newLight));
            }
        }
    }

    else{
        Light newLight;
        newLight.ambientLight.red = 1.0;
        newLight.ambientLight.green = 1.0;
        newLight.ambientLight.blue = 1.0;

        lights.push_back(make_pair("normal", newLight));
    }

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

Matrix drieDdrawings::scaleMatrix(double &scale) {

    Matrix scaleMatrix;

    scaleMatrix(1,1) = scale;
    scaleMatrix(2,2) = scale;
    scaleMatrix(3,3) = scale;


    return scaleMatrix;

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

    scaleFigure = scaleMatrix(scale);
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

                Line2D lijn = Line2D(puntenLijn[0], puntenLijn[1], figure.ambientReflection, zPoints[zPoints.size()-2], zPoints[zPoints.size()-1]);
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
                        Line2D line2D = Line2D(punt1, punt2, figure.ambientReflection, zPoint1,zPoint2);
                        lines.push_back(line2D);
                    }
                    else {
                        Point2D punt1 = doProjection(points[punten[j]-1]);
                        Point2D punt2 = doProjection(points[punten[0]-1]);
                        double zPoint1 = points[punten[j]-1].z;
                        double zPoint2 = points[punten[0]-1].z;
                        Line2D line2D = Line2D(punt1, punt2, figure.ambientReflection, zPoint1,zPoint2);
                        lines.push_back(line2D);
                    }
                }
            }
        }
//    }

//    return lines;
}

Figure drieDdrawings::createCube(bool usage) {
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
        if (!usage){
            vector3D *= rotatieMatrix * MatrixEyepoint;
        }
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
        Vector3D vector3D = Vector3D::point(punten[0][i], punten[1][i], punten[2][i]);
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
        Vector3D vector3D = Vector3D::point(punten[0][i], punten[1][i], punten[2][i]);
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
        Vector3D vector3D = Vector3D::point(punten[0][i], punten[1][i], punten[2][i]);
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
        Vector3D vector3D = Vector3D::point(punten[0][i], punten[1][i], punten[2][i]);
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

        Vector3D punti = Vector3D::point(cos((2*i*PI)/n), sin((2*i*PI)/n), 0);
        punti *= rotatieMatrix * MatrixEyepoint;
        cone.points.push_back(punti);
    }

    Vector3D bovenpunt = Vector3D::point(0,0,height);
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

        Vector3D punti = Vector3D::point(cos((2*i*PI)/n), sin((2*i*PI)/n), 0);
        punti *= rotatieMatrix * MatrixEyepoint;
        cylinder.points.push_back(punti);
    }

    for (int i = 0; i < n; ++i) {

        Vector3D punti = Vector3D::point(cos((2*i*PI)/n), sin((2*i*PI)/n), height);
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

    //////////////////////////////////////////////////////////////
    reverse(vlak.point_indexes.begin(), vlak.point_indexes.end());
    //////////////////////////////////////////////////////////////

    cylinder.faces.push_back(vlak);
    vlak.point_indexes.clear();

    for (int i = -1; i < n-1; ++i) {
        vlak.point_indexes.push_back((n-(1+i))+n);
    }

    //////////////////////////////////////////////////////////////
    reverse(vlak.point_indexes.begin(), vlak.point_indexes.end());
    //////////////////////////////////////////////////////////////

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

            Vector3D A = Vector3D::point(icosahedron.points[Ap-1]);
            Vector3D B = Vector3D::point(icosahedron.points[Bp-1]);
            Vector3D C = Vector3D::point(icosahedron.points[Cp-1]);
            Vector3D D = Vector3D::point((A+B)/2);
            Vector3D E = Vector3D::point((A+C)/2);
            Vector3D F = Vector3D::point((B+C)/2);

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
        Vector3D punt = Vector3D::point(icosahedron.points[i]);
        double xP = pow(punt.x, 2.0);
        double yP = pow(punt.y, 2.0);
        double zP = pow(punt.z, 2.0);
        double r = sqrt(xP + yP + zP);
        sphere.points.push_back(Vector3D::point((punt.x)/r,(punt.y)/r,(punt.z)/r));
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

void drieDdrawings::generateFractal(Figure &fig, Figures3D &fractal, const int nr_iterations, double scale) {

    double scale1 = 1/scale;
    Matrix scaleM = scaleMatrix(scale1);
    fractal.push_back(fig);


    for (int i = 0; i < nr_iterations; ++i) {

        Figures3D tempFiguren;

        for (auto &fractaal : fractal) {

            Figure scaledFigure;

            for (auto & point : fractaal.points) scaledFigure.points.push_back(point*scaleM);

            for (int j = 0; j < fractaal.points.size(); ++j) {

                Vector3D Pi = Vector3D::point(fractaal.points[j]);
                Vector3D P_accent = Vector3D::point(scaledFigure.points[j]);
                Vector3D moveVector = Vector3D::vector(Pi-P_accent);

                Figure newFigure;

                for (auto & point : scaledFigure.points) newFigure.points.push_back(point+moveVector);

                newFigure.faces = fractaal.faces;

                tempFiguren.push_back(newFigure);

            }
        }

        fractal = tempFiguren;
    }

}

Figure drieDdrawings::createBuckyBall() {
    Figure buckyBall;
    Figure icosahedron = createIcosahedron(true);

    //vijfhoeken

    Vector3D p0 = Vector3D::point(icosahedron.points[0]);
    Vector3D p1 = Vector3D::point(icosahedron.points[1]);
    Vector3D p2 = Vector3D::point(icosahedron.points[2]);
    Vector3D p3 = Vector3D::point(icosahedron.points[3]);
    Vector3D p4 = Vector3D::point(icosahedron.points[4]);
    Vector3D p5 = Vector3D::point(icosahedron.points[5]);
    Vector3D p6 = Vector3D::point(icosahedron.points[6]);
    Vector3D p7 = Vector3D::point(icosahedron.points[7]);
    Vector3D p8 = Vector3D::point(icosahedron.points[8]);
    Vector3D p9 = Vector3D::point(icosahedron.points[9]);
    Vector3D p10 = Vector3D::point(icosahedron.points[10]);
    Vector3D p11 = Vector3D::point(icosahedron.points[11]);

    //vijfhoek1

    Vector3D middelpoint = p0;
    buckyBall.points.push_back((middelpoint + ((p1-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 1
    buckyBall.points.push_back((middelpoint + ((p2-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 2
    buckyBall.points.push_back((middelpoint + ((p3-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 3
    buckyBall.points.push_back((middelpoint + ((p4-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 4
    buckyBall.points.push_back((middelpoint + ((p5-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 5

    Face face;
    face.point_indexes = {1,2,3,4,5};
    buckyBall.faces.push_back(face);

    //vijfhoek2
    middelpoint = p1;
    buckyBall.points.push_back((middelpoint + ((p0-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 6
    buckyBall.points.push_back((middelpoint + ((p5-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 7
    buckyBall.points.push_back((middelpoint + ((p10-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 8
    buckyBall.points.push_back((middelpoint + ((p6-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 9
    buckyBall.points.push_back((middelpoint + ((p2-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 10

    face.point_indexes.clear();
    face.point_indexes = {6,7,8,9,10};
    buckyBall.faces.push_back(face);

    //vijfhoek3
    middelpoint = p2;
    buckyBall.points.push_back((middelpoint + ((p0-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 11
    buckyBall.points.push_back((middelpoint + ((p1-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 12
    buckyBall.points.push_back((middelpoint + ((p6-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 13
    buckyBall.points.push_back((middelpoint + ((p7-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 14
    buckyBall.points.push_back((middelpoint + ((p3-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 15

    face.point_indexes.clear();
    face.point_indexes = {11,12,13,14,15};
    buckyBall.faces.push_back(face);

    //vijfhoek4
    middelpoint = p3;
    buckyBall.points.push_back((middelpoint + ((p0-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 16
    buckyBall.points.push_back((middelpoint + ((p2-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 17
    buckyBall.points.push_back((middelpoint + ((p7-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 18
    buckyBall.points.push_back((middelpoint + ((p8-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 19
    buckyBall.points.push_back((middelpoint + ((p4-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 20

    face.point_indexes.clear();
    face.point_indexes = {16,17,18,19,20};
    buckyBall.faces.push_back(face);

    //vijfhoek5
    middelpoint = p4;
    buckyBall.points.push_back((middelpoint + ((p0-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 21
    buckyBall.points.push_back((middelpoint + ((p3-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 22
    buckyBall.points.push_back((middelpoint + ((p8-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 23
    buckyBall.points.push_back((middelpoint + ((p9-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 24
    buckyBall.points.push_back((middelpoint + ((p5-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 25

    face.point_indexes.clear();
    face.point_indexes = {21,22,23,24,25};
    buckyBall.faces.push_back(face);

    //vijfhoek6
    middelpoint = p5;
    buckyBall.points.push_back((middelpoint + ((p0-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 26
    buckyBall.points.push_back((middelpoint + ((p4-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 27
    buckyBall.points.push_back((middelpoint + ((p9-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 28
    buckyBall.points.push_back((middelpoint + ((p10-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 29
    buckyBall.points.push_back((middelpoint + ((p1-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 30

    face.point_indexes.clear();
    face.point_indexes = {26,27,28,29,30};
    buckyBall.faces.push_back(face);

    //vijfhoek7
    middelpoint = p6;
    buckyBall.points.push_back((middelpoint + ((p1-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 31
    buckyBall.points.push_back((middelpoint + ((p10-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 32
    buckyBall.points.push_back((middelpoint + ((p11-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 33
    buckyBall.points.push_back((middelpoint + ((p7-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 34
    buckyBall.points.push_back((middelpoint + ((p2-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 35

    face.point_indexes.clear();
    face.point_indexes = {31,32,33,34,35};
    buckyBall.faces.push_back(face);

    //vijfhoek8
    middelpoint = p7;
    buckyBall.points.push_back((middelpoint + ((p2-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 36
    buckyBall.points.push_back((middelpoint + ((p6-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 37
    buckyBall.points.push_back((middelpoint + ((p11-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 38
    buckyBall.points.push_back((middelpoint + ((p8-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 39
    buckyBall.points.push_back((middelpoint + ((p3-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 40

    face.point_indexes.clear();
    face.point_indexes = {36,37,38,39,40};
    buckyBall.faces.push_back(face);

    //vijfhoek9
    middelpoint = p8;
    buckyBall.points.push_back((middelpoint + ((p3-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 41
    buckyBall.points.push_back((middelpoint + ((p7-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 42
    buckyBall.points.push_back((middelpoint + ((p11-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 43
    buckyBall.points.push_back((middelpoint + ((p9-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 44
    buckyBall.points.push_back((middelpoint + ((p4-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 45

    face.point_indexes.clear();
    face.point_indexes = {41,42,43,44,45};
    buckyBall.faces.push_back(face);

    //vijfhoek10
    middelpoint = p9;
    buckyBall.points.push_back((middelpoint + ((p4-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 46
    buckyBall.points.push_back((middelpoint + ((p8-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 47
    buckyBall.points.push_back((middelpoint + ((p11-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 48
    buckyBall.points.push_back((middelpoint + ((p10-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 49
    buckyBall.points.push_back((middelpoint + ((p5-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 50

    face.point_indexes.clear();
    face.point_indexes = {46,47,48,49,50};
    buckyBall.faces.push_back(face);

    //vijfhoek11
    middelpoint = p10;
    buckyBall.points.push_back((middelpoint + ((p1-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 51
    buckyBall.points.push_back((middelpoint + ((p5-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 52
    buckyBall.points.push_back((middelpoint + ((p9-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 53
    buckyBall.points.push_back((middelpoint + ((p11-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 54
    buckyBall.points.push_back((middelpoint + ((p6-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 55

    face.point_indexes.clear();
    face.point_indexes = {51,52,53,54,55};
    buckyBall.faces.push_back(face);

    //vijfhoek12
    middelpoint = p11;
    buckyBall.points.push_back((middelpoint + ((p10-middelpoint)/3))*rotatieMatrix* MatrixEyepoint); //point 56
    buckyBall.points.push_back((middelpoint + ((p9-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 57
    buckyBall.points.push_back((middelpoint + ((p8-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 58
    buckyBall.points.push_back((middelpoint + ((p7-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 59
    buckyBall.points.push_back((middelpoint + ((p6-middelpoint)/3))*rotatieMatrix * MatrixEyepoint); //point 60

    face.point_indexes.clear();
    face.point_indexes = {56,57,58,59,60};
    buckyBall.faces.push_back(face);

    //zeshoek1
    face.point_indexes.clear();
    face.point_indexes = {1,6,10,12,11,2};
    buckyBall.faces.push_back(face);

    //zeshoek2
    face.point_indexes.clear();
    face.point_indexes = {2,11,15,17,16,3};
    buckyBall.faces.push_back(face);

    //zeshoek3
    face.point_indexes.clear();
    face.point_indexes = {4,3,16,20,22,21};
    buckyBall.faces.push_back(face);

    //zeshoek4
    face.point_indexes.clear();
    face.point_indexes = {5,4,21,25,27,26};
    buckyBall.faces.push_back(face);

    //zeshoek5
    face.point_indexes.clear();
    face.point_indexes = {7,6,1,5,26,30};
    buckyBall.faces.push_back(face);

    //zeshoek6
    face.point_indexes.clear();
    face.point_indexes = {30,29,52,51,8,7};
    buckyBall.faces.push_back(face);

    //zeshoek7
    face.point_indexes.clear();
    face.point_indexes = {9,8,51,55,32,31};
    buckyBall.faces.push_back(face);

    //zeshoek8
    face.point_indexes.clear();
    face.point_indexes = {9,31,35,13,12,10};
    buckyBall.faces.push_back(face);

    //zeshoek9
    face.point_indexes.clear();
    face.point_indexes = {13,35,34,37,36,14};
    buckyBall.faces.push_back(face);

    //zeshoek10
    face.point_indexes.clear();
    face.point_indexes = {17,15,14,36,40,18};
    buckyBall.faces.push_back(face);

    //zeshoek11
    face.point_indexes.clear();
    face.point_indexes = {42,41,19,18,40,39};
    buckyBall.faces.push_back(face);

    //zeshoek12
    face.point_indexes.clear();
    face.point_indexes = {23,22,20,19,41,45};
    buckyBall.faces.push_back(face);

    //zeshoek13
    face.point_indexes.clear();
    face.point_indexes = {44,47,46,24,23,45};
    buckyBall.faces.push_back(face);

    //zeshoek14
    face.point_indexes.clear();
    face.point_indexes = {50,28,27,25,24,46};
    buckyBall.faces.push_back(face);

    //zeshoek15
    face.point_indexes.clear();
    face.point_indexes = {50,49,53,52,29,28};
    buckyBall.faces.push_back(face);

    //zeshoek16
    face.point_indexes.clear();
    face.point_indexes = {32,55,54,56,60,33};
    buckyBall.faces.push_back(face);

    //zeshoek17
    face.point_indexes.clear();
    face.point_indexes = {53,49,48,57,56,54};
    buckyBall.faces.push_back(face);

    //zeshoek18
    face.point_indexes.clear();
    face.point_indexes = {57,48,47,44,43,58};
    buckyBall.faces.push_back(face);

    //zeshoek19
    face.point_indexes.clear();
    face.point_indexes = {58,43,42,39,38,59};
    buckyBall.faces.push_back(face);

    //zeshoek20
    face.point_indexes.clear();
    face.point_indexes = {33,60,59,38,37,34};
    buckyBall.faces.push_back(face);

    return buckyBall;
}

void drieDdrawings::createMengerSponge(int nrIterations, bool zBufdriehoek, vector<double> &color, vector<double> diffuus, vector<double> spec, double reflCoeff) {
    Figure cube = createCube(false);

    double scale =  1.0/3.0;

    Matrix scaleM = scaleMatrix(scale);

    Figures3D figuren;
    figuren.push_back(cube);



    if (nrIterations == 0){
        cube.ambientReflection.red = color[0];
        cube.ambientReflection.green = color[1];
        cube.ambientReflection.blue = color[2];

        if(!diffuus.empty()){
            cube.diffuseRelfection.red = diffuus[0];
            cube.diffuseRelfection.green = diffuus[1];
            cube.diffuseRelfection.blue = diffuus[2];
        }

        if (!spec.empty()){
            cube.speculaireReflection.red = spec[0];
            cube.speculaireReflection.green = spec[1];
            cube.speculaireReflection.blue = spec[2];
            cube.reflectionCoefficient = reflCoeff;
        }



        drieDfiguren.push_back(cube);
        if(!zBufdriehoek) createDrawVector(cube);
        return;
    }


    for (int i = 0; i < nrIterations; ++i) {

        Figures3D newFigures;

        for(Figure &figuur : figuren){

            for (int j = 0; j < figuur.faces.size()-2; ++j) {

                int Ai = figuur.faces[j].point_indexes[0]-1;
                int Bi = figuur.faces[j].point_indexes[1]-1;
                int Ci = figuur.faces[j].point_indexes[2]-1;
                int Di = figuur.faces[j].point_indexes[3]-1;

                Vector3D A = Vector3D::point(figuur.points[Ai]);
                Vector3D B = Vector3D::point(figuur.points[Bi]);
                Vector3D C = Vector3D::point(figuur.points[Ci]);
                Vector3D D = Vector3D::point(figuur.points[Di]);
                Vector3D E = Vector3D::point(A + ((B-A)/3));
                Vector3D L = Vector3D::point(A - ((A-D)/3));
                Vector3D G = Vector3D::point(B + ((C-B)/3));
                Vector3D J = Vector3D::point(D - ((D-C)/3));

                Vector3D Ac = Vector3D::point(figuur.points[Ai]*scaleM);
                Vector3D Bc = Vector3D::point(figuur.points[Bi]*scaleM);
                Vector3D Cc = Vector3D::point(figuur.points[Ci]*scaleM);
                Vector3D Dc = Vector3D::point(figuur.points[Di]*scaleM);


                Vector3D mVec1 = Vector3D::vector(A-Ac);
                Vector3D mVec2 = Vector3D::vector(E-Ac);
                Vector3D mVec3 = Vector3D::vector(B-Bc);
                Vector3D mVec4 = Vector3D::vector(L-Ac);
                Vector3D mVec5 = Vector3D::vector(G-Bc);
                Vector3D mVec6 = Vector3D::vector(D-Dc);
                Vector3D mVec7 = Vector3D::vector(J-Dc);
                Vector3D mVec8 = Vector3D::vector(C-Cc);

                //voorkomt dubbele kubussen
                vector<Vector3D> moveVectors1{mVec1, mVec2, mVec3, mVec4, mVec5, mVec6, mVec7, mVec8};
                vector<Vector3D> moveVectors2{mVec2,mVec7};

                vector<Vector3D> moveVectors;

                if (j == 0 || j == 2) moveVectors = moveVectors1;
                else moveVectors = moveVectors2;

                for(Vector3D &moveVector : moveVectors){

                    Figure newFigure;

                    for (auto & point : figuur.points) newFigure.points.push_back((point*scaleM)+moveVector);

                    newFigure.faces = figuur.faces;

                    if (i == nrIterations-1){
                        newFigure.ambientReflection.red = color[0];
                        newFigure.ambientReflection.green = color[1];
                        newFigure.ambientReflection.blue = color[2];

                        if(!diffuus.empty()){
                            newFigure.diffuseRelfection.red = diffuus[0];
                            newFigure.diffuseRelfection.green = diffuus[1];
                            newFigure.diffuseRelfection.blue = diffuus[2];
                        }

                        if (!spec.empty()){
                            newFigure.speculaireReflection.red = spec[0];
                            newFigure.speculaireReflection.green = spec[1];
                            newFigure.speculaireReflection.blue = spec[2];
                            newFigure.reflectionCoefficient = reflCoeff;
                        }



                        drieDfiguren.push_back(newFigure);

                        if(!zBufdriehoek) createDrawVector(newFigure);
                    }

                    else newFigures.push_back(newFigure);
                }
            }
        }
        figuren = newFigures;
    }
}



