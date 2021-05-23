//
// Created by inte on 27.03.21.
//

#include "Z_buffer_driehoeken.h"

drieDdrawings Z_driehoek::triangulate(const drieDdrawings& oldDrawing) {

    /*
     * na het aanmaken van de figure, kan het zijn dat de figuren niet uit 3-hoeken bestaat
     * maar uit meer hoeken, we moeten dus de nieuwe figuren gaan aanmaken van de oude figure
     * En deze allemaal opgaan delen in 3-hoeken
     */
    for (auto & figure : oldDrawing.drieDfiguren) {
        Figure newFigure;

        //de points van de oude figuur blijven hetzelfde bij de nieuwe figuur als ook de kleur
        newFigure.points = figure.points;
        newFigure.ambientReflection = figure.ambientReflection;
        newFigure.diffuseRelfection = figure.diffuseRelfection;
        newFigure.speculaireReflection = figure.speculaireReflection;
        newFigure.reflectionCoefficient = figure.reflectionCoefficient;

        //we maken hier de nieuwe figuren aan
        for (const auto& currentFace : figure.faces) createNewFace(currentFace, newFigure);

        //we zetten alles 3D punten om naar 2D punten en bepalen ineens de 1/Z waarde
        for (int i = 0; i < newFigure.points.size(); ++i) {
            Point2D punt = convert3Dto2D(newFigure.points[i]);
            double _1_over_Z_ = 1/newFigure.points[i].z;
            newFigure.points2D.emplace_back(punt, _1_over_Z_);
        }

        //we voegen de figuur toe aan een lijst van nieuwe figuren
        newFaceAndPoints.drieDfiguren.push_back(newFigure);
    }
    return newFaceAndPoints;
}

void Z_driehoek::createNewFace(const Face& curFace, Figure &newFigure) {

    //als het al een driehoek is moet we het vlak niet nog eens opnieuw maken
    if (curFace.point_indexes.size() < 4) newFigure.faces.push_back(curFace);

    //anders delen we het vlak op in verschillende driehoek moet het eerst punt altijd hetzelfde "oorsprong" punt
    else{
        int startPoint = curFace.point_indexes[0];

        for (int i = 1; i < curFace.point_indexes.size()-1; ++i) {
            Face newFace;
            newFace.point_indexes.push_back(startPoint);
            newFace.point_indexes.push_back(curFace.point_indexes[i]);
            newFace.point_indexes.push_back(curFace.point_indexes[i+1]);
            newFigure.faces.push_back(newFace);

        }
    }


}

Point2D Z_driehoek::convert3Dto2D(const Vector3D &point) {

    Point2D point2D{};

    point2D.x = (point.x)/(-point.z);
    point2D.y = (point.y)/(-point.z);

    if (point2D.x < XMin) XMin = point2D.x;
    if (point2D.x > XMax) XMax = point2D.x;
    if (point2D.y < YMin) YMin = point2D.y;
    if (point2D.y > YMax) YMax = point2D.y;

    return point2D;
}

void Z_driehoek::doProjection() {

    vector<double> fitProp = moveImageToFitZb(size, XMin, XMax, YMin, YMax);
    /* fitprop vector information:
     *
     * fitProp[0] = imageX
     * fitProp[1] = imageY
     * fitProp[2] = d
     * fitProp[3] = dx
     * fitProp[4] = dy
     */

    imageX = fitProp[0];
    imageY = fitProp[1];
    d = fitProp[2];
    dx = fitProp[3];
    dy = fitProp[4];

    zBuffer = ZBuffer(fitProp[0], fitProp[1]);

}

void Z_driehoek::createLightMatrix(Vector3D &location, Light &light)
{
    //eyepointMatrix berekenen

    double r;
    double phi;
    double theta;

    Vector3D eyepoint = location;

    double xPow = pow(eyepoint.x, 2.0);
    double yPow = pow(eyepoint.y, 2.0);
    double zPow = pow(eyepoint.z, 2.0);
    double rPow = xPow + yPow + zPow;

    r = sqrt(rPow);
    theta = atan2(eyepoint.y, eyepoint.x);
    phi = acos(eyepoint.z/r);


    light.eye(1,1) = -sin(theta);
    light.eye(1,2) = -cos(theta)*cos(phi);
    light.eye(1,3) = cos(theta)*sin(phi);
    light.eye(2,1) = cos(theta);
    light.eye(2,2) = -sin(theta)*cos(phi);
    light.eye(2,3) = sin(theta)*sin(phi);
    light.eye(3,2) = sin(phi);
    light.eye(3,3) = cos(phi);
    light.eye(4,3) = -r;


}

Point2D Z_driehoek::lightConvert3Dto2D(const Vector3D &point, double &xMax, double &xMin, double &yMax, double &yMin)
{
    Point2D point2D{};

    point2D.x = (point.x)/(-point.z);
    point2D.y = (point.y)/(-point.z);

    if (point2D.x < xMin) xMin = point2D.x;
    if (point2D.x > xMax) xMax = point2D.x;
    if (point2D.y < yMin) yMin = point2D.y;
    if (point2D.y > yMax) yMax = point2D.y;

    return point2D;


}

void Z_driehoek::createZBufferLights(int shadowMask, Figures3D &figures3D)
{

    for (auto &light : lights3D)
    {
        // eerst nog punten transformeren
        Figures3D copyFigures = figures3D;

        //maal inverse
        for (auto &figure : copyFigures)
        {
            for (auto & point : figure.points)
            {
               point *= Matrix::inv(eyePointMatrix);
            }
        }

        createLightMatrix(light.second.location, light.second);

        double xMin = 0, xMax = 0, yMin = 0, yMax = 0;

        for (auto &figure : copyFigures)
        {
            figure.points2D.clear();

            for (auto & point : figure.points)
            {
                point *= light.second.eye;
                double _1_over_Z_ = 1/point.z;

                Point2D point2D = lightConvert3Dto2D(point, xMax, xMin, yMax, yMin);
                figure.points2D.emplace_back(point2D, _1_over_Z_);
            }
        }
        vector<double> fitProp  = moveImageToFitZb(shadowMask, xMin, xMax, yMin, yMax);
        /* fitprop vector information:
         *
         * fitProp[0] = imageX
         * fitProp[1] = imageY
         * fitProp[2] = d
         * fitProp[3] = dx
         * fitProp[4] = dy
         */

        light.second.d = fitProp[2];
        light.second.dx = fitProp[3];
        light.second.dy = fitProp[4];

        light.second.shadowMask = ZBuffer(fitProp[0], fitProp[1]);


        for (auto &figure : copyFigures)
        {
            for (auto & point : figure.points2D)
            {
                point.first.x = (point.first.x * light.second.d) + light.second.dx;
                point.first.y = (point.first.y * light.second.d) + light.second.dy;
            }

            for (int i = 0; i < figure.faces.size(); ++i)
            {

                //halen telkens een vlak uit de figuur en slagen die 3 punten op in een variable
                Point2D A_2D = figure.points2D[figure.faces[i].point_indexes[0] - 1].first;
                Point2D B_2D = figure.points2D[figure.faces[i].point_indexes[1] - 1].first;
                Point2D C_2D = figure.points2D[figure.faces[i].point_indexes[2] - 1].first;

                //bereken ook 1/Zg
                double _1_over_Zg = (figure.points2D[figure.faces[i].point_indexes[0]-1].second +
                                     figure.points2D[figure.faces[i].point_indexes[1]-1].second +
                                     figure.points2D[figure.faces[i].point_indexes[2]-1].second)/3;

                /*
                 * Gaan hier dzdx en dzdy berekenen
                 * Daarvoor hebben we vector W voor nodig en k
                 */
                Vector3D A_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[0] - 1]);
                Vector3D B_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[1] - 1]);
                Vector3D C_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[2] - 1]);

                Vector3D u = Vector3D::vector(B_3D-A_3D);
                Vector3D v = Vector3D::vector(C_3D-A_3D);
                Vector3D w = Vector3D::vector(Vector3D::cross(u, v));

                double k = w.x * A_3D.x + w.y * A_3D.y + w.z * A_3D.z;

                double dzdx = (w.x)/(-light.second.d*k);
                double dzdy = (w.y)/(-light.second.d*k);


                Vector3D n = Vector3D::vector(Vector3D::normalise(w));

                //we voegen dit toe aan een vlak want elk vlak heeft zijn eigen dzdx en dzdy
                figure.faces[i].dzdx_dzdy = make_pair(dzdx, dzdy);

                //bereken Ymin en Ymax van de afbeelding
                int Ymin = roundToInt(min(A_2D.y, min(B_2D.y, C_2D.y))+0.5);
                int Ymax = roundToInt(max(A_2D.y, max(B_2D.y, C_2D.y))-0.5);

                //berekenen het zwaarte punt van de afbeelding
                double Xg = (A_2D.x + B_2D.x + C_2D.x)/3;
                double Yg = (A_2D.y + B_2D.y + C_2D.y)/3;


                //maken een vector aan van de 3 lijnen waaruit de figuur bestaat
                vector<Line2D> linesFace = vector<Line2D>{Line2D(A_2D,B_2D, figure.ambientReflection), Line2D(A_2D,C_2D, figure.ambientReflection), Line2D(B_2D,C_2D, figure.ambientReflection)};

                //lopen over de hoogte van de figuur
                for (int Yi = Ymin; Yi <= Ymax; ++Yi) {

                    vector<double> XlXr;

                    //lopen over de lijnen die we hebben aangemaat
                    for (Line2D PQ : linesFace) {

                        //voorwaarde: (Yi-Yp)(Yi-Yq) <= 0 en Yp != Yq
                        if ((Yi-PQ.p1.y)*(Yi-PQ.p2.y) <= 0 && PQ.p1.y != PQ.p2.y){
                            double Xi = PQ.p1.x + (PQ.p2.x - PQ.p1.x)*((Yi-PQ.p1.y)/(PQ.p2.y - PQ.p1.y));
                            XlXr.push_back(Xi);
                        }
                    }

                    /*
                     * we hebben 2 punten gevonden die op de zelfde waarde van Y liggen, maar als we over deze horizontale gaan
                     * gelegen over de hoogte Yi snijden we 2 lijnen van de figuren, deze snijpunten zijn Xl en Xr
                     * de afstand tussen deze 2 is dan ook het aantal pixels die ingekleurd moeten worden op Yi
                     */
                    unsigned int Xl = roundToInt(min(XlXr[0], XlXr[1])+0.5);
                    unsigned int Xr = roundToInt(max(XlXr[0], XlXr[1])-0.5);


                    for (unsigned int i = Xl; i <= Xr; ++i)
                    {
                        double _1_over_Zi = _1_over_Zg + (i-Xg)*dzdx + (Yi-Yg)*dzdy;

                        if (light.second.shadowMask.zBuffer[i][Yi] > _1_over_Zi){
                            light.second.shadowMask.zBuffer[i][Yi] = _1_over_Zi;
                        }
                    }
                }
            }
        }
    }
}

img::EasyImage Z_driehoek::drawTriangle(Figures3D &figures3D, const img::Color &background, Lights3D& lights, int shadowMask) {

    //we voeren de doProjection functie uit
    doProjection();

    lights3D = lights;

    if (shadowMask != 0) createZBufferLights(shadowMask, figures3D);

    //maken een image aan met correcte formate (bepaald in doProjection)
    img::EasyImage image(roundToInt(imageX), roundToInt(imageY));
    //zetten de achtergrond op de juiste kleur
    image.clear(background);

//    auto start = high_resolution_clock::now();
//    cout << "Start Triangled" << endl << endl;

    //lopen over de figuren
    for (auto &figure : figures3D) {
//        auto start = high_resolution_clock::now();
//        cout << "Figuur aanmaken" << endl;
//

        /*
         * nu we de dx, dy en d waarde hebben bepaaldt in de doProjection functie kunnen we elk punt van de figuur
         * verschuiven naar de juiste plaats op de afbeelding zodat de figuur altijd gecentreerd staat
         */

        for (auto & point : figure.points2D) {
            point.first.x = (point.first.x * d) + dx;
            point.first.y = (point.first.y * d) + dy;
        }

        for (int i = 0; i < figure.faces.size(); ++i) {

            //halen telkens een vlak uit de figuur en slagen die 3 punten op in een variable
            Point2D A_2D = figure.points2D[figure.faces[i].point_indexes[0] - 1].first;
            Point2D B_2D = figure.points2D[figure.faces[i].point_indexes[1] - 1].first;
            Point2D C_2D = figure.points2D[figure.faces[i].point_indexes[2] - 1].first;

            //bereken ook 1/Zg
            double _1_over_Zg = (figure.points2D[figure.faces[i].point_indexes[0]-1].second +
                                figure.points2D[figure.faces[i].point_indexes[1]-1].second +
                                figure.points2D[figure.faces[i].point_indexes[2]-1].second)/3;

            /*
             * Gaan hier dzdx en dzdy berekenen
             * Daarvoor hebben we vector W voor nodig en k
             */
            Vector3D A_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[0] - 1]);
            Vector3D B_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[1] - 1]);
            Vector3D C_3D = Vector3D::point(figure.points[figure.faces[i].point_indexes[2] - 1]);

            Vector3D u = Vector3D::vector(B_3D-A_3D);
            Vector3D v = Vector3D::vector(C_3D-A_3D);
            Vector3D w = Vector3D::vector(Vector3D::cross(u, v));

            double k = w.x * A_3D.x + w.y * A_3D.y + w.z * A_3D.z;

            double dzdx = (w.x)/(-d*k);
            double dzdy = (w.y)/(-d*k);


            Vector3D n = Vector3D::vector(Vector3D::normalise(w));

            //we voegen dit toe aan een vlak want elk vlak heeft zijn eigen dzdx en dzdy
            figure.faces[i].dzdx_dzdy = make_pair(dzdx, dzdy);

            drawFaces(image, A_2D, B_2D, C_2D, figure, _1_over_Zg, figure.faces[i], n);
        }

//        // stoppen de chrono
//        auto stop = high_resolution_clock::now();
//        //eindtijd min starttijd == duration tijd
//        auto duration = duration_cast<milliseconds >(stop - start);
//
//        //heb dit gevonden op: https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
//        cout << "Figuur gemaakt : " << duration.count() << " miliseconds" << endl << endl;


    }

//    // stoppen de chrono
//    auto stop = high_resolution_clock::now();
//    //eindtijd min starttijd == duration tijd
//    auto duration = duration_cast<milliseconds >(stop - start);
//
//    //heb dit gevonden op: https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
//    cout << "Figuren gemaakt : " << duration.count() << " miliseconds" << endl << endl;

    return image;
}

void Z_driehoek::drawFaces(img::EasyImage &image, const Point2D &A, const Point2D &B, const Point2D &C, const Figure& figure, double &_1_over_Zg, const Face &face, const Vector3D &n) {



    //bereken Ymin en Ymax van de afbeelding
    int Ymin = roundToInt(min(A.y, min(B.y, C.y))+0.5);
    int Ymax = roundToInt(max(A.y, max(B.y, C.y))-0.5);

    //berekenen het zwaarte punt van de afbeelding
    double Xg = (A.x + B.x + C.x)/3;
    double Yg = (A.y + B.y + C.y)/3;

    //halen dzdx en dzdy uit de pair die deze waardes bevat
    double dzdx = face.dzdx_dzdy.first;
    double dzdy = face.dzdx_dzdy.second;

    //maken een vector aan van de 3 lijnen waaruit de figuur bestaat
    vector<Line2D> linesFace = vector<Line2D>{Line2D(A,B, figure.ambientReflection), Line2D(A,C, figure.ambientReflection), Line2D(B,C, figure.ambientReflection)};

    //lopen over de hoogte van de figuur
    for (int Yi = Ymin; Yi <= Ymax; ++Yi) {

        vector<double> XlXr;

        //lopen over de lijnen die we hebben aangemaat
        for (Line2D PQ : linesFace) {

            //voorwaarde: (Yi-Yp)(Yi-Yq) <= 0 en Yp != Yq
            if ((Yi-PQ.p1.y)*(Yi-PQ.p2.y) <= 0 && PQ.p1.y != PQ.p2.y){
                double Xi = PQ.p1.x + (PQ.p2.x - PQ.p1.x)*((Yi-PQ.p1.y)/(PQ.p2.y - PQ.p1.y));
                XlXr.push_back(Xi);
            }
        }

        /*
         * we hebben 2 punten gevonden die op de zelfde waarde van Y liggen, maar als we over deze horizontale gaan
         * gelegen over de hoogte Yi snijden we 2 lijnen van de figuren, deze snijpunten zijn Xl en Xr
         * de afstand tussen deze 2 is dan ook het aantal pixels die ingekleurd moeten worden op Yi
         */
        unsigned int Xl = roundToInt(min(XlXr[0], XlXr[1])+0.5);
        unsigned int Xr = roundToInt(max(XlXr[0], XlXr[1])-0.5);


        drawZ_bufferd_triangle(Xl, Yi, Xr, Yi, figure.ambientReflection, figure.diffuseRelfection, figure.speculaireReflection, figure.reflectionCoefficient, dzdx, dzdy, Xg, Yg, _1_over_Zg, image, n);

    }
}

void Z_driehoek::drawZ_bufferd_triangle(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, const Color &ambientReflection,
                                        const Color &diffuseReflection,const Color &specularReflection, double reflectionCoeff,
                                        double dzdx, double  dzdy, double Xg, double Yg, double _1_over_Zg, img::EasyImage & image, const Vector3D &n) {


    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());


    img::Color newColor;
    newColor.red = 0;
    newColor.blue = 0;
    newColor.green = 0;

    vector<double> colors{0.0,0.0,0.0};
    vector<double> ambientLight{0.0,0.0,0.0};
    vector<double> diffuseInfinityLight{0.0,0.0,0.0};


    for (const auto& light : lights3D){

        //maken ambientlight
        ambientLight[0]+=(light.second.ambientLight.red * ambientReflection.red);
        ambientLight[1]+=(light.second.ambientLight.green * ambientReflection.green);
        ambientLight[2]+=(light.second.ambientLight.blue * ambientReflection.blue);


        //maken infDiffuselight
        if(light.first == "infDiff" || light.first == "infSpec"){

            Vector3D l = -Vector3D::vector(light.second.ldVector);
            l = Vector3D::normalise(l);
            l*=eyePointMatrix;

            double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);

            if  (cosA > 0)
            {
                diffuseInfinityLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*cosA);
                diffuseInfinityLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*cosA);
                diffuseInfinityLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*cosA);
            }
        }

    }

//    createAmbientLight(ambientLight, ambientReflection);
//    createDiffuseInfinityLight(diffuseInfinityLight, diffuseReflection, n);

    for (unsigned int i = x0; i <= x1; ++i) {
        double _1_over_Zi = _1_over_Zg + (i-Xg)*dzdx + (y0-Yg)*dzdy;

        vector<double> diffusepointLight{0.0,0.0,0.0};
        vector<double> diffuseSpotLight{0.0,0.0,0.0};
        vector<double> speculaireLight{0.0,0.0,0.0};

        double z = 1/_1_over_Zi;
        double x = (-((double)i - dx) *z)/d ;
        double y = (-((double)y0 - dy)*z)/d;


        for (auto& light : lights3D){

            if(light.first == "fDiffPoint" || light.first == "fSpecPoint")
            {

                Vector3D location = Vector3D::point(light.second.location);
                location*=eyePointMatrix;
                Vector3D l = Vector3D::normalise(Vector3D::vector(location.x-x, location.y-y, location.z - z));

                double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);


                if (cosA > 0)
                {
                    double omzetting = 1-((1-cosA)/(1-cos(M_PI/2)));

                    if (!light.second.shadowMask.zBuffer.empty())
                    {
                        Vector3D point = Vector3D::point(x,y,z);
                        if (interpolatie(point, light.second)){
                            diffusepointLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*omzetting);
                            diffusepointLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*omzetting);
                            diffusepointLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*omzetting);
                        }
                    }

                    else
                    {
                        diffusepointLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*omzetting);
                        diffusepointLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*omzetting);
                        diffusepointLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*omzetting);
                    }


                }


                if (light.first == "fSpecPoint"){
                    Vector3D r = Vector3D::normalise(Vector3D::vector((2*n*cosA)-l));
                    Vector3D camera = Vector3D::normalise(-Vector3D::vector(x,y,z));

                    double cosB = (camera.x*r.x) + (camera.y*r.y) + (camera.z*r.z);

                    if (cosB > 0)
                    {
                        speculaireLight[0] += (light.second.specularLight.red*specularReflection.red*pow(cosB, reflectionCoeff));
                        speculaireLight[1] += (light.second.specularLight.green*specularReflection.green*pow(cosB, reflectionCoeff));
                        speculaireLight[2] += (light.second.specularLight.blue*specularReflection.blue*pow(cosB, reflectionCoeff));

                    }
                }
            }
//                createDiffusepointLight(diffusepointLight, diffuseReflection, n, _1_over_Zi, i, y0, light);


            else if(light.first == "fDiffSpot")
            {

                Vector3D location = Vector3D::point(light.second.location);
                location*=eyePointMatrix;
                Vector3D l = Vector3D::normalise(Vector3D::vector(location.x-x, location.y-y, location.z - z));

                double angle = light.second.spotAngle;
                angle = gradesToRad(angle);

                double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);

                if (cosA > cos(angle))
                {
                    double omzetting = 1-((1-cosA)/(1-cos(angle)));

                    diffusepointLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*omzetting);
                    diffusepointLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*omzetting);
                    diffusepointLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*omzetting);
                }
            }

            else if(light.first == "infSpec")
            {

                Vector3D l = -Vector3D::vector(light.second.ldVector);
                l = Vector3D::normalise(l);
                l*=eyePointMatrix;

                double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);

                Vector3D r = Vector3D::normalise(Vector3D::vector((2*n*cosA)-l));
                Vector3D camera = Vector3D::normalise(-Vector3D::vector(x,y,z));

                double cosB = (camera.x*r.x) + (camera.y*r.y) + (camera.z*r.z);

                if (cosB > 0)
                {
                    speculaireLight[0] += (light.second.specularLight.red*specularReflection.red*pow(cosB, reflectionCoeff));
                    speculaireLight[1] += (light.second.specularLight.green*specularReflection.green*pow(cosB, reflectionCoeff));
                    speculaireLight[2] += (light.second.specularLight.blue*specularReflection.blue*pow(cosB, reflectionCoeff));

                }
            }
        }


        vector<double> color{0.0,0.0,0.0};
        img::Color imgColor;
        imgColor.red = 0;
        imgColor.blue = 0;
        imgColor.green = 0;


        color[0] = ambientLight[0] + diffuseInfinityLight[0] + diffusepointLight[0] + diffuseSpotLight[0] + speculaireLight[0];
        color[1] = ambientLight[1] + diffuseInfinityLight[1] + diffusepointLight[1] + diffuseSpotLight[1] + speculaireLight[1];
        color[2] = ambientLight[2] + diffuseInfinityLight[2] + diffusepointLight[2] + diffuseSpotLight[2] + speculaireLight[2];

        if (color[0] > 1) color[0] = 1;
        if (color[1] > 1) color[1] = 1;
        if (color[2] > 1) color[2] = 1;

        imgColor.red = color[0]*255;
        imgColor.green = color[1]*255;
        imgColor.blue = color[2]*255;


//        newColor = createLight(ambientLight, diffuseInfinityLight, diffusepointLight, diffuseSpotLight, speculaireLight);

        newColor = imgColor;

        if (zBuffer.zBuffer[i][y0] > _1_over_Zi){
            zBuffer.zBuffer[i][y0] = _1_over_Zi;
            (image)(i, y0) = newColor;

        }
    }
}

bool Z_driehoek::interpolatie(Vector3D &point, Light &light) {

    point *= Matrix::inv(eyePointMatrix);

    point *= light.eye;

    double _1overZ_ = 1/point.z;

    Point2D point2D = convert3Dto2D(point);

    double x = (point2D.x *light.d) + light.dx;
    double y = (point2D.y *light.d) + light.dy;

    int LowX = floor(x);
    int HighX = floor(x);
    HighX++;

    int LowY = floor(y);
    int HighY = floor(y);
    HighY++;

    double alphaX = x - LowX;
    double alphaY = y - LowY;


    double _1_over_ZE = ((1-alphaX)*light.shadowMask.zBuffer[LowX][LowY]) + (alphaX*light.shadowMask.zBuffer[HighX][LowY]);
    double _1_over_ZF = ((1-alphaX)*light.shadowMask.zBuffer[LowX][HighY]) + (alphaX*light.shadowMask.zBuffer[HighX][HighY]);


    double _1overZi_ =  ((1-alphaY)*_1_over_ZE) + (alphaY*_1_over_ZF);

        if (abs(_1overZ_ - _1overZi_) < 0.00005 || _1overZ_ < _1overZi_) return true;
    return false;
}


//void Z_driehoek::createAmbientLight(vector<double>& ambientLight, const Color &ambientReflection) {
//
//    for (const auto& light : lights3D){
//        ambientLight[0]+=(light.second.ambientLight.red * ambientReflection.red);
//        ambientLight[1]+=(light.second.ambientLight.green * ambientReflection.green);
//        ambientLight[2]+=(light.second.ambientLight.blue * ambientReflection.blue);
//    }
//}
//
//void Z_driehoek::createDiffuseInfinityLight(vector<double> &diffuseInfinityLight, const Color &diffuseReflection, const Vector3D &n) {
//
//    for (const auto& light : lights3D){
//
//        if(light.first == "infDiff" || light.first == "infSpec"){
//
//            Vector3D l = -Vector3D::vector(light.second.ldVector);
//            l = Vector3D::normalise(l);
//            l*=eyePointMatrix;
//
//            double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);
//
//            if  (cosA > 0)
//            {
//                diffuseInfinityLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*cosA);
//                diffuseInfinityLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*cosA);
//                diffuseInfinityLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*cosA);
//            }
//        }
//    }
//}
//
//void Z_driehoek::createDiffusepointLight(vector<double> &diffusepointLight, const Color &diffuseReflection, const Vector3D &n, double _1_over_Zi,
//                                         unsigned int i, unsigned int y0, pair<string, Light> &light)
//{
//    double z = 1/_1_over_Zi;
//    double x = (-((double)i - dx) *z)/d ;
//    double y = (-((double)y0 - dy)*z)/d;
//
//    Vector3D location = Vector3D::point(light.second.location);
//    location*=eyePointMatrix;
//    Vector3D l = Vector3D::normalise(Vector3D::vector(location.x-x, location.y-y, location.z - z));
//
//
//    double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);
//    double omzetting = 1-((1-cosA)/(1-cos(M_PI/2)));
//
//    if (cosA > 0)
//    {
//        diffusepointLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*omzetting);
//        diffusepointLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*omzetting);
//        diffusepointLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*omzetting);
//    }
//}
//
//
//void Z_driehoek::createDiffuseSpotLight(vector<double> &diffusepointLight, const Color &diffuseReflection, const Vector3D &n,
//                                   double _1_over_Zi, unsigned int i, unsigned int y0, pair<string, Light> &light)
//{
//
//    double z = 1/_1_over_Zi;
//    double x = (-((double)i - dx) *z)/d ;
//    double y = (-((double)y0 - dy)*z)/d;
//
//    Vector3D location = Vector3D::point(light.second.location);
//    location*=eyePointMatrix;
//    Vector3D l = Vector3D::normalise(Vector3D::vector(location.x-x, location.y-y, location.z - z));
//
//    double angle = light.second.spotAngle;
//    angle = gradesToRad(angle);
//
//    double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);
//    double omzetting = 1-((1-cosA)/(1-cos(angle)));
//
//    if (cosA > cos(angle))
//    {
//        diffusepointLight[0] += (light.second.diffuseLight.red*diffuseReflection.red*omzetting);
//        diffusepointLight[1] += (light.second.diffuseLight.green*diffuseReflection.green*omzetting);
//        diffusepointLight[2] += (light.second.diffuseLight.blue*diffuseReflection.blue*omzetting);
//    }
//}
//
//
//void Z_driehoek::createSpeculaireLight(vector<double> &speculaireLight, const Color &speculaireReflection, const Vector3D &n,
//                                       double _1_over_Zi, unsigned int i, unsigned int y0, pair<string, Light> &light)
//{
//
//    double z = 1/_1_over_Zi;
//    double x = (-((double)i - dx) *z)/d ;
//    double y = (-((double)y0 - dy)*z)/d;
//
//    Vector3D location = Vector3D::point(light.second.location);
//    location*=eyePointMatrix;
//    Vector3D l = Vector3D::normalise(Vector3D::vector(location.x-x, location.y-y, location.z - z));
//
//    double cosA = (n.x*l.x) + (n.y*l.y) + (n.z*l.z);
//
//    Vector3D r = Vector3D::normalise(Vector3D::vector((2*n*cosA)-l));
//    Vector3D camera = Vector3D::normalise(-Vector3D::vector(x,y,z));
//
//    double cosB = (camera.x*r.x) + (camera.y*r.y) + (camera.z*r.z);
//
//    if (cosB > 0)
//    {
////        speculaireLight[0] += (light.second.specularLight.red*speculaireReflection.red*pow(cosB, reflCoeff));
////        speculaireLight[1] += (light.second.specularLight.green*speculaireReflection.green*pow(cosB, reflCoeff));
////        speculaireLight[2] += (light.second.specularLight.blue*speculaireReflection.blue*pow(cosB, reflCoeff));
//
//    }
//}
//
//
//
//img::Color Z_driehoek::createLight(vector<double> ambientLight, vector<double> diffuseInfinityLight,
//                                   vector<double> diffusepointLight, vector<double> diffuseSportLight,
//                                   vector<double> speculaireLight)
//{
//
//    vector<double> color{0.0,0.0,0.0};
//    img::Color imgColor;
//    imgColor.red = 0;
//    imgColor.blue = 0;
//    imgColor.green = 0;
//
//
//    color[0] = ambientLight[0] + diffuseInfinityLight[0] + diffusepointLight[0] + diffuseSportLight[0] + speculaireLight[0];
//    color[1] = ambientLight[1] + diffuseInfinityLight[1] + diffusepointLight[1] + diffuseSportLight[1] + speculaireLight[1];
//    color[2] = ambientLight[2] + diffuseInfinityLight[2] + diffusepointLight[2] + diffuseSportLight[2] + speculaireLight[2];
//
//
//    if (color[0] > 1) color[0] = 1;
//    if (color[1] > 1) color[1] = 1;
//    if (color[2] > 1) color[2] = 1;
//
//    imgColor.red = color[0]*255;
//    imgColor.green = color[1]*255;
//    imgColor.blue = color[2]*255;
//
//    return imgColor;
//
//}
//
