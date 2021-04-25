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
        newFigure.points = figure.points; newFigure.ambientReflection = figure.ambientReflection;

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
     * fitProp[2] = scale
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

img::EasyImage Z_driehoek::drawTriangle(Figures3D &figures3D, const img::Color &background, Lights3D& lights) {

    //we voeren de doProjection functie uit
    doProjection();

    lights3D = lights;

    //maken een image aan met correcte formate (bepaald in doProjection)
    img::EasyImage image(roundToInt(imageX), roundToInt(imageY));
    //zetten de achtergrond op de juiste kleur
    image.clear(background);

    //lopen over de figuren
    for (auto &figure : figures3D) {

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
            Vector3D A_3D = figure.points[figure.faces[i].point_indexes[0] - 1];
            Vector3D B_3D = figure.points[figure.faces[i].point_indexes[1] - 1];
            Vector3D C_3D = figure.points[figure.faces[i].point_indexes[2] - 1];

            Vector3D u = Vector3D::vector(B_3D.x - A_3D.x, B_3D.y - A_3D.y, B_3D.z - A_3D.z);
            Vector3D v = Vector3D::vector(C_3D.x - A_3D.x, C_3D.y - A_3D.y, C_3D.z - A_3D.z);
            Vector3D w = Vector3D::vector((u.y*v.z)-(u.z*v.y), (u.z*v.x)-(u.x*v.z), (u.x*v.y)-(u.y*v.x));

            double k = w.x * A_3D.x + w.y * A_3D.y + w.z * A_3D.z;

            double dzdx = (w.x)/(-d*k);
            double dzdy = (w.y)/(-d*k);

            //we voegen dit toe aan een vlak want elk vlak heeft zijn eigen dzdx en dzdy
            figure.faces[i].dzdx_dzdy = make_pair(dzdx, dzdy);

            drawFaces(image, A_2D, B_2D, C_2D, figure, _1_over_Zg, figure.faces[i]);

        }
    }

    return image;
}

void Z_driehoek::drawFaces(img::EasyImage &image, const Point2D &A, const Point2D &B, const Point2D &C, const Figure& figure, double &_1_over_Zg, const Face &face) {



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


        drawZ_bufferd_triangle(Xl, Yi, Xr, Yi, figure.ambientReflection, figure.diffuseRelfection, figure.speculaireReflection, figure.reflectionCoefficient, dzdx, dzdy, Xg, Yg, _1_over_Zg, image);

    }
}

void Z_driehoek::drawZ_bufferd_triangle(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, const Color &ambientReflection,
                                        const Color &diffuseReflection,const Color &specularReflection, double reflectionCoeff,
                                        double dzdx, double  dzdy, double Xg, double Yg, double _1_over_Zg, img::EasyImage & image) {


    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());

    img::Color newColor;

    vector<double> colors{0,0,0};

    for (auto light : lights3D) {

        vector<double> tempVector{light.ambientLight.red, light.ambientLight.green, light.ambientLight.blue};

        colors[0] += (tempVector[0]*ambientReflection.red);
        colors[1] += (tempVector[1]*ambientReflection.green);
        colors[2] += (tempVector[2]*ambientReflection.blue);

    }

    if (colors[0] > 1) colors[0] = 1;
    if (colors[1] > 1) colors[1] = 1;
    if (colors[2] > 1) colors[2] = 1;

    newColor.red = colors[0]*255;
    newColor.green = colors[1]*255;
    newColor.blue = colors[2]*255;

    for (unsigned int i = x0; i <= x1; ++i) {
        double _1_over_Zi = _1_over_Zg + (i-Xg)*dzdx + (y0-Yg)*dzdy;

        if (zBuffer.zBuffer[i][y0] > _1_over_Zi){
            zBuffer.zBuffer[i][y0] = _1_over_Zi;
            (image)(i, y0) = newColor;
        }
    }

}


