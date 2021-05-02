//
// Created by intev on 19/02/2021.
//

#include "draw2DLine.h"

//constructor
Line2D::Line2D(Point2D punt1, Point2D punt2, Color lineColor, double zp1, double zp2) {p1 = punt1; p2 = punt2; color = lineColor; z1 = zp1; z2 = zp2;}

//tekent de lijnen
img::EasyImage draw2DLines(const Lines2D &lines, const int size, const img::Color& background, bool zBuf){

    //vector die alle lijne die getekend moeten worden bijhoudt
    vector<Line2D> allLines;

    //voegt die lijnen toe aan de vector
    for (const auto & line : lines) { Line2D newLine = Line2D(line.p1, line.p2, line.color, line.z1, line.z2); allLines.push_back(newLine);}

    vector<double> fitProp = moveImageToFit(allLines, size);
    /* fitprop vector information:
     *
     * fitProp[0] = imageX
     * fitProp[1] = imageY
     * fitProp[2] = scale
     * fitProp[3] = dx
     * fitProp[4] = dy
     */

    //maakt image en geeft al lengte en breedte
    img::EasyImage image(roundToInt(fitProp[0]), roundToInt(fitProp[1]));

    //maken de zBuffer aan
    ZBuffer zBuffer(fitProp[0], fitProp[1]);


    //kleuren de achtergrond van de image in
    image.clear(background);

    //alle lijnen vermenigvuldigen met de scaleFactor en optellen met verschuivingsafstand
    for (auto & allLine : allLines) {
        //punt 1
        allLine.p1.x = (allLine.p1.x * fitProp[2]) + fitProp[3]; allLine.p1.y = (allLine.p1.y * fitProp[2]) + fitProp[4];
        //punt 2
        allLine.p2.x = (allLine.p2.x * fitProp[2]) + fitProp[3]; allLine.p2.y = (allLine.p2.y * fitProp[2]) + fitProp[4];

        //kleur van de lijn aanpassen
        img::Color lineColor(roundToInt(allLine.color.red*255),
                             roundToInt(allLine.color.green*255),
                             roundToInt(allLine.color.blue*255));


        if (zBuf){

            draw_zbuf_line(zBuffer, image, roundToInt(allLine.p1.x), roundToInt(allLine.p1.y), allLine.z1,
                           roundToInt(allLine.p2.x), roundToInt(allLine.p2.y), allLine.z2, lineColor);
        }

        else{
            //tekenen van de lijn
            image.draw_line(roundToInt(allLine.p1.x), roundToInt(allLine.p1.y),
                            roundToInt(allLine.p2.x), roundToInt(allLine.p2.y), lineColor);

        }
    }

    return image;
}
ZBuffer::ZBuffer(const int width, const int height) {
    for (int w = 0; w < width; ++w) {
        vector<double> row;

        for (int h = 0; h < height; ++h) {
            row.push_back(numeric_limits<double>::infinity());
        }
        zBuffer.push_back(row);
    }
}


void draw_zbuf_line(ZBuffer &zBuffer, img::EasyImage &image, unsigned int x0, unsigned int y0,
                    double z0, unsigned int x1, unsigned int y1, double z1,
                    const img::Color &color){

    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());

    int a; int aantalPixels;

    if(x0 == x1){
        if (y0 > y1){swap(y0,y1); swap(z0, z1);}
        a = y1-y0;
        aantalPixels = a+1;
        if (a==0) a=1;
    }

    else if(y0 == y1){
        if (x0 > x1){swap(x0,x1); swap(z0, z1);}
        a = x1-x0;
        if (a==0) a=1;
        aantalPixels = a+1;
    }

    else{

        int x = x0-x1;
        int y = y0-y1;

        //als ze negatief zijn maken we ze eerst positief
        if (x < 0){x*=-1;} if (y < 0){y*=-1;}

        //we nemen hier de lengte van y als aantalpixels
        if (x < y){a = y; aantalPixels = a+1;}

        //we nemen hier de lengte van x als aantalpixels
        else{ a = x; aantalPixels = a+1;}
    }

    for (double i = 0; i < aantalPixels; ++i) {

        double p = i/a;
        double _1_over_Zi = (p/z0) + ((1-p)/z1);

        unsigned int Xi = round(p*x0 + (1-p)*x1);
        unsigned int Yi = round(p*y0 + (1-p)*y1);

        if (zBuffer.zBuffer[Xi][Yi] > _1_over_Zi){
            zBuffer.zBuffer[Xi][Yi] = _1_over_Zi;
            (image)(Xi, Yi) = color;
        }
    }

}




