#include "easy_image.h"
#include "ini_configuration.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include "draw2DLine.h"
#include "_2DLsystems.h"
#include <time.h>
#include <algorithm>
#include <chrono>
#include "vector3d.h"
#include "drieDdrawings.h"
#include "Z_buffer_driehoeken.h"
using namespace std;
using namespace std::chrono;

img::EasyImage testFunctionDraw2DLine(const ini::Configuration &configuration){

    int size = configuration["ImageProperties"]["size"].as_int_or_die();
    vector<double> kleur = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();

    double p1X = configuration["LineProperties"]["p1X"].as_double_or_die();
    double p1Y = configuration["LineProperties"]["p1Y"].as_double_or_die();
    double p2X = configuration["LineProperties"]["p2X"].as_double_or_die();
    double p2Y = configuration["LineProperties"]["p2Y"].as_double_or_die();
    vector<double> achtergrondKleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    img::Color backColor;
    backColor.red = achtergrondKleur[0]*255;
    backColor.green = achtergrondKleur[1]*255;
    backColor.blue = achtergrondKleur[2]*255;


    Lines2D lijnen;
    Point2D p1{p1X, p1Y}, p2{p2X, p2Y};

    Color color1{kleur[0], kleur[1], kleur[2]};

    Line2D lijn1 = Line2D(p1, p2, color1);
    lijnen.push_back(lijn1);

    p1X = configuration["LineProperties"]["p1Xx"].as_double_or_die();
    p1Y = configuration["LineProperties"]["p1Yx"].as_double_or_die();
    p2X = configuration["LineProperties"]["p2Xx"].as_double_or_die();
    p2Y = configuration["LineProperties"]["p2Yx"].as_double_or_die();

    p1.x = p1X;
    p1.y = p1Y;

    p2.x = p2X;
    p2.y = p2Y;

    Line2D lijn2 = Line2D(p1, p2, color1);
    lijnen.push_back(lijn2);

    return draw2DLines(lijnen, size, backColor, false);


}

img::EasyImage generate_2DLSystem(const ini::Configuration &configuration){

    img::Color backColor;
    TweeDLSystem system;

    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> kleur = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
    vector<double> achtergrondKleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color lineColor{kleur[0], kleur[1], kleur[2]};

    backColor.red = achtergrondKleur[0]*255;
    backColor.green = achtergrondKleur[1]*255;
    backColor.blue = achtergrondKleur[2]*255;

    system.parse2DL(configuration["2DLSystem"]["inputfile"]);
    Lines2D drawLines = system.createDrawVector(lineColor);

    return draw2DLines(drawLines, size, backColor, false);
}

img::EasyImage generate_3Ddrawing(const ini::Configuration &configuration){


    img::Color backColor;
    Lines2D drawlines;

    drieDdrawings drawing;

    vector<double> eyepointVec = configuration["General"]["eye"].as_double_tuple_or_die();

    Vector3D eyepoint = Vector3D::point( eyepointVec[0], eyepointVec[1],eyepointVec[2]);

    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> achtergrondKleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    backColor.red = achtergrondKleur[0]*255;
    backColor.green = achtergrondKleur[1]*255;
    backColor.blue = achtergrondKleur[2]*255;

    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    drawing.eyePointTrans(eyepoint);
    for (int i = 0; i < nrFigures; ++i) {
        drawing.parse3Ddrawing(configuration, i, false);
    }

    drawlines = drawing.getLines();
    
    return draw2DLines(drawlines, size, backColor, false);

}

img::EasyImage generate_3DdrawingWithZbufferingWireFrame(const ini::Configuration &configuration){


    img::Color backColor;
    Lines2D drawlines;

    drieDdrawings drawing;

    vector<double> eyepointVec = configuration["General"]["eye"].as_double_tuple_or_die();

    Vector3D eyepoint = Vector3D::point( eyepointVec[0], eyepointVec[1],eyepointVec[2]);

    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> achtergrondKleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    backColor.red = achtergrondKleur[0]*255;
    backColor.green = achtergrondKleur[1]*255;
    backColor.blue = achtergrondKleur[2]*255;

    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    drawing.eyePointTrans(eyepoint);
    for (int i = 0; i < nrFigures; ++i) {
        drawing.parse3Ddrawing(configuration, i, false);
    }

    drawlines = drawing.getLines();

    return draw2DLines(drawlines, size, backColor, true);

}

img::EasyImage generate_3DdrawingWithZbuffering(const ini::Configuration &configuration){


    img::Color backColor;
    Lines2D drawlines;

    drieDdrawings drawing;

    vector<double> eyepointVec = configuration["General"]["eye"].as_double_tuple_or_die();

    Vector3D eyepoint = Vector3D::point( eyepointVec[0], eyepointVec[1],eyepointVec[2]);

    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> achtergrondKleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    backColor.red = achtergrondKleur[0]*255;
    backColor.green = achtergrondKleur[1]*255;
    backColor.blue = achtergrondKleur[2]*255;

    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    drawing.eyePointTrans(eyepoint);
    for (int i = 0; i < nrFigures; ++i) {
        drawing.parse3Ddrawing(configuration, i, true);
    }

    Z_driehoek zDriehoek;
    zDriehoek.size = size;

    drieDdrawings newFaceAndPoints = zDriehoek.triangulate(drawing);

    return zDriehoek.drawTriangle(newFaceAndPoints.drieDfiguren, backColor);

}



img::EasyImage generate_image(const ini::Configuration &configuration){

    //we veranderen de seed number altijd van het rand getal anders is de rand altijd hetzelfde
    srand (time(NULL));

    //2DLSystem
    if(configuration["General"]["type"].as_string_or_die() == "2DLSystem"){ return generate_2DLSystem(configuration);}

    else if (configuration["General"]["type"].as_string_or_die() == "Wireframe"){return generate_3Ddrawing(configuration);}

    else if (configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){return generate_3DdrawingWithZbufferingWireFrame(configuration);}

    else if (configuration["General"]["type"].as_string_or_die() == "ZBuffering"){return generate_3DdrawingWithZbuffering(configuration);}


        //hardcoded om lijnen te tekenen (test versie fucntie drawline)
    else{ return testFunctionDraw2DLine(configuration);}

}


int main(int argc, char const* argv[])
{

    //hiermee kan ik de runtime meten voor het maken van alle foto's
    auto start = high_resolution_clock::now();

    int retVal = 0;
    try
    {
        std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string fileName;
            while (std::getline(fileIn, fileName)) {
                args.push_back(fileName);
            }
        }
        for(std::string fileName : args)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }

    // stoppen de chrono
    auto stop = high_resolution_clock::now();
    //eindtijd min starttijd == duration tijd
    auto duration = duration_cast<milliseconds >(stop - start);

    //heb dit gevonden op: https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
    cout << "Time for making pictures: " << duration.count() << " miliseconds" << endl;
    return retVal;
}