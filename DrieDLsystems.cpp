//
// Created by inte on 14.03.21.
//

#include "DrieDLsystems.h"
#include "drieDdrawings.h"

DrieDLsystems::DrieDLsystems(Matrix &rotMatrix) {rotatieMatrix = rotMatrix;}


void DrieDLsystems::parse3DL(const string &L3DinputFile) {

    std::filebuf fb;
    if (fb.open (L3DinputFile,std::ios::in))
    {
        std::istream is(&fb);
        LParser::LSystem3D input3DLSystem = LParser::LSystem3D(is);

        Angle = input3DLSystem.get_angle();
        Alfabet = input3DLSystem.get_alphabet();
        Initiator = input3DLSystem.get_initiator();
        Iterations = input3DLSystem.get_nr_iterations();

        for (char letter : Alfabet) {Replacements[letter] = input3DLSystem.get_replacement(letter);Draw[letter] = input3DLSystem.draw(letter);}
        fb.close();
    }
}

Figure DrieDLsystems::create3DLsystem() {

    //maken van de volledige string
    for (int i = 0; i < Iterations; ++i) { Initiator = createDrawString(Initiator);}
    pair<vector<Vector3D>,vector<Face>> returnPair = create3DPointsVector();
    Figure DLfigure;
    DLfigure.faces = returnPair.second;
    DLfigure.points = returnPair.first;

    return DLfigure;
}

string DrieDLsystems::createDrawString(string &initiator) {

    string newString;

    //lopen over elk character in de string
    for (char IniChar : initiator) {

        //als het een +, -, (, ), ^, &, \, / of | is (=operator), dan kan het al niet meer in het alfabet zitten
        //dus voegen we het gewoon toe aan de string
        if (isOperator3DL(IniChar)){newString += IniChar;}

        //geen stoch replacement rule want lengte vector is 1
        else if(Replacements.find(IniChar)->second.size() == 1){
            newString += Replacements.find(IniChar)->second[0].first;
        }

    }

    return newString;

}

pair<vector<Vector3D>,vector<Face>> DrieDLsystems::create3DPointsVector() {

    vector<Vector3D> points;
    vector<Face> faces;

    Vector3D punt1 = Vector3D::point(0,0,0);

    Vector3D H = Vector3D::vector(1,0,0);
    Vector3D L = Vector3D::vector(0,1,0);
    Vector3D U = Vector3D::vector(0,0,1);

    Angle = gradesToRad(Angle);


    for (char command : Initiator) {

        if (command == '+'){
            Vector3D Hnew = Vector3D::vector(H*cos(Angle) + L*sin((Angle)));
            Vector3D Lnew = Vector3D::vector(-H*sin((Angle)) + L*cos((Angle)));
            H=Hnew; L=Lnew;
        }
        else if (command == '-'){
            Vector3D Hnew = Vector3D::vector(H*cos((Angle)) + L*sin((-Angle)));
            Vector3D Lnew = Vector3D::vector(-H*sin((-Angle)) + L*cos((Angle)));
            H=Hnew; L=Lnew;
        }
        else if (command == '^'){
            Vector3D Hnew = Vector3D::vector(H*cos((Angle)) + U*sin((Angle)));
            Vector3D Unew = Vector3D::vector(-H*sin((Angle)) + U*cos((Angle)));
            H=Hnew; U=Unew;
        }
        else if (command == '&'){
            Vector3D Hnew = Vector3D::vector(H*cos((Angle)) + U*sin((-Angle)));
            Vector3D Unew = Vector3D::vector(-H*sin((-Angle)) + U*cos((Angle)));
            H=Hnew; U=Unew;
        }
        else if (command == 92){
            Vector3D Lnew = Vector3D::vector(L*cos((Angle)) + U*sin((-Angle)));
            Vector3D Unew = Vector3D::vector(-L*sin((-Angle)) + U*cos((Angle)));
            L=Lnew; U=Unew;
        }
        else if (command == '/'){
            Vector3D Lnew = Vector3D::vector(L*cos((Angle)) + U*sin((Angle)));
            Vector3D Unew = Vector3D::vector(-L*sin((Angle)) + U*cos((Angle)));
            L=Lnew; U=Unew;
        }
        else if (command == '|'){
            Vector3D Hnew = Vector3D::vector(H*cos(gradesToRad(180)) + L*sin(gradesToRad(180)));
            Vector3D Lnew = Vector3D::vector(-H*sin(gradesToRad(180)) + L*cos(gradesToRad(180)));
            H=Hnew; L=Lnew;
        }

        else if(command == '('){

            vector<Vector3D> newVec;
            newVec.push_back(punt1);
            newVec.push_back(H);
            newVec.push_back(L);
            newVec.push_back(U);
            stack.push_back(newVec);topStack += 1;
        }

        else if(command == ')'){

            vector<Vector3D> newVec = stack[topStack];
            punt1 = newVec[0];
            H = newVec[1];
            L = newVec[2];
            U = newVec[3];

            stack.pop_back();
            topStack-=1;
        }

        else{
            bool canDraw;
            canDraw = Draw.find(command)->second;
            Vector3D punt2 = punt1;
            punt1 += H;

            if (canDraw){
                Face newFace;
                points.push_back(punt2*rotatieMatrix);
                points.push_back(punt1*rotatieMatrix);
                newFace.point_indexes.push_back(points.size()-2);
                newFace.point_indexes.push_back(points.size()-1);
                faces.push_back(newFace);
            }
        }

    }


    return make_pair(points, faces);
}

