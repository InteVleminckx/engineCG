//
// Created by intev on 21/02/2021.
//

#include "_2DLsystems.h"

//lijst het 2DL bestand uit
void TweeDLSystem::parse2DL(const string& L2DinputFile){

    std::filebuf fb;
    if (fb.open (L2DinputFile,std::ios::in))
    {
        std::istream is(&fb);
        LParser::LSystem2D input2DLSystem = LParser::LSystem2D(is);

        Angle = input2DLSystem.get_angle();
        Alfabet = input2DLSystem.get_alphabet();
        Initiator = input2DLSystem.get_initiator();
        StartingAngle = input2DLSystem.get_starting_angle();
        Iterations = input2DLSystem.get_nr_iterations();

        for (char letter : Alfabet) {Replacements[letter] = input2DLSystem.get_replacement(letter);Draw[letter] = input2DLSystem.draw(letter);}
        fb.close();
    }

}

//maakt een Lines2D aan, waarin alle lijnen zitten
Lines2D TweeDLSystem::createDrawVector(Color lineColor) {

    Lines2D drawVector; Point2D p1{}, p2{};

    double currentAngle = gradesToRad(StartingAngle);

    //we starten op punt 0,0
    p1.x = 0;  p1.y = 0;

    //maken van de volledige string
    for (int i = 0; i < Iterations; ++i) { Initiator = createDrawString(Initiator);}



    for (char command : Initiator) {

        //hoek verkleinen
        if (command == '-'){currentAngle -= gradesToRad(Angle);}

        //hoe vergroten
        else if(command == '+'){currentAngle += gradesToRad(Angle);}

        //We voegen een pair toe die dan de coor van het punt
        // en de hoek opslagen
        else if(command == '('){stack.emplace_back(p1, currentAngle);topStack += 1;}

        //we stellen p1 terug gelijk aan de waarde die we hebben opgeslagen
        //als ook voor de hoek en poppen het van de stack
        else if(command == ')'){p1 = stack[topStack].first; currentAngle = stack[topStack].second; stack.pop_back(); topStack-=1;}

        //lijn tekenen of punt verplaatsten
        else{
            bool canDraw;

            canDraw = Draw.find(command)->second;

            if (canDraw){
                //berekenen van het 2de punt
                p2.x = p1.x + cos(currentAngle); p2.y = p1.y + sin(currentAngle);

                //we maken de lijn aan
                Line2D drawLine = Line2D(p1, p2, lineColor);

                //voegen de lijn toe aan de vector // p1 wordt p2 we hebben onze "pen" verplaatst naar p2 en willen vanaf deze locatie verder tekenen
                drawVector.push_back(drawLine); p1 = p2;
            }
            //als we niet kunnen tekene wordt het punt gwn verplaatste
            else{p2.x = p1.x + cos(currentAngle); p2.y = p1.y + sin(currentAngle); p1 = p2;}
        }
    }
    return drawVector;
}

//maakt de string aan en doet de replacements
string TweeDLSystem::createDrawString(const string& Initia) {

    string newString;
    //lopen over elk character in de string
    for (char IniChar : Initia) {
        //als het een +, -, ( of ) is (=operator), dan kan het al niet meer in het alfabet zitten
        //dus voegen we het gewoon toe aan de string
        if (isOperator(IniChar)){newString += IniChar;}

        //als het de input wel tot het alfabet behoort
        //we voegen het replacement toe aan de string

        //geen stoch replacement rule want lengte vector is 1
        else if(Replacements.find(IniChar)->second.size() == 1){
            newString += Replacements.find(IniChar)->second[0].first;
        }

        else {
            double randompercentage = (double) rand() / RAND_MAX;
            if (randompercentage < Replacements.find(IniChar)->second[0].second) {newString += Replacements.find(IniChar)->second[0].first;}

            else {
                double currentPer = Replacements.find(IniChar)->second[0].second;

                for (int i = 0; i < Replacements.find(IniChar)->second.size()-1; ++i) {
                    if (currentPer < randompercentage < currentPer+Replacements.find(IniChar)->second[i+1].second){
                        newString += Replacements.find(IniChar)->second[i].first;
                        break;
                    }
                    currentPer += Replacements.find(IniChar)->second[i+1].second;
                }
            }
        }
    }
    return newString;
}
