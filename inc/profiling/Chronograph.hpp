#pragma once

#include"0_General.hpp"

#include <string>

#include <sstream>

#include <iostream> 


auto get_Chronograph(
    simpulationVariable& simu,
    clockstruct& chron
){


    std::string A;
    std::string t = ", ";
    std::string n = "\n";

    chron.beginNew.stop();
    chron.beginMainloop.stop();


    static auto s = [&](auto a){ return std::to_string(a);};
    static auto ws = [&](auto a){ std::ostringstream ss; ss << a; return ss.str();};

    
    if (simu.loop == 1)
    {
        A += "MainCount, ";
        A += "ALL, , ";
        A += "SolvePossionEquation, , ";
        A += "# iters, error, ";
        A += "Con & Dif, ";
        A += "T1 to T3, ";
        A += "check SS (L2), ";
        A += "T3 to T0, ";
        A += "BC";
        A += n;
    }

        A +=  ( s(simu.loop)+t);                                 //MainCount
        A +=  ( s(chron.beginNew.elapsedTime())+t+t);            //ALL
        A +=  ( s(chron.p.elapsedTime())+t+t);             //SolvePossionEquation
        A +=  ( s(simu.iters)+t);                                //#iters 
        A +=  ( ws(simu.error)+t);                               //error
        A +=  ( s(chron.convectionDifussion.elapsedTime())+t);   //Convection_Difussion
        A +=  ( s(chron.updateT1toT3.elapsedTime())+t);          //update_UandF_seq
        A +=  ( s(chron.checkL2norm.elapsedTime())+t) ;          //checkL2norm
        A +=  ( s(chron.updateT3toT0.elapsedTime())+t);          //T3toT0
        A +=  ( s(chron.BC.elapsedTime())+t);                    //BC
        A +=  n;
                        // Time_file << chron.beginMainloop.elapsedTime()<< ", ";    //Mainloop


    chron.beginNew.stop();
    chron.beginMainloop.stop();

    return A;
}



auto get_Chronograph_DAT(
    simpulationVariable& simu,
    clockstruct& chron
){
    std::string A;
    std::string t = ", ";
    std::string tab = " ";
    std::string n = "\n";
    // * ---------------------------------------init 
    if (simu.loop == 1)
    {
        std::vector<std::string> variables;
        variables.push_back("simulation time");
        variables.push_back("Calculation time (pressure)");
        variables.push_back("#iters");

                    A +=  "TITLE     = \"\"\n";
                    A += "VARIABLES = \"";
                    A += variables.at(0);
                    A += "\",\"";
                    A += variables.at(1);
                    A += "\",\"";
                    A += variables.at(2);
                    A += "\"\n";
                    A += "ZONE T=\"";
                    A += simu.ZONE();
                    A += "\"\n";
    }
    // * ---------------------------------------init 


    static auto s = [&](auto a){ return std::to_string(a);};
    static auto ws = [&](auto a){ std::ostringstream s; s << a; return s.str();};

    A +=  ( ws(simu.get_simuT())+tab);                           //MainCount
    A +=  ( s(chron.p.elapsedTime())+tab);             //SolvePossionEquation
    A +=  ( s(simu.iters)+tab+tab);             //SolvePossionEquation
    A +=  n;

    return A;
}


auto recorderTime_unit(
    std::string N,
    simpulationVariable& simu,
    clockstruct& chron
){

    std::ofstream Time_file;
    std::string fileName = "Information/Chronograph/" + N;
    fileName += std::to_string(simu.TID);
    fileName += ".csv";
    Time_file.open (fileName, std::ios::out|ios::app);
    Time_file << get_Chronograph(simu, chron);
    Time_file.close();

    return true;
}

auto recorderTime_unit_DAT(
    std::string N,
    simpulationVariable& simu,
    clockstruct& chron
){

    std::ofstream Time_file;

    std::string fileName = "Information/Chronograph/" + N;

    fileName += std::to_string(simu.TID);

    fileName += ".dat";

    Time_file.open (fileName, std::ios::out|ios::app);

    Time_file << get_Chronograph_DAT(simu, chron);

    Time_file.close();

    return true;
}




void recorderTime(
    clockstruct& sec_chronograph,
    simpulationVariable& simu,
    shareMenory& ShareM
)
{
    // recorderTime_unit("Information/ATime", sec_chronograph , simu);

    recorderTime_unit("ATime",simu, sec_chronograph);
    recorderTime_unit_DAT("ATime",simu, sec_chronograph);

    if ((simu.loop)%100 == 1){
        recorderTime_unit("BTime",simu, sec_chronograph);
        recorderTime_unit_DAT("BTime",simu, sec_chronograph);
    }

    if ((simu.loop)%1000 == 1){
        recorderTime_unit("CTime",simu, sec_chronograph);
        recorderTime_unit_DAT("CTime",simu, sec_chronograph);
    }

    if ((simu.loop)%10000 == 1){
        recorderTime_unit("DTime",simu, sec_chronograph);
        recorderTime_unit_DAT("DTime",simu, sec_chronograph);
    }

}
