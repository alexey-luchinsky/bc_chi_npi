/*******************************************************************************
* Copyright 2019 Alexey Luchinsky                                              *
*   email: alexey.luchinsky@gmail.com                                          *
*                                                                              *
* This program is free software: you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License as published by         *
* the Free Software Foundation, either version 3 of the License, or            *
* (at your option) any later version.                                          *
*                                                                              *
* It is distributed in the hope that it will be useful,                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see   <https://www.gnu.org/licenses/>.     *
*******************************************************************************/

#include "segR.h"
#include "EvtGenModels/EvtBcVHad.hh"
#include "EvtGenBase/EvtModel.hh"

using std::endl;
using std::ofstream;
using std::cout;




int main(int argc, char* argv[]) {
    TStopwatch timer;
    timer.Start();

    EvtModel &modellist = EvtModel::instance();
    modellist.registerModel(new EvtBcVHad());
    
    segR(argc, argv);
    
    //out.close();
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);

    return 1;
}
