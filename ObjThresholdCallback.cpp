#include "ObjThresholdCallback.h"

#include <iostream>

using namespace std;

/*
runs at least 50 seconds.
After that it stops when at least one biclique with objective > 1 was found.
*/
void ObjThresholdCallback::callback() {
    //called on every new found solution
    if (where == GRB_CB_MIPSOL) {
        //starting time set
        if (!timeSet)startTime = getDoubleInfo(GRB_CB_RUNTIME);
        timeSet = true;
        double currentTime = getDoubleInfo(GRB_CB_RUNTIME);
        double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
        //increase count of good solutions
        if (obj > 1.000001) {
            cnt++;
            greater1 = true;
        }
        //check whether we can terminate
        if (cnt>0 && currentTime-startTime>50) {
            GRBCallback::abort();  // ends optimize()
        }
    }
    //called on progress (e.g. also reduced upper bound and similar)
    if (where == GRB_CB_MIP) {
        //set timer
        if (!timeSet)startTime = getDoubleInfo(GRB_CB_RUNTIME);
        timeSet = true;
        double currentTime = getDoubleInfo(GRB_CB_RUNTIME);
        double bestBound = getDoubleInfo(GRB_CB_MIP_OBJBND);
        double obj = getDoubleInfo(GRB_CB_MIP_OBJBST);
        // if bestBound is not good enough we can stop.
        if (bestBound <= 1) {
            cout << "MIP OBJBND = Bestbound = " << bestBound << endl;
            GRBCallback::abort(); 
        }
        //if we found a valid solution bigger than 1 we are happy
        else if (obj > 1.000001 && currentTime - startTime > 50) {
            cout << "MIP OBJBST = Objective = " << obj << endl;
            GRBCallback::abort();
        }
    }
}
