#pragma once

#include "gurobi_c++.h"

/*
class for Gurobi interruptions
*/
class ObjThresholdCallback : public GRBCallback {
private:
    int cnt = 0;
    bool greater1 = false;
    bool timeSet = false;
    double startTime;

public:
    void callback() override;
};
