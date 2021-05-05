/*
Copyright (C) 2018-2021 Gábor Péterffy <peterffy95@gmail.com>

This file is part of sdddstEQ.

sdddstEQ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sdddstEQ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sdddstEQ. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef PRECISION_HANDLER_H
#define PRECISION_HANDLER_H

#include <utility>
#include <vector>

namespace sdddst {

class PrecisionHandler
{
public:
    PrecisionHandler();
    ~PrecisionHandler();

    void setSize(unsigned int size);
    unsigned int getSize();

    void reset();
    void updateTolerance(const double & distanceSqr, const unsigned int &ID);

    void updateError(const double & error, const unsigned int &ID);

    double getNewStepSize(const double & oldStepSize);

    double getMinPrecisity() const;
    void setMinPrecisity(double value);

    double getMaxErrorRatioSqr() const;

private:
    std::vector<std::pair<double, double> > toleranceAndError;
    double minPrecisity;
    double minPrecisitySqr;
    double maxErrorRatioSqr;
    unsigned int selectedID;
};

}

#endif
