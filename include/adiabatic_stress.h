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

#ifndef ADIABATIC_STRESS_H
#define ADIABATIC_STRESS_H

#include <vector>
#include <memory>

struct AdiabaticStressEntry
{
    double rate;
    double threshold;
    double baseStress;
    double startTime;
    bool triggered ;
    double triggerTime;
    double strainIncrement;
    bool unreadAvalancheData;
    double endTriggerTime;
    double deformation;
    double springConstant;
};

class AdiabaticStress
{
public:
    AdiabaticStress();

    unsigned int createNew(double rate, double threshold, double baseStress, double startTime, double currentDeformation, double D);
    unsigned int createNewBasedOn(const unsigned int & id);
    void deleteCase(const unsigned int & id);
    void deleteAll();

    void updateEntry(const unsigned int & id, const double & t, const double & avgabssp);
    double getStress(const unsigned int & id,const double & t);

    bool triggered(const unsigned int & id);

    double getStrainIncrement(const unsigned int & id);
    void addStrainIncrement(const unsigned int & id, const double & strainIncrement);

    double getTriggerTime(const unsigned int & id);

    double getBaseStress(const unsigned int & id);

    bool avalanche(const unsigned int & id);

    double getEndTriggerTime(const unsigned int & id);
private:
    std::vector<AdiabaticStressEntry> entries;
};

#endif // ADIABATIC_STRESS_H
