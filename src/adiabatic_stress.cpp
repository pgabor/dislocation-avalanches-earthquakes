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

#include "adiabatic_stress.h"

#include <fstream>
#include <iomanip>
#include <cstring>

AdiabaticStress::AdiabaticStress()
{
    //Nothing to do
}

unsigned int AdiabaticStress::createNew(double rate, double threshold, double baseStress, double startTime, double currentDeformation, double D)
{
    entries.push_back(AdiabaticStressEntry{rate, threshold, baseStress, startTime, false, 0, 0, false, 0, currentDeformation, D});
    return entries.size()-1;
}

unsigned int AdiabaticStress::createNewBasedOn(const unsigned int &id)
{
    entries.push_back(AdiabaticStressEntry{entries[id].rate, entries[id].threshold, entries[id].baseStress, entries[id].startTime,
                                           entries[id].triggered, entries[id].triggerTime, entries[id].strainIncrement,
                                           entries[id].unreadAvalancheData, entries[id].endTriggerTime, entries[id].deformation, entries[id].springConstant});
    return entries.size()-1;
}

void AdiabaticStress::deleteCase(const unsigned int &id)
{
    entries.erase(entries.begin() + id);
}

void AdiabaticStress::deleteAll()
{
    entries.resize(0);
}

void AdiabaticStress::updateEntry(const unsigned int &id, const double &t, const double &avgabssp)
{
    entries[id].deformation += entries[id].strainIncrement;
    entries[id].strainIncrement = 0;
    if (entries[id].threshold < avgabssp)
    {
        if (!entries[id].triggered)
        {
            entries[id].baseStress = 0;
            entries[id].triggered = true;
            entries[id].triggerTime = t;
        }
        entries[id].startTime = t;
    }
    else
    {
        if (entries[id].triggered)
        {
            entries[id].triggered = false;
            entries[id].endTriggerTime = t;
            entries[id].unreadAvalancheData = true;
        }
    }
}

double AdiabaticStress::getStress(const unsigned int &id, const double &t)
{
    return entries[id].baseStress + ((t) * entries[id].rate - entries[id].deformation)*entries[id].springConstant;
}

bool AdiabaticStress::triggered(const unsigned int &id)
{
    return entries[id].triggered;
}

double AdiabaticStress::getStrainIncrement(const unsigned int &id)
{
    return entries[id].strainIncrement;
}

void AdiabaticStress::addStrainIncrement(const unsigned int &id, const double &strainIncrement)
{
    entries[id].strainIncrement += strainIncrement;
}

double AdiabaticStress::getTriggerTime(const unsigned int &id)
{
    return entries[id].triggerTime;
}

double AdiabaticStress::getBaseStress(const unsigned int &id)
{
    return entries[id].baseStress;
}

bool AdiabaticStress::avalanche(const unsigned int &id)
{
    if (entries[id].unreadAvalancheData)
    {
        entries[id].unreadAvalancheData = false;
        return true;
    }
    return false;
}

double AdiabaticStress::getEndTriggerTime(const unsigned int &id)
{
    return entries[id].endTriggerTime;
}
