//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
// 

#include <simplebattery/SimpleBattery.h>

namespace inet {
namespace power {

Define_Module(inet::power::SimpleBattery);

void SimpleBattery::initialize(int stage) {
    if (stage == INITSTAGE_LOCAL) {
        initialCapacity = par("initialCapacity");
        fullCapacity = par("nominalCapacity");
        batteryLevel = initialCapacity;
        updateInterval = par("updateInterval");

        autoMsg = new cMessage("batteryLevelUpdate");
        scheduleAt(simTime() + updateInterval, autoMsg);

        WATCH(batteryLevel);
        WATCH(bState);
    }

}

void SimpleBattery::handleMessage(cMessage *msg) {
    if (msg->isSelfMessage()) {
        if (msg == autoMsg){
            updateBatteryLevel();

            scheduleAt(simTime() + updateInterval, autoMsg);
        }
    }
}

void SimpleBattery::updateBatteryLevel(void) {

    if (bState == DISCHARGING) {
        batteryLevel -= 1.0 * updateInterval;

        if(batteryLevel < 0) {
            batteryLevel = 0;
        }
    }
    else if (bState == CHARGING){
        batteryLevel += 5.0 * updateInterval;

        if(batteryLevel > fullCapacity) {
            batteryLevel = fullCapacity;
        }
    }
}

} /* namespace power */
} /* namespace inet */
