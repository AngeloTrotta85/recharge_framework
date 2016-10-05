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

#ifndef INET_POWER_SIMPLEBATTERY_SIMPLEBATTERY_H_
#define INET_POWER_SIMPLEBATTERY_SIMPLEBATTERY_H_


#include "inet/common/INETDefs.h"
#include "inet/common/Units.h"

namespace inet {
namespace power {

class INET_API SimpleBattery : public cSimpleModule {

public:
    typedef enum {
        CHARGING,
        DISCHARGING,
        UNDEFINED_STATE
    } batteryState;

    friend std::ostream& operator<<( std::ostream& os, const batteryState bs )
    {
        if (bs == CHARGING) {
            os << "CHARGING";
        }
        else if (bs == DISCHARGING) {
            os << "DISCHARGING";
        }
        else {
            os << "UNDEFINED_STATE";
        }

        return os;
    }

protected:
  virtual int numInitStages() const override { return NUM_INIT_STAGES; }
  virtual void initialize(int stage) override;
  virtual void handleMessage(cMessage *msg) override;

  virtual void updateBatteryLevel(void);

public:
    SimpleBattery(){}
    virtual ~SimpleBattery() {}

    double getBatteryLevelAbs(void) {updateBatteryLevel(); return batteryLevel;}
    double getBatteryLevelPerc(void) {updateBatteryLevel(); return batteryLevel/fullCapacity;}

    bool isFull(void) {return (batteryLevel >= fullCapacity);}

    void setState(batteryState bs);
    batteryState getState(void) const {return bState;}

    bool isCharging(void) {return bState == CHARGING;}

    double getChargingFactorVal() const {
        return chargingFactor;
    }

    double getDischargingFactorVal() const {
        return dischargingFactor;
    }

    double getChargingFactor(double time) const {
        return (chargingFactor * time);
    }

    double getDischargingFactor(double time) const {
        return (dischargingFactor * time);
    }

    double getSwapLoose(void);

    double getFullCapacity() const {
        return fullCapacity;
    }

private:

    cMessage *autoMsg = nullptr;
    batteryState bState = DISCHARGING;

    simtime_t lastBatteryCheck;

    //parameters
    double batteryLevel;
    double initialCapacity;
    double fullCapacity;
    double updateInterval;

    double chargingFactor;
    double dischargingFactor;
    double flightHeight;
    double swapHeightFactor;

};

} /* namespace power */
} /* namespace inet */

#endif /* INET_POWER_SIMPLEBATTERY_SIMPLEBATTERY_H_ */
