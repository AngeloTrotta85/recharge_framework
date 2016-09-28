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

#ifndef INET_APPLICATIONS_UDPAPP_UDPBASICRECHARGE_H_
#define INET_APPLICATIONS_UDPAPP_UDPBASICRECHARGE_H_

#include "inet/common/INETDefs.h"

#include <simplebattery/SimpleBattery.h>
#include <vector>
#include <map>
#include <list>
#include <iomanip>      // std::setprecision

#include "inet/applications/udpapp/UDPBasicApp.h"

#include "inet/mobility/single/VirtualSpringMobility.h"
#include "inet/applications/base/ApplicationPacketRecharge_m.h"

namespace inet {

class INET_API UDPBasicRecharge : public UDPBasicApp {

public:
    typedef struct {
        L3Address addr;
        int appAddr;
        double rcvPow;
        double rcvSnr;
        Coord pos;
    } nodeInfo_t;

    typedef struct {
        int addr;
        int assignedRecharge;
        int executedRecharge;
    } nodeAlgo_t;

    typedef struct {
        int chargingAppAddr;
        std::list<nodeAlgo_t> nodeList;
    } groupInfo_t;

protected:
  virtual int numInitStages() const override { return NUM_INIT_STAGES; }
  virtual void initialize(int stage) override;
  virtual void handleMessageWhenUp(cMessage *msg) override;

  virtual void sendPacket();
  //virtual void processStart();
  virtual void processPacket(cPacket *msg);

  virtual double calculateInterDistance(double radious);
  virtual void updateVirtualForces(void);

  virtual double calculateRechargeProb(void);
  virtual void checkRecharge(void);
  virtual void checkCentralizedRecharge(void);
  virtual void initCentralizedRecharge(void);

  virtual void decideRechargeSceduling(void);

  int getNodeWithMaxEnergy(groupInfo_t *gi, double &battVal);
  int getNodeWithMinEnergy(groupInfo_t *gi, double &battVal);

public:
    UDPBasicRecharge() {}
    virtual ~UDPBasicRecharge();


private:
    L3Address myAddr;
    int myAppAddr;

    cMessage *autoMsgRecharge = nullptr;
    cMessage *autoMsgCentralizedRecharge = nullptr;

    VirtualSpringMobility *mob = nullptr;
    power::SimpleBattery *sb = nullptr;

    std::map<int, nodeInfo_t> neigh;
    std::list<groupInfo_t> groupList;

    //parameters
    double checkRechargeTimer;
    double sensorRadious;
    Coord rebornPos;

    bool isCentralized;
    int chargingStationNumber;
};

} /* namespace inet */

#endif /* INET_APPLICATIONS_UDPAPP_UDPBASICRECHARGE_H_ */
