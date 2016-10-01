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

#include <UDPBasicRecharge.h>

#include <stdio.h>
#include <stdlib.h>

#include <iomanip>      // std::setprecision

#include "inet/networklayer/common/L3AddressResolver.h"

#include "inet/applications/base/ApplicationPacket_m.h"
#include "inet/transportlayer/contract/udp/UDPDataIndicationExt_m.h"

namespace inet {

Define_Module(UDPBasicRecharge)

UDPBasicRecharge::~UDPBasicRecharge() {
    cancelAndDelete(autoMsgRecharge);
    cancelAndDelete(autoMsgCentralizedRecharge);
}

void UDPBasicRecharge::initialize(int stage)
{
    UDPBasicApp::initialize(stage);

    if (stage == INITSTAGE_LOCAL) {

        mob = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getSubmodule("mobility"));
        sb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getSubmodule("battery"));

        myAppAddr = this->getParentModule()->getIndex();

        double rx = (mob->getConstraintAreaMax().x - mob->getConstraintAreaMin().x) / 2.0;
        double ry = (mob->getConstraintAreaMax().y - mob->getConstraintAreaMin().y) / 2.0;
        rebornPos = Coord(rx, ry);

        EV << "Reborn pos: " << rebornPos << endl;

        checkRechargeTimer = par("checkRechargeTimer");
        sensorRadious = par("sensorRadious");
        isCentralized = par("isCentralized").boolValue();
        chargingStationNumber = par("chargingStationNumber");

        autoMsgRecharge = new cMessage("msgRecharge");
        if (isCentralized) {
            if (myAppAddr == 0) {
                autoMsgCentralizedRecharge = new cMessage("msgCentralizedRecharge");
                scheduleAt(simTime() + checkRechargeTimer, autoMsgCentralizedRecharge);
            }
        }
        else {
            scheduleAt(simTime() + checkRechargeTimer + (dblrand() - 0.5), autoMsgRecharge);

        }
    }
    else if (stage == INITSTAGE_LAST) {
        myAddr = L3AddressResolver().resolve(this->getParentModule()->getFullPath().c_str());
        //myAppAddr = this->getParentModule()->getIndex();
        EV << "[" << myAppAddr << "] My address is: " << myAddr << std::endl;

        this->getParentModule()->getDisplayString().setTagArg("t", 0, myAddr.str().c_str());

        if (isCentralized && (myAppAddr == 0)) {
            initCentralizedRecharge();
        }

        EV
        << "BATTERY STEP CHARGE: " << sb->getChargingFactor(checkRechargeTimer)
        << " STEP DISCHARGE: " << sb->getDischargingFactor(checkRechargeTimer)
        << " SWAP LOOSE: " << sb->getSwapLoose()
        << endl;
    }
}

void UDPBasicRecharge::finish(void) {
    recordScalar("LIFETIME", simTime());
}


void UDPBasicRecharge::handleMessageWhenUp(cMessage *msg)
{
    if ((msg->isSelfMessage()) && (msg == autoMsgRecharge)) {
        if ((sb->getState() == power::SimpleBattery::CHARGING) && (sb->isFull())){
            sb->setState(power::SimpleBattery::DISCHARGING);

            mob->clearVirtualSpringsAndsetPosition(rebornPos);
        }

        if (sb->getState() == power::SimpleBattery::DISCHARGING){
            checkRecharge();
            scheduleAt(simTime() + checkRechargeTimer + (dblrand() - 0.5), msg);
        }
        else {
            scheduleAt(simTime() + 0.5, msg);
        }
    }
    else if ((msg->isSelfMessage()) && (msg == autoMsgCentralizedRecharge)) {
        checkCentralizedRecharge();
        scheduleAt(simTime() + checkRechargeTimer, msg);
    }
    else {
        UDPBasicApp::handleMessageWhenUp(msg);
    }
}
/*
void UDPBasicRecharge::processStart()
{
    socket.setOutputGate(gate("udpOut"));
    const char *localAddress = par("localAddress");
    EV << "Local ADDR: '" << localAddress << "'" << endl;
    socket.bind(*localAddress ? L3AddressResolver().resolve(localAddress) : L3Address(), localPort);
    //socket.bind(localPort);
    setSocketOptions();

    const char *destAddrs = par("destAddresses");
    cStringTokenizer tokenizer(destAddrs);
    const char *token;

    while ((token = tokenizer.nextToken()) != nullptr) {
        L3Address result;
        L3AddressResolver().tryResolve(token, result);
        if (result.isUnspecified())
            EV_ERROR << "cannot resolve destination address: " << token << endl;
        else
            destAddresses.push_back(result);
    }

    if (!destAddresses.empty()) {
        selfMsg->setKind(SEND);
        processSend();
    }
    else {
        if (stopTime >= SIMTIME_ZERO) {
            selfMsg->setKind(STOP);
            scheduleAt(stopTime, selfMsg);
        }
    }
}
*/
void UDPBasicRecharge::processPacket(cPacket *pk)
{
    //EV << "RECEIVED PACKET: " << pk->getName() << endl;
    if (sb->getState() == power::SimpleBattery::DISCHARGING) {
        ApplicationPacketRecharge *aPkt = check_and_cast<ApplicationPacketRecharge *> (pk);
        if (myAddr != aPkt->getAddr()) {

            cObject *c = pk->getControlInfo();
            UDPDataIndicationExt *di = check_and_cast<UDPDataIndicationExt *>(c);

            //EV_DEBUG << "Received recharge packet " << aPkt->getName() << " with " << di->getFullName() << endl;

            if (neigh.count(aPkt->getAppAddr()) == 0) {
                nodeInfo_t newInfo;
                newInfo.addr = aPkt->getAddr();
                newInfo.appAddr = aPkt->getAppAddr();

                neigh[aPkt->getAppAddr()] = newInfo;
            }

            nodeInfo_t *node = &(neigh[aPkt->getAppAddr()]);
            node->pos = aPkt->getPos();
            node->rcvPow = di->getPow();
            node->rcvSnr = di->getSnr();

            updateVirtualForces();

        }

        emit(rcvdPkSignal, pk);
        //EV_INFO << "Received packet: " << UDPSocket::getReceivedPacketInfo(pk) << endl;
        numReceived++;
    }

    delete pk;
}


void UDPBasicRecharge::sendPacket()
{
    if (sb->getState() == power::SimpleBattery::DISCHARGING) {
        std::ostringstream str;
        str << packetName << "-" << numSent;
        ApplicationPacketRecharge *payload = new ApplicationPacketRecharge(str.str().c_str());
        payload->setByteLength(par("messageLength").longValue());
        payload->setSequenceNumber(numSent);

        payload->setPos(mob->getCurrentPosition());
        payload->setAddr(myAddr);
        payload->setAppAddr(myAppAddr);
        payload->setBatteryLevelAbs(sb->getBatteryLevelAbs());
        payload->setBatteryLevelAbs(sb->getBatteryLevelPerc());

        L3Address destAddr = chooseDestAddr();
        emit(sentPkSignal, payload);
        socket.sendTo(payload, destAddr, destPort);
        numSent++;
    }
}


double UDPBasicRecharge::calculateInterDistance(double radious) {
    return (sqrt(3.0) * radious);
}

void UDPBasicRecharge::updateVirtualForces(void) {
    Coord myPos = mob->getCurrentPosition();

    // clear everything
    mob->clearVirtualSprings();

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);
        double preferredDistance = calculateInterDistance(sensorRadious);

        double distance = act->pos.distance(myPos);
        double springDispl = preferredDistance - distance;

        Coord uVec = Coord(1, 1);
        if (distance == 0)  uVec = Coord(dblrand(), dblrand());
        else  uVec = act->pos - myPos;
        uVec.normalize();

        //EV << "Setting force with displacement: " << springDispl << " (distance: " << distance << ")" << endl;
        mob->addVirtualSpring(uVec, preferredDistance, springDispl);
    }

    // add the force towards the center rebornPos
    if (rebornPos.distance(myPos) > 10) {
        Coord uVec = rebornPos - myPos;
        //Coord uVec = myPos - rebornPos;
        uVec.normalize();
        mob->addVirtualSpring(uVec, rebornPos.distance(myPos), -3);
    }
}

double UDPBasicRecharge::calculateRechargeProb(void) {
    double ris = 0.1;

    return ris;
}

void UDPBasicRecharge::checkRecharge(void) {
    double prob = calculateRechargeProb();

    if (dblrand() < prob) {
        sb->setState(power::SimpleBattery::CHARGING);

        //stop the node
        neigh.clear();
        //mob->clearVirtualSprings();
        mob->clearVirtualSpringsAndsetPosition(rebornPos);

        //mob->addVirtualSpring(Coord(1,1,0), 1, 1000);
    }
}

void UDPBasicRecharge::initCentralizedRecharge(void) {
    // create the groups
    int numberNodes = this->getParentModule()->getVectorSize();
    int numG = numberNodes / chargingStationNumber;

    if ((numberNodes % chargingStationNumber) > 0) numG++;

    int actG = 0;
    for (int i = 0; i < numberNodes; i++) {


        if (actG == 0) {
            groupInfo_t newG;

            newG.chargingAppAddr = -1;
            newG.swapNumber = 0;
            //newG.nodeList.push_front(newNodeInfo);

            groupList.push_front(newG);
        }
        //else {
        //    groupList.front().nodeList.push_front(i);
        //}

        nodeAlgo_t newNodeInfo;
        newNodeInfo.addr = i;
        newNodeInfo.executedRecharge = 0;
        newNodeInfo.assignedRecharge = 0;
        newNodeInfo.isCharging = false;
        newNodeInfo.energy = 0;

        groupList.front().nodeList.push_front(newNodeInfo);

        actG++;
        if (actG >= numG) {
            actG = 0;
        }
    }



    /*
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        int minChargeAddr = -1;
        double minChargeVal = 0;

        EV << "GRUPPO: ";
        for (auto it2 = actGI->nodeList.begin(); it2 != actGI->nodeList.end(); it2++) {
            int actAddress = it2->addr;
            power::SimpleBattery *battNeigh = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", actAddress)->getSubmodule("battery"));
            double actPow = battNeigh->getBatteryLevelAbs();

            EV << actAddress << "|" << actPow << " ";

            if ((minChargeAddr < 0) || (actPow < minChargeVal)) {
                minChargeAddr = actAddress;
                minChargeVal = actPow;
            }
        }

        actGI->chargingAppAddr = minChargeAddr;
        for (auto it2 = actGI->nodeList.begin(); it2 != actGI->nodeList.end(); it2++) {
            if (it2->addr == minChargeAddr) {
                it2->isCharging = true;
                break;
            }
        }

        EV << " --- : minChargeAddr: " << minChargeAddr << endl;

        // setting the min node in charging
        VirtualSpringMobility *mobMin = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getParentModule()->getSubmodule("host", minChargeAddr)->getSubmodule("mobility"));
        power::SimpleBattery *battMin = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", minChargeAddr)->getSubmodule("battery"));
        UDPBasicRecharge *nodeMin = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", minChargeAddr)->getSubmodule("udpApp", 0));
        battMin->setState(power::SimpleBattery::CHARGING);
        mobMin->clearVirtualSpringsAndsetPosition(rebornPos);
        nodeMin->neigh.clear();

    }
    */
    decideRechargeSceduling();

    checkCentralizedRecharge();

}

int UDPBasicRecharge::getNodeWithMaxEnergy(groupInfo_t *gi, double &battVal) {
    int ris = -1;
    double maxE = -1;

    for (auto it = gi->nodeList.begin(); it != gi->nodeList.end(); it++) {
        nodeAlgo_t *actN = &(*it);
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", actN->addr)->getSubmodule("battery"));

        if (battN->getBatteryLevelAbs() > maxE) {
            ris = actN->addr;
            maxE = battN->getBatteryLevelAbs();
        }
    }

    battVal = maxE;
    return ris;
}

int UDPBasicRecharge::getNodeWithMinEnergy(groupInfo_t *gi, double &battVal) {
    int ris = -1;
    double minE = std::numeric_limits<double>::max();

    for (auto it = gi->nodeList.begin(); it != gi->nodeList.end(); it++) {
        nodeAlgo_t *actN = &(*it);
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", actN->addr)->getSubmodule("battery"));

        if (battN->getBatteryLevelAbs() < minE) {
            ris = actN->addr;
            minE = battN->getBatteryLevelAbs();
        }
    }

    battVal = minE;
    return ris;

}

void UDPBasicRecharge::updateBatteryVals(std::list<nodeAlgo_t> *list) {
    for (auto itn = list->begin(); itn != list->end(); itn++){
        nodeAlgo_t *actNO = &(*itn);
        power::SimpleBattery *batt = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", actNO->addr)->getSubmodule("battery"));
        actNO->energy = batt->getBatteryLevelAbs();
    }
}

// comparison, energy.
//bool UDPBasicRecharge::compare_energy (const nodeAlgo_t& first, const nodeAlgo_t& second) {
bool compare_energy (const inet::UDPBasicRecharge::nodeAlgo_t& first, const inet::UDPBasicRecharge::nodeAlgo_t& second) {
    //power::SimpleBattery *batt1 = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", first.addr)->getSubmodule("battery"));
    //power::SimpleBattery *batt2 = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", second.addr)->getSubmodule("battery"));

    //return (batt1->getBatteryLevelAbs() < batt2->getBatteryLevelAbs());
    //EV << std::scientific << std::setprecision(20) << "First " << first.energy << ":" << first.isCharging << " - Second " << second.energy << ":" << second.isCharging << endl;

    if (fabs(first.energy - second.energy) < EPSILON) {
        if (first.isCharging) {
            return true;
        }
        return false;
    }
    else {
        return first.energy < second.energy;
    }
}

bool compare_charge (const inet::UDPBasicRecharge::nodeAlgo_t& first, const inet::UDPBasicRecharge::nodeAlgo_t& second) {
    if (first.isCharging) {
        return true;
    }
    else if (second.isCharging) {
        return false;
    }
    else {
        return first.energy < second.energy;
    }
}

bool UDPBasicRecharge::decideRechargeSceduling(void) {
    bool ris = true;
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        bool groupRis = decideRechargeScedulingGroup(actGI);

        ris = ris && groupRis;
    }

    return ris;
}

void UDPBasicRecharge::decideRechargeScedulingGroupLast(groupInfo_t *actGI) {
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);
        actNO->assignedRecharge = 0;
    }

    updateBatteryVals(&(actGI->nodeList));
    actGI->nodeList.sort(compare_energy);

    // schedule to make sure "change when dying"
    int sumDischargingSteps = 0;
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);
        auto itn2 = itn;
        itn2++;

        if (itn2 != actGI->nodeList.end()) {
            nodeAlgo_t *nextNO = &(*itn2);

            // to calculate the possible steps check the next node's energy (minus 1), remove all the previous discharging steps
            // and finally divide on the discharge factor

            //EV << "Next energy: " << nextNO->energy << "; steps done: " << sumDischargingSteps;
            actNO->assignedRecharge = (nextNO->energy - 1 - (sumDischargingSteps * sb->getDischargingFactor(checkRechargeTimer))) / sb->getDischargingFactor(checkRechargeTimer);
            //EV << "; assigned: " << actNO->assignedRecharge << endl;
        }
        else {
            actNO->assignedRecharge = (actNO->energy - 1 - (sumDischargingSteps * sb->getDischargingFactor(checkRechargeTimer))) / sb->getDischargingFactor(checkRechargeTimer);
        }

        if (actNO->assignedRecharge < 0) {
            actNO->assignedRecharge = 0;
        }

        sumDischargingSteps += actNO->assignedRecharge;
    }
}

bool UDPBasicRecharge::decideRechargeScedulingGroup(groupInfo_t *actGI) {
    double maxE;
    bool ris = true;

    getNodeWithMaxEnergy(actGI, maxE);

    int numSteps = (maxE - 1) / sb->getDischargingFactor(checkRechargeTimer);
    int numChargeSlots = numSteps / (actGI->nodeList.size() - 1);
    //int plusSteps = ((int) maxE) % ((int)sb->getDischargingFactor());
    int plusSteps = numSteps - (numChargeSlots * (actGI->nodeList.size() - 1));


    updateBatteryVals(&(actGI->nodeList));
    printChargingInfo("BEFORE SORT -> ");
    actGI->nodeList.sort(compare_energy);
    printChargingInfo("AFTER SORT -> ");

    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        actNO->executedRecharge = 0;
        actNO->assignedRecharge = numChargeSlots;
        if (plusSteps > 0) {
            actNO->assignedRecharge++;
            plusSteps--;
        }
    }

    printChargingInfo("BEFORE CHECKING FEASIBILITY -> ");

    //check the feasibility of the schedule
    bool feasible = true;
    int requestSteps = 0;
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);
        requestSteps += actNO->assignedRecharge;
    }

    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        //int possibleSteps = (actNO->energy - 1 + (actNO->assignedRecharge * sb->getChargingFactor())) / sb->getDischargingFactor();
        int possibleSteps = ((actNO->energy - 1) / sb->getDischargingFactor(checkRechargeTimer)) + actNO->assignedRecharge;

        EV << "RequestSteps: " << requestSteps << " - PossibleSteps: " << possibleSteps << endl;

        if ((requestSteps < ((int)actGI->nodeList.size())) || (possibleSteps < requestSteps)) {
            EV << "NOT FEASIBLE!!!" << endl;
            feasible = false;
            break;
        }
    }

    if (!feasible) {
        decideRechargeScedulingGroupLast(actGI);
        ris = false;
    }

    printChargingInfo("AFTER Scheduling decision -> ");

    return ris;
}

void UDPBasicRecharge::printChargingInfo(const char *str) {
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        EV << str;
        for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
            nodeAlgo_t *actNO = &(*itn);

            EV << "[" << actNO->addr << "]" << actNO->energy << "|" << actNO->executedRecharge << "/" << actNO->assignedRecharge << "|" << actNO->isCharging << "  ";
        }
        EV << endl;
    }
}

void UDPBasicRecharge::printChargingInfo(void) {
    printChargingInfo("BATTERY INFO -> ");
}

void UDPBasicRecharge::checkAliveGroup(groupInfo_t *actGI) {
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        if ((actNO->energy < sb->getDischargingFactor(checkRechargeTimer)) && (!actNO->isCharging)) {
            EV << "System dead... finishing the simulation" << endl;
            printChargingInfo("BATTERY FINAL STATE -> ");
            endSimulation();
        }
    }
}

void UDPBasicRecharge::checkCentralizedRecharge(void) {
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        checkAliveGroup(actGI);

        checkCentralizedRechargeGroup(actGI);
    }
    printChargingInfo();
}

void UDPBasicRecharge::checkCentralizedRechargeGroup(groupInfo_t *actGI) {
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        updateBatteryVals(&(actGI->nodeList));

        if (actGI->chargingAppAddr < 0) {
            // no-body is charging (this is the first time)
            actNO->isCharging = true;
            actGI->chargingAppAddr = actNO->addr;

            putNodeInCharging(actNO->addr);

            EV << "FIRST SCHEDULER STEP" << endl;

            // PUT THE OTHERS IN DISCHARGE
            itn++;
            for (; itn != actGI->nodeList.end(); itn++){
                nodeAlgo_t *actNO3 = &(*itn);
                putNodeInDischarging(actNO3->addr);
            }

            break;
        }
        else if (actNO->isCharging) {

            // recharging slot executed
            actNO->executedRecharge++;

            if (actNO->executedRecharge == actNO->assignedRecharge) {
                // recharge fase finished... swap the charging uav with the next one
                //actNO->isCharging = false;

                itn++;
                if (itn != actGI->nodeList.end()) {
                    nodeAlgo_t *actNO2 = &(*itn);

                    if (actNO2->assignedRecharge > 0) {
                        actNO->isCharging = false;
                        putNodeInDischarging(actNO->addr);

                        actNO2->isCharging = true;
                        putNodeInCharging(actNO2->addr);

                        actGI->chargingAppAddr = actNO2->addr;
                        actGI->swapNumber++;
                    }
                }
                else {
                    // I'm the last one, so start again
                    decideRechargeScedulingGroup(actGI);

                    if (!actGI->nodeList.begin()->isCharging) {

                        EV << "THE FIRST ONE IS NOT IN CHARGING !!!!! WHY???" << endl;
                        printChargingInfo("BEFORE SORT -> ");
                        actGI->nodeList.sort(compare_charge);
                        printChargingInfo("AFTER SORT -> ");

                        /*if (nextPossible) {
                            EV << "THE FIRST ONE IS NOT IN CHARGING !!!!! WHY???" << endl;

                            actGI->chargingAppAddr = -1;
                            for (itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
                                nodeAlgo_t *actNO3 = &(*itn);
                                actNO3->isCharging = false;
                                putNodeInDischarging(actNO3->addr);
                            }

                            checkCentralizedRechargeGroup(actGI);
                        }*/
                    }

                    //actGI->chargingAppAddr = -1;
                    //actNO->isCharging = false;
                    //putNodeInDischarging(actNO->addr);

                    //checkCentralizedRechargeGroup(actGI);
                }

                break;
            }
            //else {
                //nothing to do
            //}
        }
    }
}

void UDPBasicRecharge::putNodeInCharging(int addr) {
    // setting the node in charging
    VirtualSpringMobility *mobN = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getParentModule()->getSubmodule("host", addr)->getSubmodule("mobility"));
    power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", addr)->getSubmodule("battery"));
    UDPBasicRecharge *nodeN = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", addr)->getSubmodule("udpApp", 0));

    battN->setState(power::SimpleBattery::CHARGING);
    mobN->clearVirtualSpringsAndsetPosition(rebornPos);
    nodeN->neigh.clear();
}

void UDPBasicRecharge::putNodeInDischarging(int addr) {
    // setting the node in charging
    VirtualSpringMobility *mobN = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getParentModule()->getSubmodule("host", addr)->getSubmodule("mobility"));
    power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", addr)->getSubmodule("battery"));

    battN->setState(power::SimpleBattery::DISCHARGING);
    mobN->clearVirtualSpringsAndsetPosition(rebornPos);
}

} /* namespace inet */
