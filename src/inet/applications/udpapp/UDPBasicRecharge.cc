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
    }
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
    EV << "RECEIVED PACKET: " << pk->getName() << endl;
    if (sb->getState() == power::SimpleBattery::DISCHARGING) {
        ApplicationPacketRecharge *aPkt = check_and_cast<ApplicationPacketRecharge *> (pk);
        if (myAddr != aPkt->getAddr()) {

            cObject *c = pk->getControlInfo();
            UDPDataIndicationExt *di = check_and_cast<UDPDataIndicationExt *>(c);

            EV_DEBUG << "Received recharge packet " << aPkt->getName() << " with " << di->getFullName() << endl;

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
        EV_INFO << "Received packet: " << UDPSocket::getReceivedPacketInfo(pk) << endl;
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

        EV << "Setting force with displacement: " << springDispl << " (distance: " << distance << ")" << endl;
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

        groupList.front().nodeList.push_front(newNodeInfo);

        actG++;
        if (actG >= numG) {
            actG = 0;
        }
    }

    decideRechargeSceduling();

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

void UDPBasicRecharge::decideRechargeSceduling(void) {
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);
        double maxE;

        int maxNode = getNodeWithMaxEnergy(actGI, maxE);
    }
}

void UDPBasicRecharge::checkCentralizedRecharge(void) {
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        //groupInfo_t *actGI = &(*it);
    }
}

} /* namespace inet */
