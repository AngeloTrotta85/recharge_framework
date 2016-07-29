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

void UDPBasicRecharge::initialize(int stage)
{
    UDPBasicApp::initialize(stage);

    if (stage == INITSTAGE_LOCAL) {

        mob = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getSubmodule("mobility"));
        sb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getSubmodule("battery"));

        double rx = (mob->getConstraintAreaMax().x - mob->getConstraintAreaMin().x) / 2.0;
        double ry = (mob->getConstraintAreaMax().y - mob->getConstraintAreaMin().y) / 2.0;
        rebornPos = Coord(rx, ry);

        EV << "Reborn pos: " << rebornPos << endl;

        checkRechargeTimer = par("checkRechargeTimer");
        sensorRadious = par("sensorRadious");

        autoMsgRecharge = new cMessage("msgRecharge");
        scheduleAt(simTime() + checkRechargeTimer + (dblrand() - 0.5), autoMsgRecharge);
    }
    else if (stage == INITSTAGE_LAST) {
        myAddr = L3AddressResolver().resolve(this->getParentModule()->getFullPath().c_str());
        myAppAddr = this->getParentModule()->getIndex();
        EV << "[" << myAppAddr << "] My address is: " << myAddr << std::endl;

        this->getParentModule()->getDisplayString().setTagArg("t", 0, myAddr.str().c_str());
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

} /* namespace inet */
