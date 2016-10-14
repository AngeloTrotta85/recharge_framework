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

#include "inet/common/geometry/common/Coord.h"
namespace inet {

Define_Module(UDPBasicRecharge)

UDPBasicRecharge::~UDPBasicRecharge() {
    cancelAndDelete(autoMsgRecharge);
    cancelAndDelete(autoMsgCentralizedRecharge);
    cancelAndDelete(stat1sec);
    cancelAndDelete(stat5sec);
    cancelAndDelete(dischargeTimer);
    cancelAndDelete(goToCharge);
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

        lastRechargeTimestamp = simTime();

        EV << "Reborn pos: " << rebornPos << endl;

        checkRechargeTimer = par("checkRechargeTimer");
        sensorRadious = par("sensorRadious");
        //isCentralized = par("isCentralized").boolValue();
        chargingStationNumber = par("chargingStationNumber");
        stimulusExponent = par("stimulusExponent");
        roundrobinRechargeSize = par("roundrobinRechargeSize");
        numRechargeSlotsProbabilistic = par("numRechargeSlotsProbabilistic");
        makeLowEnergyFactorCurves = par("makeLowEnergyFactorCurves").boolValue();
        timeFactorMultiplier = par("timeFactorMultiplier");
        godCheckIfRechargeStationFree = par("godCheckIfRechargeStationFree").boolValue();
        numRechargeSlotsStimulusZeroNeigh = par("numRechargeSlotsStimulusZeroNeigh");
        returnBackAfterRecharge = par("returnBackAfterRecharge").boolValue();
        stationANDnodeKNOWN = par("stationANDnodeKNOWN").boolValue();

        //logFile = par("analticalLogFile").str();
        printAnalticalLog = par("printAnalticalLog").boolValue();
        snprintf(logFile, sizeof(logFile), "%s", par("analticalLogFile").stringValue());
        remove(logFile);

        firstRecharge = true;
        lastPosBeforeCharge = rebornPos;
        rechargeLostAccess = 0;

        std::string schedulingType = par("schedulingType").stdstringValue();
        //ANALYTICAL, ROUNDROBIN, STIMULUS
        if (schedulingType.compare("ANALYTICAL") == 0) {
            st = ANALYTICAL;
        }
        else if (schedulingType.compare("ROUNDROBIN") == 0) {
            st = ROUNDROBIN;
        }
        else if (schedulingType.compare("STIMULUS") == 0) {
            st = STIMULUS;
        }
        else if (schedulingType.compare("PROBABILISTIC") == 0) {
            st = PROBABILISTIC;
        }
        else {
            error("Wrong \"schedulingType\" parameter");
        }

        isCentralized = false;
        if ((st == ANALYTICAL) || (st == ROUNDROBIN)) {
            isCentralized = true;
        }

        autoMsgRecharge = new cMessage("msgRecharge");
        if (isCentralized) {
            if (myAppAddr == 0) {
                autoMsgCentralizedRecharge = new cMessage("msgCentralizedRecharge");
                scheduleAt(simTime() + checkRechargeTimer, autoMsgCentralizedRecharge);
            }
        }
        else {
            //scheduleAt(simTime() + checkRechargeTimer + (dblrand() - 0.5), autoMsgRecharge);
            //scheduleAt(simTime() + (checkRechargeTimer * dblrand()), autoMsgRecharge);
            scheduleAt(simTime() + checkRechargeTimer, autoMsgRecharge);
            //scheduleAt(simTime() + checkRechargeTimer + ((dblrand() - 0.5) / 2.0), autoMsgRecharge);
        }

        stat1sec = new cMessage("stat1secMsg");
        scheduleAt(simTime(), stat1sec);

        stat5sec = new cMessage("stat5secMsg");
        scheduleAt(simTime(), stat5sec);

        dischargeTimer = new cMessage("dischargeTimer");
        goToCharge = new cMessage("goToCharge");

        //personalUniqueCoverageVector.setName("personalCoverage");
        totalCoverageVector.setName("totalCoverage");
        activeNodesVector.setName("activeNodes");
        rechargingNodesVector.setName("rechargingNodes");
        stimulusVector.setName("StimulusVal");
        thresholdVector.setName("ThresholdVal");
        responseVector.setName("ResponseVal");
        degreeVector.setName("DegreeVal");
        timeFactorVector.setName("TimeFactorVal");
        energyFactorVector.setName("EnergyFactorVal");
        energyVector.setName("EnergyVal");

        WATCH(st);
    }
    else if (stage == INITSTAGE_LAST) {
        myAddr = L3AddressResolver().resolve(this->getParentModule()->getFullPath().c_str());
        //myAppAddr = this->getParentModule()->getIndex();
        EV << "[" << myAppAddr << "] My address is: " << myAddr << std::endl;

        this->getParentModule()->getDisplayString().setTagArg("t", 0, myAddr.str().c_str());

        if (isCentralized && (myAppAddr == 0)) {
            initCentralizedRecharge();
        }

        if(!isCentralized){
            sb->setState(power::SimpleBattery::DISCHARGING);
        }

        EV
        << "BATTERY STEP CHARGE: " << sb->getChargingFactor(checkRechargeTimer)
        << " STEP DISCHARGE: " << sb->getDischargingFactor(checkRechargeTimer)
        << " SWAP LOOSE: " << sb->getSwapLoose()
        << endl;
    }
}

void UDPBasicRecharge::finish(void) {
    if (myAppAddr == 0) {

        if (isCentralized) {

            double sumEnergy = 0.0;
            for (auto it = groupList.begin(); it != groupList.end(); it++) {
                groupInfo_t *actGI = &(*it);
                for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
                    nodeAlgo_t *actNO = &(*itn);

                    sumEnergy += actNO->energy;
                }
            }
            recordScalar("FINALENERGY", sumEnergy);

            int g = 1;
            for (auto it = groupList.begin(); it != groupList.end(); it++) {
                char buff[64];
                groupInfo_t *actGI = &(*it);

                snprintf(buff, sizeof(buff), "GSWAP%d", g++);
                recordScalar(buff, actGI->swapNumber);
            }
        }
        else {
            int numberNodes = this->getParentModule()->getVectorSize();

            double sumEnergy = 0.0;
            for (int i = 0; i < numberNodes; i++) {
                power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));
                sumEnergy += battN->getBatteryLevelAbs();
            }
            recordScalar("FINALENERGY", sumEnergy);
        }

        recordScalar("LIFETIME", simTime());
    }
}


void UDPBasicRecharge::handleMessageWhenUp(cMessage *msg)
{
    if ((msg->isSelfMessage()) && (msg == goToCharge)) {

        firstRecharge = false;

        if (checkRechargingStationFree()) {
            double cTime = calculateRechargeTime(true);

            rechargeLostAccess = 0;

            lastPosBeforeCharge = mob->getCurrentPosition();

            sb->setState(power::SimpleBattery::CHARGING);

            fprintf(stderr, "[%d] - Going in charging: %f\n", myAppAddr, cTime);fflush(stderr);

            scheduleAt(simTime() + cTime, dischargeTimer);

            if (printAnalticalLog) {
                FILE *f = fopen(logFile, "a");
                if (f) {
                    std::stringstream ss;
                    printDistributedChargingInfo(ss, "BATTERY INFO -> ");
                    fwrite(ss.str().c_str(), 1 , ss.str().length(), f);
                    fclose(f);
                }
                else {
                    fprintf(stderr, "Error opening file: %s \n", logFile); fflush(stderr);
                    error("Error writing on file\n");
                }
            }
        }
        else {
            rechargeLostAccess++;
        }

        //stop the node
        neigh.clear();
        //mob->clearVirtualSprings();
        mob->clearVirtualSpringsAndsetPosition(rebornPos);
    }
    else if ((msg->isSelfMessage()) && (msg == dischargeTimer)) {
        sb->setState(power::SimpleBattery::DISCHARGING);
        neigh.clear();
        if (returnBackAfterRecharge) {
            mob->clearVirtualSpringsAndsetPosition(lastPosBeforeCharge);
        }
        else {
            mob->clearVirtualSpringsAndsetPosition(rebornPos);
        }
        lastRechargeTimestamp = simTime();

        if (printAnalticalLog) {
            FILE *f = fopen(logFile, "a");
            if (f) {
                std::stringstream ss;
                printDistributedChargingInfo(ss, "BATTERY INFO -> ");
                fwrite(ss.str().c_str(), 1 , ss.str().length(), f);
                fclose(f);
            }
            else {
                fprintf(stderr, "Error opening file: %s \n", logFile); fflush(stderr);
                error("Error writing on file\n");
            }
        }
    }
    else if ((msg->isSelfMessage()) && (msg == stat5sec)) {
        make5secStats();
        scheduleAt(simTime() + 5, msg);
    }
    else if ((msg->isSelfMessage()) && (msg == stat1sec)) {
        updateNeighbourhood();
        if (myAppAddr == 0) {
            checkAliveDistributed();
        }

        make1secStats();
        scheduleAt(simTime() + 1, msg);
    }
    else if ((msg->isSelfMessage()) && (msg == autoMsgRecharge)) {
        //if ((sb->getState() == power::SimpleBattery::CHARGING) && (sb->isFull())){
        //    sb->setState(power::SimpleBattery::DISCHARGING);
        //
        //    mob->clearVirtualSpringsAndsetPosition(rebornPos);
        //}

        if (sb->getState() == power::SimpleBattery::CHARGING){
            checkDischarge();
        }
        else if (sb->getState() == power::SimpleBattery::DISCHARGING){
            checkRecharge();
            //scheduleAt(simTime() + checkRechargeTimer + (dblrand() - 0.5), msg);
        }
        //else {
        //    scheduleAt(simTime() + 0.5, msg);
        //}

        //scheduleAt(simTime() + checkRechargeTimer + ((dblrand() - 0.5) / 2.0), msg);
        //scheduleAt(simTime() + (checkRechargeTimer / (((double) rechargeLostAccess) + 1.0)) + ((dblrand() - 0.5) / 2.0), msg);

        scheduleAt(simTime() + checkRechargeTimer, msg);
    }
    else if ((msg->isSelfMessage()) && (msg == autoMsgCentralizedRecharge)) {
        checkCentralizedRecharge();
        scheduleAt(simTime() + checkRechargeTimer, msg);
    }
    else {
        UDPBasicApp::handleMessageWhenUp(msg);
    }
}

void UDPBasicRecharge::make5secStats(void) {
    if (myAppAddr == 0) {
        //totalCoverageVector.record(getFullCoverage());
    }
}

void UDPBasicRecharge::make1secStats(void) {
    //double persCoverage = getMyCoverageActual() / getMyCoverageMax();
    //personalUniqueCoverageVector.record(persCoverage);
    if (myAppAddr == 0) {
        // TODO
        //totalCoverageVector.record(getFullCoverage());

        int nnodesActive = 0;
        int nnodesRecharging = 0;

        int numberNodes = this->getParentModule()->getVectorSize();

        for (int i = 0; i < numberNodes; i++) {
            power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));

            if (battN->isCharging()) {
                nnodesRecharging++;
            }
            else {
                nnodesActive++;
            }
        }

        activeNodesVector.record(nnodesActive);
        rechargingNodesVector.record(nnodesRecharging);
    }
    if (!isCentralized) {
        stimulusVector.record(calculateRechargeStimuli());
        thresholdVector.record(calculateRechargeThreshold());
        responseVector.record(calculateRechargeProb());
        degreeVector.record(calculateNodeDegree());
        timeFactorVector.record(calculateRechargeStimuliTimeFactor());
        energyFactorVector.record(calculateRechargeStimuliEnergyFactor());
    }
    energyVector.record(sb->getBatteryLevelAbs());
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

            if (aPkt->getGoingToRecharge()) {
                if (neigh.count(aPkt->getAppAddr()) != 0) {
                    neigh.erase(aPkt->getAppAddr());
                }
                lastRechargeTimestamp = simTime();
                firstRecharge = false;


            }
            else {

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

                node->timestamp = simTime();

                node->batteryLevelAbs = aPkt->getBatteryLevelAbs();
                node->batteryLevelPerc = aPkt->getBatteryLevelPerc();
                node->coveragePercentage = aPkt->getCoveragePercentage();
                node->leftLifetime = aPkt->getLeftLifetime();
                node->nodeDegree = aPkt->getNodeDegree();
            }

            updateVirtualForces();

            EV << "[" << myAppAddr << "] - NEW PACKET arrived. Neigh size: " << neigh.size() << endl;

        }

        emit(rcvdPkSignal, pk);
        //EV_INFO << "Received packet: " << UDPSocket::getReceivedPacketInfo(pk) << endl;
        numReceived++;
    }

    delete pk;
}


void UDPBasicRecharge::sendRechargeMessage(void) {
    std::ostringstream str;
    str << packetName << "-RECHARGE-" << numSent;

    //fprintf(stderr, "sendRechargeMessage OK 1\n");fflush(stderr);

    ApplicationPacketRecharge *payload = new ApplicationPacketRecharge(str.str().c_str());
    payload->setByteLength(par("messageLength").longValue());
    payload->setSequenceNumber(numSent);

    payload->setPos(mob->getCurrentPosition());
    payload->setAddr(myAddr);
    payload->setAppAddr(myAppAddr);
    payload->setBatteryLevelAbs(sb->getBatteryLevelAbs());
    payload->setBatteryLevelPerc(sb->getBatteryLevelPerc());
    //payload->setCoveragePercentage(getMyCoverageActual() / getMyCoverageMax());
    payload->setCoveragePercentage(0);
    payload->setLeftLifetime(sb->getBatteryLevelAbs() / sb->getDischargingFactor(checkRechargeTimer));
    payload->setNodeDegree(calculateNodeDegree());
    payload->setGoingToRecharge(true);

    L3Address destAddr = chooseDestAddr();
    emit(sentPkSignal, payload);
    socket.sendTo(payload, destAddr, destPort);
    numSent++;
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
        payload->setBatteryLevelPerc(sb->getBatteryLevelPerc());
        //payload->setCoveragePercentage(getMyCoverageActual() / getMyCoverageMax());
        payload->setCoveragePercentage(0);
        payload->setLeftLifetime(sb->getBatteryLevelAbs() / sb->getDischargingFactor(checkRechargeTimer));
        payload->setNodeDegree(calculateNodeDegree());
        payload->setGoingToRecharge(false);

        L3Address destAddr = chooseDestAddr();
        emit(sentPkSignal, payload);
        socket.sendTo(payload, destAddr, destPort);
        numSent++;
    }
}

int UDPBasicRecharge::calculateNodeDegree(void) {
    /*int actDegree = 0;
    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);
        if ((act->pos.distance(mob->getCurrentPosition())) < (2.5 * sensorRadious)) {
            actDegree++;
        }
    }

    return actDegree;*/

    std::map<int, nodeInfo_t> filteredNeigh;
    getFilteredNeigh(filteredNeigh);
    return filteredNeigh.size();
}

void UDPBasicRecharge::updateNeighbourhood(void) {
    bool removed;
    do {
        removed = false;
        for (auto it = neigh.begin(); it != neigh.end(); it++) {
            nodeInfo_t *act = &(it->second);

            if ((simTime() - act->timestamp) > (2.0 * par("sendInterval").doubleValue())) {

                EV << "[" << myAppAddr << "] - UPDATE NEIGHBOURHOOD. Removing: " << it->second.appAddr << endl;

                neigh.erase(it);
                removed = true;
                break;
            }
        }
    } while(removed);
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

bool UDPBasicRecharge::checkRechargingStationFree(void) {
    bool ris = true;
    int numberNodes = this->getParentModule()->getVectorSize();
    int nodesInCharging = 0;

    for (int i = 0; i < numberNodes; i++) {
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));

        if (battN->isCharging()) {
            nodesInCharging++;
        }
    }

    if (nodesInCharging >= chargingStationNumber) {
        ris = false;
    }

    return ris;
}

double UDPBasicRecharge::calculateRechargeProb(void) {

    if (st == STIMULUS) {
        double stim = calculateRechargeStimuli();
        double tetha = calculateRechargeThreshold();

        return (pow(stim, stimulusExponent) / (pow(stim, stimulusExponent) + pow(tetha, stimulusExponent)));
    }
    else {//if (st == PROBABILISTIC) {
        return (1 - sb->getBatteryLevelPerc());
    }

}

double UDPBasicRecharge::calculateRechargeStimuliTimeFactor(void) {
    //double averageE = 0;
    double sumE = sb->getBatteryLevelAbs();
    double maxE = sb->getBatteryLevelAbs();

    if (sb->getState() != power::SimpleBattery::DISCHARGING){
        return 0;
    }

    std::map<int, nodeInfo_t> filteredNeigh;
    getFilteredNeigh(filteredNeigh);

    for (auto it = filteredNeigh.begin(); it != filteredNeigh.end(); it++) {
    //for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);

        sumE += act->batteryLevelAbs;
        if (act->batteryLevelAbs > maxE)
            maxE = act->batteryLevelAbs;
    }
    //averageE = sumE / ((double) (neigh.size() + 1));
    //averageE = sumE / ((double) (filteredNeigh.size() + 1));

    double rechargeEstimation = ((sb->getFullCapacity() - sb->getBatteryLevelAbs()) / sb->getChargingFactor(checkRechargeTimer)) * checkRechargeTimer;
    if (filteredNeigh.size() > 0) {
    //if (neigh.size() > 0) {
        //rechargeEstimation = (averageE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;
        //rechargeEstimation = (maxE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;
        rechargeEstimation = (maxE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) filteredNeigh.size()))) * checkRechargeTimer;
    }

    //double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * timeFactorMultiplier);
    double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * ((double) filteredNeigh.size()));
    //timeFactor = timeFactor / timeFactorMultiplier; //TODO
    if (timeFactor > 1) timeFactor = 1;

    return timeFactor;

}

double UDPBasicRecharge::calculateRechargeStimuliEnergyFactor(void) {
    double maxE = sb->getBatteryLevelAbs();
    double minE = maxE;

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);
        double actBatt = act->batteryLevelAbs;

        if (actBatt > maxE) {
            maxE = actBatt;
        }
        if (actBatt < minE) {
            minE = actBatt;
        }

        //fprintf(stderr, "BATT: %f, MAX: %f, MIN: %f\n", actBatt, maxE, minE);fflush(stderr);
    }
    if (maxE == minE) {
        return 1;
    }
    if (sb->getState() != power::SimpleBattery::DISCHARGING){
        return 1;
    }
    double ris = (sb->getBatteryLevelAbs() - minE) / (maxE - minE);

    //fprintf(stderr, "[%d] - Energy Factor. MAX: %f; MIN: %f, MY: %f, Ris: %f\n", myAppAddr, maxE, minE, sb->getBatteryLevelAbs(), ris);fflush(stderr);

    return ris;

}

double UDPBasicRecharge::calculateRechargeStimuli(void) {
    /*double averageE = 0;
    double maxE = sb->getBatteryLevelAbs();
    double sumE = sb->getBatteryLevelAbs();

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);

        sumE += act->batteryLevelAbs;

        if (act->batteryLevelAbs > maxE) {
            maxE = act->batteryLevelAbs;
        }
    }
    averageE = sumE / ((double) (neigh.size() + 1));

    double rechargeEstimation = ((sb->getFullCapacity() - sb->getBatteryLevelAbs()) / sb->getChargingFactor(checkRechargeTimer)) * checkRechargeTimer;
    if (neigh.size() > 0) {
        rechargeEstimation = averageE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()));
    }
    double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / rechargeEstimation;
    if (timeFactor > 1) timeFactor = 1;
    double energyFactor = sb->getBatteryLevelAbs() / maxE;

    return pow(timeFactor, energyFactor);*/
    if (sb->getState() != power::SimpleBattery::DISCHARGING){
        return 0;
    }
    if (firstRecharge) {
        return dblrand();
    }
    else {
        if (makeLowEnergyFactorCurves) {
            return pow(calculateRechargeStimuliTimeFactor(), (1.0 / (1.0 - calculateRechargeStimuliEnergyFactor())));
        }
        else {
            return pow(calculateRechargeStimuliTimeFactor(), calculateRechargeStimuliEnergyFactor());
        }
    }
}

/*
double UDPBasicRecharge::calculateRechargeStimuli(void) {
    double s = 1.0 - (getMyCoverageActual() / getMyCoverageMax());

    if (s < 0) s = 0;
    else if (s > 1) s = 1;

    return s;
}
*/

double UDPBasicRecharge::calculateRechargeThreshold(void) {
    int myDegree = calculateNodeDegree();
    /*int maxDegree = myDegree;
    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);
        if (maxDegree < act->nodeDegree) {
            maxDegree = act->nodeDegree;
        }
    }
    if (maxDegree > 0) {
        return (((double)myDegree) / ((double)maxDegree));
    }
    else {
        return 0;
    }*/
    double ris = ((double) myDegree) / 6.0;
    if (ris > 1) ris = 1;
    return ris;
}

/*
double UDPBasicRecharge::calculateRechargeThreshold(void) {
    double t = 0;
    double myleft = sb->getBatteryLevelAbs() / sb->getDischargingFactor(checkRechargeTimer);
    double maxLeftTime = myleft;


    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);

        if (act->leftLifetime > maxLeftTime) {
            maxLeftTime = act->leftLifetime;
        }
    }

    t = myleft / maxLeftTime;

    if (t < 0) t = 0;
    else if (t > 1) t = 1;

    return t;

}*/

void UDPBasicRecharge::getFilteredNeigh(std::map<int, nodeInfo_t> &filteredNeigh){
    std::list<VirtualSpringMobility::NodeBasicInfo> nodesToFilter;
    std::list<VirtualSpringMobility::NodeBasicInfo> nodeFiltered;

    filteredNeigh.clear();

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        VirtualSpringMobility::NodeBasicInfo newinfo;
        newinfo.id = it->first;
        newinfo.position = it->second.pos;
        nodesToFilter.push_back(newinfo);
    }
    mob->filterNodeListAcuteAngleTest(nodesToFilter, nodeFiltered);

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);

        for (auto it2 = nodeFiltered.begin(); it2 != nodeFiltered.end(); it2++) {
            if (act->appAddr == it2->id) {
                filteredNeigh[it->first] = it->second;
                break;
            }
        }
    }
}

double UDPBasicRecharge::calculateRechargeTime(bool log) {

    double recTime = 0;
    std::stringstream ss;

    if (st == STIMULUS) {

        // default charge until full charge
        //double tt = ((sb->getFullCapacity() - sb->getBatteryLevelAbs()) / sb->getChargingFactor(checkRechargeTimer)) * checkRechargeTimer;
        double tt = numRechargeSlotsStimulusZeroNeigh * checkRechargeTimer;

        if (log) ss << "RECHARGETIME STIMULUS: Default charge time: " << tt << " - Neigh size: " << neigh.size() << endl;

        if (neigh.size() > 0) {

            std::map<int, nodeInfo_t> filteredNeigh;
            getFilteredNeigh(filteredNeigh);

            double sumE = sb->getBatteryLevelAbs();
            double maxE = sb->getBatteryLevelAbs();
            for (auto it = filteredNeigh.begin(); it != filteredNeigh.end(); it++) {
            //for (auto it = neigh.begin(); it != neigh.end(); it++) {
                nodeInfo_t *act = &(it->second);

                sumE += act->batteryLevelAbs;
                if (act->batteryLevelAbs > maxE)
                    maxE = act->batteryLevelAbs;

                break;
            }

            //double averageE = sumE / (((double) neigh.size()) + 1.0);
            //double averageE = sumE / (((double) nodeFiltered.size()) + 1.0);
            //double averageE = sumE / (((double) filteredNeigh.size()) + 1.0);

            if (log) ss << "RECHARGETIME STIMULUS: Max Energy: " << maxE
            //if (log) ss << "RECHARGETIME STIMULUS: Average Energy: " << averageE << ", Max Energy: " << maxE
                    << " - Discharging Factor: " << sb->getDischargingFactor(checkRechargeTimer)
                    << " - SwapLoose Factor: " << sb->getSwapLoose()
                    << " - NodeFiltered size: " << filteredNeigh.size()
                    << " - Neigh size: " << neigh.size()
                    << " - checkRechargeTimer: " << checkRechargeTimer
                    << endl;

            double numSteps = (maxE - (2.0 * sb->getSwapLoose())) / sb->getDischargingFactor(checkRechargeTimer);
            //int actualNeigh = neigh.size();
            //if (actualNeigh > 7) actualNeigh = 7;
            //double numChargeSlots = numSteps / ((double) neigh.size());
            //double numChargeSlots = numSteps / ((double) actualNeigh);
            //double numChargeSlots = numSteps / ((double) nodeFiltered.size());
            double numChargeSlots;
            if (stationANDnodeKNOWN) {
                int numberNodes = this->getParentModule()->getVectorSize();
                numChargeSlots = numSteps / (((double) numberNodes) / ((double) chargingStationNumber));
            }
            else {
                numChargeSlots = numSteps / ((double) filteredNeigh.size() + 1.0);
            }
            //double numChargeSlots = numSteps / ((double) filteredNeigh.size() + 1.0);
            tt = numChargeSlots * checkRechargeTimer;

            //tt = (averageE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;
            //tt = (maxE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;
        }

        recTime = tt;
    }
    else { //if (st == PROBABILISTIC) {
        recTime = numRechargeSlotsProbabilistic * checkRechargeTimer;
    }

    if (log) {
        ss << "RECHARGETIME Final decision charge time: " << recTime << endl;

        fprintf(stderr, "%s", ss.str().c_str());
    }

    return recTime;
}

double UDPBasicRecharge::calculateSendBackoff(void){
    double cw = 1;
    double ris = 0;

    cw = cw / (((double)rechargeLostAccess) + 1.0);

    ris = cw * dblrand();
    ris += 0.01;

    return ris;

}

void UDPBasicRecharge::checkRecharge(void) {
    double prob = calculateRechargeProb();

    if (dblrand() < prob) {
        double backoff = calculateSendBackoff();
        if (godCheckIfRechargeStationFree) {
            if (checkRechargingStationFree()) {
                sendRechargeMessage();
                scheduleAt(simTime() + backoff, goToCharge);
            }
        }
        else {
            sendRechargeMessage();
            scheduleAt(simTime() + backoff, goToCharge);
        }
    }
    else {
        rechargeLostAccess = 0;
    }


    //if (checkRechargingStationFree() && (dblrand() < prob)) {
        //fprintf(stderr, "checkRecharge OK 1\n");fflush(stderr);
        //sendRechargeMessage();
        //fprintf(stderr, "checkRecharge OK 2\n");fflush(stderr);
        //scheduleAt(simTime() + 0.01, goToCharge);
        //fprintf(stderr, "checkRecharge OK 3\n");fflush(stderr);
        /*
        sb->setState(power::SimpleBattery::CHARGING);

        //stop the node
        neigh.clear();
        //mob->clearVirtualSprings();
        mob->clearVirtualSpringsAndsetPosition(rebornPos);

        scheduleAt(simTime() + calculateRechargeTime(), dischargeTimer);
        */

        //mob->addVirtualSpring(Coord(1,1,0), 1, 1000);

        //godCheckIfRechargeStationFree
    //}
}

double UDPBasicRecharge::calculateDischargeProb(void) {
    if (sb->isFull()) {
        return 1;
    }
    else {
        return 0;
    }
}

void UDPBasicRecharge::checkDischarge(void) {
    double prob = calculateDischargeProb();

    if (dblrand() < prob) {
        sb->setState(power::SimpleBattery::DISCHARGING);

        //stop the node
        neigh.clear();
        //mob->clearVirtualSprings();
        mob->clearVirtualSpringsAndsetPosition(rebornPos);

        lastRechargeTimestamp = simTime();

        if (printAnalticalLog) {
            FILE *f = fopen(logFile, "a");
            if (f) {
                std::stringstream ss;
                printDistributedChargingInfo(ss, "BATTERY INFO -> ");
                fwrite(ss.str().c_str(), 1 , ss.str().length(), f);
                fclose(f);
            }
            else {
                fprintf(stderr, "Error opening file: %s \n", logFile); fflush(stderr);
                error("Error writing on file\n");
            }
        }

        if (dischargeTimer->isScheduled()) {
            cancelEvent(dischargeTimer);
        }
    }

}

void UDPBasicRecharge::checkAliveDistributed(void) {
    int numberNodes = this->getParentModule()->getVectorSize();

    for (int i = 0; i < numberNodes; i++) {
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));

        if (battN->getBatteryLevelAbs() <= 0) {
            endSimulation();
        }
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
        bool groupRis = true;

        if (st == ANALYTICAL) {
            groupRis = decideRechargeScedulingGroup(actGI);
        }
        else if (st == ROUNDROBIN){
            groupRis = decideRechargeScedulingGroupRR(actGI);
        }

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
    actGI->nodeList.sort(compare_charge);

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

        if (actNO->assignedRecharge == 0) break;
    }
}

bool UDPBasicRecharge::decideRechargeScedulingGroupRR(groupInfo_t *actGI) {
    bool ris = true;

    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        actNO->executedRecharge = 0;
        actNO->assignedRecharge = roundrobinRechargeSize;
    }

    return ris;
}

bool UDPBasicRecharge::checkScheduleFeasibilityGroup(groupInfo_t *actGI) {
    //check the feasibility of the schedule
    bool feasible = true;
    int rechargeSteps = 0;

    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);
        rechargeSteps += actNO->assignedRecharge;
    }

    int beforeME = 0;
    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

        //int othersRechargeSteps = rechargeSteps - actNO->assignedRecharge;

        //int possibleSteps = (actNO->energy - 1 + (actNO->assignedRecharge * sb->getChargingFactor())) / sb->getDischargingFactor();
        //int possibleSteps = ((actNO->energy - 1 - (2.0 * sb->getSwapLoose())) / sb->getDischargingFactor(checkRechargeTimer))
        //        + actNO->assignedRecharge;


        // check if I can reach my recharge slots
        if (feasible) {
            double energyBeforeRecharge = actNO->energy -1 - (beforeME * sb->getDischargingFactor(checkRechargeTimer)) - sb->getSwapLoose();
            if (energyBeforeRecharge <= 0) {
                feasible = false;
                EV << "I cannot reach the recharge steps. EnergyBeforeRecharge: " << energyBeforeRecharge << endl;
            }
            else {
                EV << "I can reach the recharge steps. EnergyBeforeRecharge: " << energyBeforeRecharge << endl;
            }
        }

        if (feasible) {
            int nextME = rechargeSteps - beforeME - actNO->assignedRecharge;
            double energyAfterRecharge = actNO->energy - 1
                    - (beforeME * sb->getDischargingFactor(checkRechargeTimer))
                    - (2 * sb->getSwapLoose())
                    + (actNO->assignedRecharge * sb->getChargingFactor(checkRechargeTimer));
            if ((energyAfterRecharge - (nextME * sb->getDischargingFactor(checkRechargeTimer))) <= 0 ) {
                feasible = false;
                EV << "I cannot reach the final steps. EnergyAtTheEnd: " << (energyAfterRecharge - (nextME * sb->getDischargingFactor(checkRechargeTimer))) << endl;
            }
            else {
                EV << "I can reach the final steps. EnergyAtTheEnd: " << (energyAfterRecharge - (nextME * sb->getDischargingFactor(checkRechargeTimer))) << endl;
            }
        }

        //int requestSteps = rechargeSteps;

        //EV << "RequestSteps: " << rechargeSteps << " - PossibleSteps: " << possibleSteps << endl;

        //if ((rechargeSteps < ((int)actGI->nodeList.size())) || (possibleSteps < rechargeSteps)) {
        //    EV << "NOT FEASIBLE!!!" << endl;
        //    feasible = false;
        //   break;
        //}

        beforeME += actNO->assignedRecharge;
    }

    return feasible;
}

bool UDPBasicRecharge::decideRechargeScedulingGroup(groupInfo_t *actGI) {
    double maxE;
    bool ris = true;

    getNodeWithMaxEnergy(actGI, maxE);

    //int numSteps = (maxE - 1) / sb->getDischargingFactor(checkRechargeTimer);
    int numSteps = (maxE - (2.0 * sb->getSwapLoose())) / sb->getDischargingFactor(checkRechargeTimer);
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

    if (!checkScheduleFeasibilityGroup(actGI)) {

        for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
            nodeAlgo_t *actNO = &(*itn);

            actNO->executedRecharge = 0;
            actNO->assignedRecharge = 1;
        }
        actGI->nodeList.sort(compare_charge);

        printChargingInfo("BEFORE SECOND CHECKING FEASIBILITY -> ");

        if (!checkScheduleFeasibilityGroup(actGI)) {
            decideRechargeScedulingGroupLast(actGI);
            ris = false;
        }
    }

    printChargingInfo("AFTER Scheduling decision -> ");

    return ris;
}

void UDPBasicRecharge::printDistributedChargingInfo(std::ostream &ss, const char *str) {
    int numberNodes = this->getParentModule()->getVectorSize();
    ss << ((int) simTime().dbl()) << " - ";
    ss << str;
    int nInCharge = 0;
    ss << "{";
    for (int i = 0; i < numberNodes; i++) {
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));
        if (battN->isCharging()) {
            nInCharge++;
            ss << i << " ";
        }
    }
    ss << "| " << nInCharge << "]}|||";
    for (int i = 0; i < numberNodes; i++) {
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));

        ss << "[" << i << "]" << battN->getBatteryLevelAbs() << "|" << battN->isCharging() << "  ";
    }
    ss << endl;
}

void UDPBasicRecharge::printChargingInfo(std::ostream &ss, const char *str) {
    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);
        ss << ((int) simTime().dbl()) << " - ";
        ss << str;
        for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
            nodeAlgo_t *actNO = &(*itn);

            ss << "[" << actNO->addr << "]" << actNO->energy << "|" << actNO->executedRecharge << "/" << actNO->assignedRecharge << "|" << actNO->isCharging << "  ";
        }
        ss << endl;
    }
}

void UDPBasicRecharge::printChargingInfo(const char *str) {
    /*for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        EV << str;
        for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
            nodeAlgo_t *actNO = &(*itn);

            EV << "[" << actNO->addr << "]" << actNO->energy << "|" << actNO->executedRecharge << "/" << actNO->assignedRecharge << "|" << actNO->isCharging << "  ";
        }
        EV << endl;
    }*/
    printChargingInfo(EV, str);
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
        updateBatteryVals(&(actGI->nodeList));
    }

    for (auto it = groupList.begin(); it != groupList.end(); it++) {
        groupInfo_t *actGI = &(*it);

        checkCentralizedRechargeGroup(actGI);

        //checkAliveGroup(actGI);
    }
    printChargingInfo();

    if (printAnalticalLog) {
        FILE *f = fopen(logFile, "a");
        if (f) {
            std::stringstream ss;
            printChargingInfo(ss, "BATTERY INFO -> ");
            fwrite(ss.str().c_str(), 1 , ss.str().length(), f);
            fclose(f);
        }
        else {
            fprintf(stderr, "Error opening file: %s \n", logFile); fflush(stderr);
            perror("Error writing on file");
            error("Error writing on file\n");
        }
    }
}

void UDPBasicRecharge::checkCentralizedRechargeGroup(groupInfo_t *actGI) {
    //updateBatteryVals(&(actGI->nodeList));

    for (auto itn = actGI->nodeList.begin(); itn != actGI->nodeList.end(); itn++){
        nodeAlgo_t *actNO = &(*itn);

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
                    //decideRechargeScedulingGroup(actGI);
                    if (st == ANALYTICAL) {
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
                    else if (st == ROUNDROBIN){
                        actGI->chargingAppAddr = -1;
                        actNO->isCharging = false;
                        putNodeInDischarging(actNO->addr);
                        decideRechargeScedulingGroupRR(actGI);
                        checkCentralizedRechargeGroup(actGI);
                    }
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

double UDPBasicRecharge::getFullCoverage(void) {

    // create the groups
    int activeNodes = 0;
    int numberNodes = this->getParentModule()->getVectorSize();
    std::vector< std::vector<bool> > matrixVal;

    //Grow rows by matrixsideSize
    matrixVal.resize(mob->getConstraintAreaMax().x);
    for(int i = 0 ; i < (int)matrixVal.size() ; ++i) {
        //Grow Columns by matrixsideSize
        matrixVal[i].resize(mob->getConstraintAreaMax().y);
        for(int j = 0 ; j < (int)matrixVal[i].size() ; ++j) {      //modify matrix
            matrixVal[i][j] = false;
        }
    }

    double actArea = 0.0;
    for (int i = 0; i < numberNodes; i++) {
        VirtualSpringMobility *mobN = check_and_cast<VirtualSpringMobility *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("mobility"));
        power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", i)->getSubmodule("battery"));

        if (battN->getState() == power::SimpleBattery::DISCHARGING) {
            activeNodes++;
            for(int i = 0 ; i < (int)matrixVal.size() ; ++i) {
                for(int j = 0 ; j < (int)matrixVal[i].size() ; ++j) {      //modify matrix
                    if (matrixVal[i][j] == false) {
                        Coord point = Coord(i,j);
                        if(point.distance(mobN->getCurrentPosition()) <= sensorRadious) {
                            matrixVal[i][j] = true;
                            actArea += 1.0;
                        }
                    }
                }
            }
        }
    }
    //printMatrix(matrixVal);

    //double actArea = 0.0;
    //for(int i = 0 ; i < (int)matrixVal.size() ; ++i) {
    //    for(int j = 0 ; j < (int)matrixVal[i].size() ; ++j) {      //modify matrix
    //        if (matrixVal[i][j] == true) {
    //            actArea += 1.0;
    //        }
    //    }
    //}
    //double maxArea = ((double) numberNodes) * ((sensorRadious*sensorRadious) * (3.0 / 2.0) * sqrt(3.0));
    //double maxArea = ((double) numberNodes) * ((sensorRadious*sensorRadious) * 2.598076211);
    double maxArea = ((double) (numberNodes - chargingStationNumber)) * ((sensorRadious*sensorRadious) * 2.598076211);

    double ratio = actArea / maxArea;

    if (ratio > 1.0) ratio = 1.0;

    return ratio;
}

double UDPBasicRecharge::getMyCoverageMax(void) {
    double maxCov = 0;
    int matrixsideSize = (sensorRadious * 2.0) + 2.0;
    //std::vector< std::vector<bool> > matrixVal;

    //Grow rows by matrixsideSize
    //matrixVal.resize(matrixsideSize);
    //for(int i = 0 ; i < matrixsideSize ; ++i) {
        //Grow Columns by matrixsideSize
    //    matrixVal[i].resize(matrixsideSize);
    //}

    for(int i = 0 ; i < matrixsideSize ; ++i) {
        for(int j = 0 ; j < matrixsideSize ; ++j) {      //modify matrix
            Coord center = Coord::ZERO;
            Coord point = Coord(i-(sensorRadious+1),j-(sensorRadious+1));

            if(center.distance(point) <= sensorRadious) {
                //matrixVal[i][j] = true;
                maxCov++;
            }
            //else{
                //matrixVal[i][j] = false;
            //}
        }

    }

    return maxCov;
}

void UDPBasicRecharge::printMatrix(std::vector< std::vector<bool> > &matrix) {
    for(int i = 0 ; i < (int)(matrix.size()) ; ++i) {
        for(int j = 0 ; j < (int)(matrix[i].size()) ; ++j) {
            if (matrix[i][j]){
                EV << "1 ";
            }
            else {
                EV << "0 ";
            }
            //EV << matrix[i][j] ? "1 " : "0 ";
        }
        EV << endl;
    }
}

double UDPBasicRecharge::getMyCoverageActual(void) {
    double countCov = 0;
    int matrixsideSize = (sensorRadious * 2.0) + 2.0;
    std::vector< std::vector<bool> > matrixVal;

    //Grow rows by matrixsideSize
    matrixVal.resize(matrixsideSize);
    for(int i = 0 ; i < matrixsideSize ; ++i) {
        //Grow Columns by matrixsideSize
        matrixVal[i].resize(matrixsideSize);
    }

    for(int i = 0 ; i < matrixsideSize ; ++i) {
        for(int j = 0 ; j < matrixsideSize ; ++j) {      //modify matrix
            Coord center = Coord::ZERO;
            Coord point = Coord(i-(sensorRadious+1),j-(sensorRadious+1));
            if(center.distance(point) <= sensorRadious) {
                matrixVal[i][j] = true;
            }
            else{
                matrixVal[i][j] = false;
            }
        }
    }

    //EV << "FULL " << endl;
    //printMatrix(matrixVal);

    for (auto it = neigh.begin(); it != neigh.end(); it++) {
        nodeInfo_t *act = &(it->second);

        for(int i = 0 ; i < matrixsideSize ; ++i) {
            for(int j = 0 ; j < matrixsideSize ; ++j) {      //modify matrix
                Coord node = act->pos - mob->getCurrentPosition();
                Coord point = Coord(i-(sensorRadious+1),j-(sensorRadious+1));
                if(node.distance(point) <= sensorRadious) {
                    matrixVal[i][j] = false;
                }
            }
        }
    }

    //EV << "CUT " << endl;
    //printMatrix(matrixVal);

    for(int i = 0 ; i < matrixsideSize ; ++i) {
        for(int j = 0 ; j < matrixsideSize ; ++j) {
            if(matrixVal[i][j] == true) {
                countCov++;
            }
        }

    }

    EV << "NUM OF POINT COVERED UNIQUE: " << countCov;

    return countCov;

}

} /* namespace inet */
