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
        reinforcementRechargeTime = par("reinforcementRechargeTime").boolValue();
        reinforcementRechargeAlpha = par("reinforcementRechargeAlpha");
        reinforcementRechargeAlphaFinal = par("reinforcementRechargeAlphaFinal");
        chargeTimeOthersNodeFactor = par("chargeTimeOthersNodeFactor");
        makeCoverageLog = par("makeCoverageLog").boolValue();
        //developingStimuli = par("developingStimuli").boolValue();
        constantTheta = par("constantTheta");
        dicountminLINEAR4 = par("dicountminLINEAR4");
        temp_factorProbDischarge = par("temp_factorProbDischarge");
        exponential_dischargeProb_decay = par("exponential_dischargeProb_decay");
        const_c_dischargeProb = par("const_c_dischargeProb");

        //logFile = par("analticalLogFile").str();
        printAnalticalLog = par("printAnalticalLog").boolValue();
        snprintf(logFile, sizeof(logFile), "%s", par("analticalLogFile").stringValue());
        remove(logFile);

        firstRecharge = true;
        lastPosBeforeCharge = rebornPos;
        rechargeLostAccess = 0;
        reinforcementVal = -1;
        inRechargingTime = 0;

        failedAttemptCount = 0;

        startRecharge = simTime();

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

        std::string chargeLengthType = par("chargeLengthType").stdstringValue();
        if (chargeLengthType.compare("MIN") == 0) {
            rlt = MIN_VAL;
        }
        else if (chargeLengthType.compare("MAX") == 0) {
            rlt = MAX_VAL;
        }
        else if (chargeLengthType.compare("AVG") == 0) {
            rlt = AVG_VAL;
        }
        else {
            error("Wrong \"chargeLengthType\" parameter");
        }


        std::string dischargeProbEnergyToUseType = par("dischargeProbEnergyToUse").stdstringValue();
        if (dischargeProbEnergyToUseType.compare("ENERGYMIN") == 0) {
            dischargeProbEnergyToUse = ENERGYMIN;
        }
        else if (dischargeProbEnergyToUseType.compare("ENERGYMAX") == 0) {
            dischargeProbEnergyToUse = ENERGYMAX;
        }
        else if (dischargeProbEnergyToUseType.compare("ENERGYAVG") == 0) {
            dischargeProbEnergyToUse = ENERGYAVG;
        }
        else {
            error("Wrong \"dischargeProbEnergyToUse\" parameter");
        }


        std::string gameTheoryKnowledgeType_str = par("gameTheoryKnowledgeType").stdstringValue();
        if (gameTheoryKnowledgeType_str.compare("LOCAL_KNOWLEDGE") == 0) {
            gameTheoryKnowledgeType = LOCAL_KNOWLEDGE;
        }
        else if (gameTheoryKnowledgeType_str.compare("GLOBAL_KNOWLEDGE") == 0) {
            gameTheoryKnowledgeType = GLOBAL_KNOWLEDGE;
        }
        else {
            error("Wrong \"gameTheoryKnowledgeType\" parameter");
        }



        std::string stimType = par("stimulusType").stdstringValue();
        //"CONST_C", "VAR_C_P0", "VAR_C_VAR_P"
        if (stimType.compare("STIM_OLD") == 0) {
            stim_type = STIM_OLD;
        }
        else if (stimType.compare("CONST_C") == 0) {
            stim_type = CONST_C;
        }
        else if (stimType.compare("VAR_C_P1") == 0) {
            stim_type = VAR_C_P1;
        }
        else if (stimType.compare("VAR_C_VAR_P") == 0) {
            stim_type = VAR_C_VAR_P;
        }
        else {
            error("Wrong \"stimulusType\" parameter");
        }

        std::string constType = par("varConstantType").stdstringValue();
        //"SIGMOID", "LINEAR1", "LINEAR2"
        if (constType.compare("SIGMOID") == 0) {
            constant_type = SIGMOID;
        }
        else if (constType.compare("LINEAR1") == 0) {
            constant_type = LINEAR1;
        }
        else if (constType.compare("LINEAR2") == 0) {
            constant_type = LINEAR2;
        }
        else if (constType.compare("LINEAR3") == 0) {
            constant_type = LINEAR3;
        }
        else if (constType.compare("LINEARDISCOUNT") == 0) {
            constant_type = LINEARDISCOUNT;
        }
        else if (constType.compare("SIGMOIDDISCOUNT") == 0) {
            constant_type = SIGMOIDDISCOUNT;
        }
        else if (constType.compare("LINEARINCREASE") == 0) {
            constant_type = LINEARINCREASE;
        }
        else if (constType.compare("SIGMOIDINCREASE") == 0) {
            constant_type = SIGMOIDINCREASE;
        }
        else {
            error("Wrong \"varConstantType\" parameter");
        }

        std::string probType = par("varProbabilityType").stdstringValue();
        //"ONE_OVER_FORMULAPAPER"
        if (probType.compare("ONE_OVER_FORMULAPAPER") == 0) {
            probability_type = ONE_OVER_FORMULAPAPER;
        }
        else {
            error("Wrong \"varProbabilityType\" parameter");
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
        failedAttemptVector.setName("FailedAttemptVal");
        dischargeProbVector.setName("DischargeProbVal");

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


            rechargeLostAccess = 0;

            //inRechargingTime += cTime;
            startRecharge = simTime();

            lastPosBeforeCharge = mob->getCurrentPosition();

            sb->setState(power::SimpleBattery::CHARGING);

            if ((st == STIMULUS) && (stim_type != STIM_OLD)) {
                // do nothing
                inRechargingTime += checkRechargeTimer;
            }
            else {
                double cTime = calculateRechargeTime(true);

                inRechargingTime += cTime;

                fprintf(stderr, "[%d] - Going in charging: %f\n", myAppAddr, cTime);fflush(stderr);

                scheduleAt(simTime() + cTime, dischargeTimer);
            }

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
            failedAttemptCount++;

            sb->setDoubleSwapPenality();
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
            if (sb->getState() == power::SimpleBattery::CHARGING) {
                inRechargingTime += checkRechargeTimer;
            }
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
        if (makeCoverageLog){
            totalCoverageVector.record(getFullCoverage());
        }
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
        if (st == STIMULUS) {
            stimulusVector.record(calculateRechargeStimuli());
        }
        else if (st == PROBABILISTIC) {
            stimulusVector.record(1 - sb->getBatteryLevelPercInitial());
        }

        thresholdVector.record(calculateRechargeThreshold());
        responseVector.record(calculateRechargeProb(true));
        degreeVector.record(calculateNodeDegree());
        timeFactorVector.record(calculateRechargeStimuliTimeFactor());
        energyFactorVector.record(calculateRechargeStimuliEnergyFactor());
        failedAttemptVector.record(failedAttemptCount);

        if (sb->isCharging()) {
            dischargeProbVector.record(calculateNodeDischargeProb());
        }
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

                //lastRechargeTimestamp += calculateRechargeTime(false) / ((double) chargingStationNumber);
                //if (lastRechargeTimestamp > simTime())
                //    lastRechargeTimestamp = simTime();

                firstRecharge = false;

                if (reinforcementVal >= 0) {
                    reinforcementVal = (reinforcementRechargeAlpha * aPkt->getGoingToRechargeTime()) +
                            ((1.0 - reinforcementRechargeAlpha) * reinforcementVal);
                }
                else {
                    reinforcementVal = aPkt->getGoingToRechargeTime();
                }
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
                node->inRechargeT = aPkt->getInRecharge();
                node->gameTheoryC = aPkt->getGameTheoryC();
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
    payload->setInRecharge(inRechargingTime);
    payload->setGoingToRecharge(true);
    payload->setGoingToRechargeTime(calculateRechargeTime(false));
    payload->setGameTheoryC(getGameTheoryC());

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
        payload->setInRecharge(inRechargingTime);
        payload->setGoingToRecharge(false);
        payload->setGoingToRechargeTime(0);
        payload->setGameTheoryC(getGameTheoryC());

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

double UDPBasicRecharge::calculateRechargeProb(bool useDishargeProbIfTheCase) {

    if (st == STIMULUS) {
        double ris = 0;
        int numberNodes = this->getParentModule()->getVectorSize();
        double s, c, t;
        double nmeno1SquareRoot, unomenoCi;
        long double produttoria;
        long double dischargeP = -1;

        //if (rechargeLostAccess > 0) {
        //    return 1;
        //}
        //else {
        switch (stim_type) {
            case STIM_OLD:
                s = calculateRechargeStimuli();
                t = calculateRechargeThreshold();

                ris = pow(s, stimulusExponent) / (pow(s, stimulusExponent) + pow(t, stimulusExponent));

                break;
            case CONST_C:
                //c = (2.0 * sb->getSwapLoose()) / (sb->getDischargingFactor(checkRechargeTimer) + (sb->getChargingFactor(checkRechargeTimer)));
                c = (getGamma()+getTheta()) / (getAlpha() + getBeta());
                c = c / const_c_dischargeProb;
                s = pow(c, 1.0 / (((double) numberNodes) - 1.0));

                ris = 1.0 - s;

                //fprintf(stderr, "DEVSTIM: calculateRechargeProb_STATIC, c: %lf; s: %lf -> %lf\n", c, s, ris); fflush(stderr);

                break;
            case VAR_C_P1:
            case VAR_C_VAR_P:
                produttoria = 1.0;

                if (stim_type == VAR_C_VAR_P) {
                    if (useDishargeProbIfTheCase) {
                        if (gameTheoryKnowledgeType == GLOBAL_KNOWLEDGE) {
                            for (int j = 0; j < numberNodes; j++) {
                                UDPBasicRecharge *hostj = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("udpApp", 0));
                                power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));

                                if (battN->getState() == power::SimpleBattery::CHARGING) {
                                    dischargeP = hostj->calculateNodeDischargeProb();
                                }
                            }
                        }
                        else if (gameTheoryKnowledgeType == LOCAL_KNOWLEDGE) {
                            // use mine
                            dischargeP = calculateNodeDischargeProb();
                        }
                        else {
                            error("Wrong knowledge scope");
                        }
                    }
                }
                unomenoCi = 1.0 - getGameTheoryC();

                if (unomenoCi > 0){

                    if (gameTheoryKnowledgeType == GLOBAL_KNOWLEDGE) {
                        for (int j = 0; j < numberNodes; j++) {
                            UDPBasicRecharge *hostj = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("udpApp", 0));
                            //power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));
                            double hostC = hostj->getGameTheoryC();
                            long double ppp = (1.0 - hostC) / unomenoCi;
                            produttoria = produttoria * ppp;

                            //fprintf(stderr, "%Lf ", ppp);

                        }
                    }
                    else if (gameTheoryKnowledgeType == LOCAL_KNOWLEDGE) {
                        for (auto it = neigh.begin(); it != neigh.end(); it++) {
                            nodeInfo_t *act = &(it->second);

                            double hostC = act->gameTheoryC;
                            long double ppp = (1.0 - hostC) / unomenoCi;
                            produttoria = produttoria * ppp;
                        }
                    }
                    else {
                        error("Wrong knowledge scope");
                    }

                    //fprintf(stderr, "\n");
                    if (dischargeP > 0) {
                        produttoria = produttoria * (unomenoCi / dischargeP);
                    }
                    else {
                        produttoria = produttoria * unomenoCi;
                    }

                    if (gameTheoryKnowledgeType == GLOBAL_KNOWLEDGE) {
                        nmeno1SquareRoot = powl(produttoria, 1.0 / (((long double) numberNodes) - 1.0));
                    }
                    else if (gameTheoryKnowledgeType == LOCAL_KNOWLEDGE){
                        nmeno1SquareRoot = powl(produttoria, 1.0 / ((long double) neigh.size()));
                    }
                    else {
                        error("Wrong knowledge scope");
                    }

                    s = nmeno1SquareRoot;

                    if (s > 1) s = 1;
                    if (s < 0) s = 0;
                }
                else {
                    s = 0.0;
                }

                /*
                for (int j = 0; j < numberNodes; j++) {
                    UDPBasicRecharge *hostj = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("udpApp", 0));
                    power::SimpleBattery *battN = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));
                    double hostC = hostj->getGameTheoryC();
                    long double ppp = 1.0 - hostC;
                    produttoria = produttoria * ppp;

                    //fprintf(stderr, "DEVSTIM: produttoriaTMP: %Lf; ppp: %Lf\n", produttoria, ppp); fflush(stderr);
                }
                if (dischargeP > 0) {
                    nmeno1SquareRoot = powl(produttoria * (1.0 / dischargeP), 1.0 / (((long double) numberNodes) - 1.0));
                }
                else {
                    nmeno1SquareRoot = powl(produttoria, 1.0 / (((long double) numberNodes) - 1.0));
                }

                if (unomenoCi > 0){
                    s = nmeno1SquareRoot / unomenoCi;
                    if (s > 1) s = 1;
                    if (s < 0) s = 0;
                }
                else {
                    s = 1.0;
                }
                */

                ris = 1.0 - s;

                //fprintf(stderr, "DEVSTIM. NN: %d; produttoria: %Lf; nmeno1S: %lf; myC: %lf; unomenoCi: %lf; s: %lf\n",
                //        numberNodes, produttoria, nmeno1SquareRoot, getGameTheoryC(), unomenoCi, s); fflush(stderr);

                //fprintf(stderr, "DEVSTIM: calculateRechargeProb_DYNAMIC, ris: %lf\n\n", ris); fflush(stderr);

                break;

            default:
                break;
        }

        return ris;


        /*if (developingStimuli) {
            int numberNodes = this->getParentModule()->getVectorSize();
            double ris = 0;

            if (stim_type ==  CONST_C) {
                double c = (2.0 * sb->getSwapLoose()) / (sb->getDischargingFactor(checkRechargeTimer) + (sb->getChargingFactor(checkRechargeTimer)));
                double s = pow(c, 1.0 / (((double) numberNodes) - 1.0));

                fprintf(stderr, "DEVSTIM: calculateRechargeProb_STATIC, c: %lf; s: %lf -> %lf\n", c, s, (1.0 - s)); fflush(stderr);
            }
            else if (stim_type ==  VAR_C_P1) {

            }

            double produttoria = 1.0;

            for (int j = 0; j < numberNodes; j++) {
                UDPBasicRecharge *hostj = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("udpApp", 0));

                produttoria *= (1.0 - hostj->getGameTheoryC());
            }
            double nmeno1SquareRoot = pow(produttoria, 1.0 / (((double) numberNodes) - 1.0));
            double unomenoCi = 1.0 - getGameTheoryC();

            double s = nmeno1SquareRoot / unomenoCi;

            fprintf(stderr, "DEVSTIM: calculateRechargeProb_DYNAMIC, ris: %lf\n", (1.0 - s)); fflush(stderr);

            return (1.0 - s);

            return ris;
        }
        else {
            double stim = calculateRechargeStimuli();
            double tetha = calculateRechargeThreshold();

            return (pow(stim, stimulusExponent) / (pow(stim, stimulusExponent) + pow(tetha, stimulusExponent)));
        }*/
    }
    else {//if (st == PROBABILISTIC) {
        return (1 - sb->getBatteryLevelPercInitial());
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

    rechargeEstimation = calculateRechargeTime(false);

    //double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * timeFactorMultiplier);
    //double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * ((double) filteredNeigh.size()));
    //timeFactor = timeFactor / timeFactorMultiplier; //TODO

    double timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * ((double) filteredNeigh.size()));
    if (stationANDnodeKNOWN) {
        int numberNodes = this->getParentModule()->getVectorSize();
        timeFactor = (simTime() - lastRechargeTimestamp).dbl() / (rechargeEstimation * (((double) numberNodes) / ((double) chargingStationNumber)));
    }
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

    if (stim_type != STIM_OLD) {

        //TODO per ora ritorno un valore costante
        return 0;
    }
    else {
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
    double ris = 0;

    if (stim_type != STIM_OLD) {
        //TODO per ora ritorno un valore costante
        ris = 0;
    }
    else {
        ris = ((double) myDegree) / 6.0;
        if (ris > 1) ris = 1;
    }
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

double UDPBasicRecharge::reinforceTimeVal(double val) {
    double ris = val;

    if (reinforcementVal > 0) {
        ris = (reinforcementRechargeAlphaFinal * val) + ((1.0 - reinforcementRechargeAlphaFinal) * reinforcementVal);
    }

    return ris;
}

double UDPBasicRecharge::calculateChargeDiff (double myChoice) {
    double ris = myChoice;

    if (neigh.size() > 0) {
        std::map<int, nodeInfo_t> filteredNeigh;
        getFilteredNeigh(filteredNeigh);

        double maxC = -1;
        for (auto it = filteredNeigh.begin(); it != filteredNeigh.end(); it++) {
            nodeInfo_t *act = &(it->second);

            if (act->inRechargeT > maxC) {
                maxC = act->inRechargeT;
            }
        }

        if ((maxC > 0) && (maxC > (inRechargingTime + myChoice))) {
            ris = maxC - inRechargingTime;
        }
    }

    return ris;
}

double UDPBasicRecharge::calculateSwapPenalitiesEstimationCount(double estimatedSteps) {
    double ris = 0;
    int nSteps = estimatedSteps * ((((double) this->getParentModule()->getVectorSize()) / ((double) chargingStationNumber)) - 1.0);

    for (int i = 1; i <= nSteps; i++) {
        ris += ((double) i) / ((double) nSteps);
    }

    return ris;
}

double UDPBasicRecharge::calculateRechargeTime(bool log) {

    double recTime = 0;
    std::stringstream ss;

    if (st == STIMULUS) {

        if (stim_type != STIM_OLD) {
            return checkRechargeTimer;
        }
        else {

            // default charge until full charge
            //double tt = ((sb->getFullCapacity() - sb->getBatteryLevelAbs()) / sb->getChargingFactor(checkRechargeTimer)) * checkRechargeTimer;
            double tt = numRechargeSlotsStimulusZeroNeigh * checkRechargeTimer;

            if (log) ss << "RECHARGETIME STIMULUS: Default charge time: " << tt << " - Neigh size: " << neigh.size() << endl;

            if (neigh.size() > 0) {

                std::map<int, nodeInfo_t> filteredNeigh;
                getFilteredNeigh(filteredNeigh);

                double sumE = sb->getBatteryLevelAbs();
                double maxE = sb->getBatteryLevelAbs();
                double minE = sb->getBatteryLevelAbs();
                if (log) ss << "RECHARGETIME STIMULUS: my battery: " << sb->getBatteryLevelAbs() << endl;
                for (auto it = filteredNeigh.begin(); it != filteredNeigh.end(); it++) {
                //for (auto it = neigh.begin(); it != neigh.end(); it++) {
                    nodeInfo_t *act = &(it->second);

                    sumE += act->batteryLevelAbs;
                    if (act->batteryLevelAbs > maxE)
                        maxE = act->batteryLevelAbs;
                    if (act->batteryLevelAbs < minE)
                        minE = act->batteryLevelAbs;

                    if (log) ss << "RECHARGETIME STIMULUS: others battery: " << act->batteryLevelAbs << endl;

                    //TODO remove
                    //break;
                }

                //double averageE = sumE / (((double) neigh.size()) + 1.0);
                //double averageE = sumE / (((double) nodeFiltered.size()) + 1.0);
                double averageE = sumE / (((double) filteredNeigh.size()) + 1.0);

                if (log) ss << "RECHARGETIME STIMULUS: Average Energy: " << averageE << ", Max Energy: " << maxE
                //if (log) ss << "RECHARGETIME STIMULUS: Average Energy: " << averageE << ", Max Energy: " << maxE
                        << " - Discharging Factor: " << sb->getDischargingFactor(checkRechargeTimer)
                        << " - SwapLoose Factor: " << sb->getSwapLoose()
                        << " - NodeFiltered size: " << filteredNeigh.size()
                        << " - Neigh size: " << neigh.size()
                        //<< " - stationANDnodeKNOWN: " << stationANDnodeKNOWN
                        //<< " - reinforcementRechargeTime: " << reinforcementRechargeTime
                        << " - checkRechargeTimer: " << checkRechargeTimer
                        << endl;

                //double numSteps = (maxE - (2.0 * sb->getSwapLoose())) / sb->getDischargingFactor(checkRechargeTimer);
                double valToUse = 1;
                switch (rlt) {
                case MIN_VAL:
                    valToUse = minE;
                    break;

                case MAX_VAL:
                    valToUse = maxE;
                    break;

                case AVG_VAL:
                    valToUse = averageE;
                    break;

                default:
                    error("Wring rlt value");
                    break;
                }
                double tmpnumSteps = (valToUse - (2.0 * sb->getSwapLoose())) / sb->getDischargingFactor(checkRechargeTimer);
                double swapPenalitiesEstimation = calculateSwapPenalitiesEstimationCount(tmpnumSteps/checkRechargeTimer) * (2.0 * sb->getSwapLoose());
                double numSteps = (valToUse - (2.0 * sb->getSwapLoose()) - swapPenalitiesEstimation) / sb->getDischargingFactor(checkRechargeTimer);

                if (log) ss << "RECHARGETIME STIMULUS: swapPenalitiesEstimation: " << swapPenalitiesEstimation <<
                        ", tmpnumSteps: " << tmpnumSteps <<
                        ", numSteps: " << numSteps << endl;

                //int actualNeigh = neigh.size();
                //if (actualNeigh > 7) actualNeigh = 7;
                //double numChargeSlots = numSteps / ((double) neigh.size());
                //double numChargeSlots = numSteps / ((double) actualNeigh);
                //double numChargeSlots = numSteps / ((double) nodeFiltered.size());
                double numChargeSlots;
                if (stationANDnodeKNOWN) {
                    int numberNodes = this->getParentModule()->getVectorSize();
                    numChargeSlots = numSteps / ((((double) numberNodes) / ((double) chargingStationNumber)) - 1.0);

                    if (log) ss << "RECHARGETIME STIMULUS: numSteps: " << numSteps <<
                            ", numberNodes: " << numberNodes <<
                            ", chargingStationNumber: " << chargingStationNumber <<
                            ", n-1: " << ((((double) numberNodes) / ((double) chargingStationNumber)) - 1.0) <<
                            ", numChargeSlots: " << numChargeSlots <<
                            endl;
                }
                else {
                    numChargeSlots = numSteps / ((double) filteredNeigh.size() + 1.0);
                }
                //double numChargeSlots = numSteps / ((double) filteredNeigh.size() + 1.0);
                tt = numChargeSlots * checkRechargeTimer;

                //tt = (averageE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;
                //tt = (maxE / ((sb->getDischargingFactor(checkRechargeTimer)) * ((double) neigh.size()))) * checkRechargeTimer;

                if (reinforcementRechargeTime) {
                    tt = reinforceTimeVal(tt);
                }

                if (chargeTimeOthersNodeFactor > 0) {
                    double diffOthers = calculateChargeDiff(tt);
                    if (diffOthers > 0) {
                        if (log) ss << "RECHARGETIME STIMULUS: my tt: " << tt << endl;
                        tt = (chargeTimeOthersNodeFactor * diffOthers) + ((1.0 - chargeTimeOthersNodeFactor) * tt);
                        if (log) ss << "RECHARGETIME STIMULUS: diffOthers: " << diffOthers
                                << " chargeTimeOthersNodeFactor: " << chargeTimeOthersNodeFactor
                                << endl;
                    }
                }
            }

            recTime = tt;
        }
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
    double cw = 3;
    double ris = 0;

    //cw = cw / (((double)rechargeLostAccess) + 1.0);
    cw = cw / pow(2.0, ((double) rechargeLostAccess));

    ris = cw * dblrand();
    ris += 0.01;

    return ris;

}

void UDPBasicRecharge::checkRecharge(void) {
    double prob = calculateRechargeProb(true);

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
    //else {
    //    rechargeLostAccess = 0;
    //}


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

double UDPBasicRecharge::calculateNodeDischargeProb(void) {
    if (sb->isFull()) {
        return 1;
    }
    else {
        if ((st == STIMULUS) && (stim_type != STIM_OLD)) {
            double ris = 1;

            if (stim_type == VAR_C_VAR_P) {

                int numberNodes = this->getParentModule()->getVectorSize();
                double estimatedTimeInRecharging;
                //double energyToUse = getEavg(true);
                //double energyToUse = getEavg(false);
                //double energyToUse = getEmin(false);
                double energyToUse;
                //bool isThereAnyCharging = false;
                double timeCalcNum, timeCalcDen1, timeCalcDen2;
                double gPLUSt = getGamma() + getTheta();

                energyToUse = sb->getBatteryLevelAbs();
                switch (dischargeProbEnergyToUse) {
                case ENERGYMIN:
                default:
                    energyToUse = getEmin(false, GLOBAL_KNOWLEDGE);
                    break;
                case ENERGYMAX:
                    energyToUse = getEmax(false, GLOBAL_KNOWLEDGE);
                    break;
                case ENERGYAVG:
                    energyToUse = getEavg(false, GLOBAL_KNOWLEDGE);
                    break;
                }

                /*
                for (int j = 0; j < numberNodes; j++) {
                    power::SimpleBattery *hostjsb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));

                    if (hostjsb->isCharging()) {
                        isThereAnyCharging = true;
                        break;
                    }
                }

                if (isThereAnyCharging) {
                    estimatedTimeInRecharging = (energyToUse - getGamma() - getTheta()) / (getAlpha() * ((double)(numberNodes - 1.0)));
                }
                else {
                    estimatedTimeInRecharging = (energyToUse - getGamma() - getTheta()) / (getAlpha() * ((double)(numberNodes)));
                }
                */

                timeCalcNum = energyToUse - gPLUSt;
                timeCalcDen1 = getAlpha() * ((double)(numberNodes - 1.0));
                timeCalcDen2 = 0.0;

                if (true) {
                    //timeCalcDen2 = ((double)(numberNodes - 1.0)) * (getRechargeProbMax(false) * (getGamma() + getTheta()));// / checkRechargeTimer;
                    //timeCalcDen2 = ((double)(numberNodes - 1.0)) * ((getGamma() + getTheta()));// / checkRechargeTimer;
                    timeCalcDen2 = ((double)(numberNodes - 1.0)) * gPLUSt;// / checkRechargeTimer;
                    //timeCalcDen2 = (1 * (getGamma() + getTheta()));// / checkRechargeTimer;
                    //timeCalcDen2 = (0 * (getGamma() + getTheta())) / checkRechargeTimer;

                    //fprintf(stderr, "timeCalcDen1: %lf and timeCalcDen2 = %lf\n", timeCalcDen1, timeCalcDen2); fflush(stderr);
                }

                estimatedTimeInRecharging = timeCalcNum / (timeCalcDen1 + timeCalcDen2);
                //estimatedTimeInRecharging = (energyToUse - getGamma() - getTheta()) / (getAlpha() * ((double)(numberNodes)));
                //estimatedTimeInRecharging = estimatedTimeInRecharging * checkRechargeTimer;
                //estimatedTimeInRecharging = estimatedTimeInRecharging / temp_factorProbDischarge;

                if (exponential_dischargeProb_decay == 0) {
                    estimatedTimeInRecharging = estimatedTimeInRecharging / temp_factorProbDischarge;
                    ris = 1.0 / estimatedTimeInRecharging;
                }
                else {
                    double timeInCharge = (simTime() - startRecharge).dbl();

                    estimatedTimeInRecharging = estimatedTimeInRecharging * checkRechargeTimer;

                    //fprintf(stderr, "timeInCharge: %lf and estimatedTimeInRecharging = %lf\n", timeInCharge, estimatedTimeInRecharging); fflush(stderr);

                    if (timeInCharge >= estimatedTimeInRecharging){
                        ris = 1.0;
                    }
                    else {
                        ris = pow(timeInCharge / estimatedTimeInRecharging, exponential_dischargeProb_decay);
                    }
                }

            }
            else if (stim_type == CONST_C){
                ris = const_c_dischargeProb;
            }
            else {
                ris = 1.0;
            }

            //fprintf(stderr, "calculateNodeDischargeProb = %lf\n", ris); fflush(stderr);

            if (ris < 0) ris = 0;
            if (ris > 1) ris = 1;

            return ris;
        }
        else {
            // in the older case, this function is not used for discharging method
            return 0;
        }
    }
}

void UDPBasicRecharge::checkDischarge(void) {
    double prob = 0;

    if ((st == STIMULUS) && (stim_type != STIM_OLD)) {
        prob = calculateNodeDischargeProb();
    }
    else {
        prob = calculateDischargeProb();
    }

    if (dblrand() <= prob) {
        sb->setState(power::SimpleBattery::DISCHARGING);

        //stop the node
        neigh.clear();
        //mob->clearVirtualSprings();
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

double UDPBasicRecharge::getGameTheoryC_Sigmoid(void) {
    double ris = 0;
    double eavg = getEavg(false, gameTheoryKnowledgeType);
    //double e = pow (eavg / sb->getBatteryLevelAbs(), 2.0);
    double e = exp (1.0 - (sb->getBatteryLevelAbs() / eavg) );
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();

    ris = ((a / e) + b - t - g) / (a + b);

    return ris;
}

double UDPBasicRecharge::getGameTheoryC_Linear1(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double e = (eMAX - sb->getBatteryLevelAbs()) / eMAX;

    ris = ((a / e) + b - t - g) / (a + b);

    return ris;
}

double UDPBasicRecharge::getGameTheoryC_Linear2(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double eMIN = getEmin(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double e = (sb->getBatteryLevelAbs() - eMIN) / (eMAX - eMIN);

    ris = (a + b - t - g) / ((e * (a + t + g)) + b - t - g);

    return ris;
}

double UDPBasicRecharge::getGameTheoryC_Linear3(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double e = (eMAX - sb->getBatteryLevelAbs()) / eMAX;

    ris = 1.0 - ( (t + g) / (e * (a + b)) );

    return ris;
}

double UDPBasicRecharge::getGameTheoryC_LinearDiscount(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double eMIN = getEmin(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double myE = sb->getBatteryLevelAbs();
    //double e = (((myE - eMIN) / (eMAX - eMIN)) * (1.0 - dicountminLINEAR4)) + dicountminLINEAR4;
    double e = (((eMAX - myE) / (eMAX - eMIN)) * (1.0 - dicountminLINEAR4)) + dicountminLINEAR4;

    if ((eMAX - eMIN) == 0) {
        e = 1;
    }

    ris = (a + b - t - g) / ((e * (a + t + g)) + b - t - g);

    //fprintf(stderr, "getGameTheoryC_LinearDiscount: %lf; alpha: %lf; beta: %lf; gamma: %lf; theta: %lf; myE: %lf; eMIN: %lf; eMAX: %lf; e: %lf\n",
    //        ris, a, b, g, t, myE, eMIN, eMAX, e); fflush(stderr);

    return ris;
}


double UDPBasicRecharge::getGameTheoryC_SigmoidDiscount(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double sig1, sig2;
    double div = eMAX / 50.0;

    if (div < 1.0) div = 1.0;

    sig1 = (0.5 / (1.0 + exp(((eMAX / 4.0)-sb->getBatteryLevelAbs())/div)));
    sig2 = (0.5 / (1.0 + exp((((eMAX * 3.0) / 4.0)-sb->getBatteryLevelAbs())/div)));

    double sig = 1.0 - (sig1 + sig2);
    double e = (sig * (1.0 - dicountminLINEAR4)) + dicountminLINEAR4;

    ris = (a + b - t - g) / ((e * (a + t + g)) + b - t - g);

    return ris;
}

double UDPBasicRecharge::getGameTheoryC_LinearIncrease(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double eMIN = getEmin(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double myE = sb->getBatteryLevelAbs();
    //double e = (((myE - eMIN) / (eMAX - eMIN)) * (1.0 - dicountminLINEAR4)) + 1.0;
    double e = (((eMAX - myE) / (eMAX - eMIN)) * (1.0 - dicountminLINEAR4)) + 1.0;

    if ((eMAX - eMIN) == 0) {
        e = 1;
    }

    ris = (a + b - t - g) / ((e * (a + t + g)) + b - t - g);

    //fprintf(stderr, "getGameTheoryC_LinearDiscount: %lf; alpha: %lf; beta: %lf; gamma: %lf; theta: %lf; myE: %lf; eMIN: %lf; eMAX: %lf; e: %lf\n",
    //        ris, a, b, g, t, myE, eMIN, eMAX, e); fflush(stderr);

    return ris;
}


double UDPBasicRecharge::getGameTheoryC_SigmoidIncrease(void) {
    double ris = 0;
    double eMAX = getEmax(false, gameTheoryKnowledgeType);
    double a = getAlpha();
    double b = getBeta();
    double t = getTheta();
    double g = getGamma();
    double sig1, sig2;
    double div = eMAX / 50.0;

    if (div < 1.0) div = 1.0;

    sig1 = (0.5 / (1.0 + exp(((eMAX / 4.0)-sb->getBatteryLevelAbs())/div)));
    sig2 = (0.5 / (1.0 + exp((((eMAX * 3.0) / 4.0)-sb->getBatteryLevelAbs())/div)));

    double sig = sig1 + sig2;
    double e = (sig * (1.0 - dicountminLINEAR4)) + 1.0;

    ris = (a + b - t - g) / ((e * (a + t + g)) + b - t - g);

    return ris;
}


double UDPBasicRecharge::getGameTheoryC(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        switch (constant_type) {
        case SIGMOID:
        default:
            ris = getGameTheoryC_Sigmoid();
            break;
        case LINEAR1:
            ris = getGameTheoryC_Linear1();
            break;
        case LINEAR2:
            ris = getGameTheoryC_Linear2();
            break;
        case LINEAR3:
            ris = getGameTheoryC_Linear3();
            break;
        case LINEARDISCOUNT:
            ris = getGameTheoryC_LinearDiscount();
            break;
        case SIGMOIDDISCOUNT:
            ris = getGameTheoryC_SigmoidDiscount();
            break;
        case LINEARINCREASE:
            ris = getGameTheoryC_LinearIncrease();
            break;
        case SIGMOIDINCREASE:
            ris = getGameTheoryC_SigmoidIncrease();
            break;
        }
        /*
        double eavg = getEavg();
        //double e = pow (eavg / sb->getBatteryLevelAbs(), 2.0);
        double e = exp (1.0 - (sb->getBatteryLevelAbs() / eavg) );
        double a = getAlpha();
        double b = getBeta();
        double t = getTheta();
        double g = getGamma();

        ris = ((a / e) + b - t - g) / (a + b);
        */

        //fprintf(stderr, "DEVSTIM: myE: %lf; avgE: %lf; ratio: %lf\n", sb->getBatteryLevelAbs(), eavg, eavg / sb->getBatteryLevelAbs()); fflush(stderr);
        //fprintf(stderr, "DEVSTIM: getGameTheoryC - e^2: %lf theta: %lf; gamma: %lf; alpha: %lf; beta: %lf; ris: %lf\n", e, t, g, a, b, ris); fflush(stderr);
    }
    else {
        double d = (getTheta() + getGamma()) / (getAlpha() + getBeta());

        //fprintf(stderr, "DEVSTIM: theta: %lf; gamma: %lf; alpha: %lf; beta: %lf\n", getTheta(), getGamma(), getAlpha(), getBeta()); fflush(stderr);

        ris = 1.0 - d;
    }

    if (ris > 1) ris = 1;
    if (ris < 0) ris = 0;

    return ris;
}

double UDPBasicRecharge::getTheta(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        ris = sb->getSwapLoose();
    }
    else {
        ris = sb->getSwapLoose();
    }

    return ris;
}

double UDPBasicRecharge::getGamma(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        ris = sb->getSwapLoose();
    }
    else {
        ris = sb->getSwapLoose();
    }

    return ris;
}

double UDPBasicRecharge::getAlpha(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        ris = sb->getDischargingFactor(checkRechargeTimer);
    }
    else {
        ris = sb->getDischargingFactor(checkRechargeTimer);
    }

    return ris;
}

double UDPBasicRecharge::getBeta(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        ris = sb->getChargingFactor(checkRechargeTimer);
    }
    else {
        ris = sb->getChargingFactor(checkRechargeTimer);
    }

    return ris;
}

double UDPBasicRecharge::getP(void) {
    double ris = 0;

    if ( (stim_type == VAR_C_P1) || (stim_type == VAR_C_VAR_P) ) {
        ris = 1;
    }
    else {
        ris = 1;
    }

    return ris;
}

double UDPBasicRecharge::getEavg(bool activeOnly, GameTheoryKnowledge_Type scope) {
    int numberNodes = this->getParentModule()->getVectorSize();
    double sum = 0;
    double nn;
    if (scope == GLOBAL_KNOWLEDGE) {
        for (int j = 0; j < numberNodes; j++) {
            power::SimpleBattery *hostjsb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));

            if (activeOnly && hostjsb->isCharging()) continue;

            sum += hostjsb->getBatteryLevelAbs();
        }
        nn = numberNodes;
    }
    else if (scope == LOCAL_KNOWLEDGE){
        sum = sb->getBatteryLevelAbs();
        for (auto it = neigh.begin(); it != neigh.end(); it++) {
            nodeInfo_t *act = &(it->second);
            sum += act->batteryLevelAbs;
        }
        nn = neigh.size() + 1;
    }
    else {
        error("Wrong knowledge scope");
    }

    return (sum / nn);
}

double UDPBasicRecharge::getEmax(bool activeOnly, GameTheoryKnowledge_Type scope) {
    int numberNodes = this->getParentModule()->getVectorSize();
    //double max = 0;
    double max = sb->getBatteryLevelAbs();
    if (scope == GLOBAL_KNOWLEDGE) {
        for (int j = 0; j < numberNodes; j++) {
            power::SimpleBattery *hostjsb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));

            if (activeOnly && hostjsb->isCharging()) continue;

            if (hostjsb->getBatteryLevelAbs() > max){
                max = hostjsb->getBatteryLevelAbs();
            }
        }
    }
    else if (scope == LOCAL_KNOWLEDGE){
        for (auto it = neigh.begin(); it != neigh.end(); it++) {
            nodeInfo_t *act = &(it->second);
            double actBatt = act->batteryLevelAbs;
            if (actBatt > max) {
                max = actBatt;
            }
        }
    }
    else {
        error("Wrong knowledge scope");
    }

    return max;
}

double UDPBasicRecharge::getEmin(bool activeOnly, GameTheoryKnowledge_Type scope) {
    int numberNodes = this->getParentModule()->getVectorSize();
    //double min = 1000000000;
    double min = sb->getBatteryLevelAbs();
    if (scope == GLOBAL_KNOWLEDGE) {
        for (int j = 0; j < numberNodes; j++) {
            power::SimpleBattery *hostjsb = check_and_cast<power::SimpleBattery *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("battery"));

            if (activeOnly && hostjsb->isCharging()) continue;

            if (hostjsb->getBatteryLevelAbs() < min){
                min = hostjsb->getBatteryLevelAbs();
            }
        }
    }
    else if (scope == LOCAL_KNOWLEDGE){
        for (auto it = neigh.begin(); it != neigh.end(); it++) {
            nodeInfo_t *act = &(it->second);
            double actBatt = act->batteryLevelAbs;
            if (actBatt < min) {
                min = actBatt;
            }
        }
    }
    else {
        error("Wrong knowledge scope");
    }

    return min;
}

double UDPBasicRecharge::getRechargeProbMax(bool useDishargeProbIfTheCase) {
    int numberNodes = this->getParentModule()->getVectorSize();
    double max = 0;
    for (int j = 0; j < numberNodes; j++) {
        UDPBasicRecharge *nodeN = check_and_cast<UDPBasicRecharge *>(this->getParentModule()->getParentModule()->getSubmodule("host", j)->getSubmodule("udpApp", 0));
        double val = nodeN->calculateRechargeProb(useDishargeProbIfTheCase);

        if (val > max){
            max = val;
        }
    }

    return max;
}

} /* namespace inet */
