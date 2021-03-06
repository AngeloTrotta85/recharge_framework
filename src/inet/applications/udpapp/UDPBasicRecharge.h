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

    typedef enum {
        ANALYTICAL,
        ROUNDROBIN,
        STIMULUS,
        PROBABILISTIC
    } Scheduling_Type;

    typedef enum {
        STIM_OLD,
        CONST_C,
        VAR_C_P1,
        VAR_C_VAR_P
    } Stimulus_Type;

    typedef enum {
        SIGMOID,
        LINEAR1,
        LINEAR2,
        LINEAR3,
        LINEARDISCOUNT,
        SIGMOIDDISCOUNT,
        LINEARINCREASE,
        SIGMOIDINCREASE,
        LINEARINCREASECONSISTENT1
    } VarConstant_Type;

    typedef enum {
        ONE_OVER_FORMULAPAPER
    } VarProbability_Type;

    friend std::ostream& operator<<( std::ostream& os, const Scheduling_Type sstt )
    {
        if (sstt == ANALYTICAL) {
            os << "ANALYTICAL";
        }
        else if (sstt == ROUNDROBIN) {
            os << "ROUNDROBIN";
        }
        else if (sstt == STIMULUS) {
            os << "STIMULUS";
        }
        else {
            os << "UNDEFINED_SCHEDULER";
        }

        return os;
    }

    typedef enum {
        MIN_VAL,
        MAX_VAL,
        AVG_VAL
    } RechargeLength_Type;

    typedef enum {
        ENERGYMIN,
        ENERGYMAX,
        ENERGYAVG
    } DischargeProbEnergyToUse_Type;

    typedef enum {
        LOCAL_KNOWLEDGE,
        GLOBAL_KNOWLEDGE
    } GameTheoryKnowledge_Type;

    friend std::ostream& operator<<( std::ostream& os, const RechargeLength_Type rrtt )
    {
        if (rrtt == MIN_VAL) {
            os << "MINIMUM ALGO";
        }
        else if (rrtt == MAX_VAL) {
            os << "MAXIMUM ALGO";
        }
        else if (rrtt == AVG_VAL) {
            os << "AVERAGE ALGO";
        }
        else {
            os << "UNDEFINED ALGO";
        }

        return os;
    }

    typedef struct {
        L3Address addr;
        int appAddr;
        double rcvPow;
        double rcvSnr;
        Coord pos;

        simtime_t timestamp;

        double batteryLevelAbs;
        double batteryLevelPerc;
        double coveragePercentage;
        double leftLifetime;
        int nodeDegree;
        double inRechargeT;
        double gameTheoryC;
    } nodeInfo_t;

    typedef struct {
        int addr;
        int assignedRecharge;
        int executedRecharge;
        double energy;
        bool isCharging;
    } nodeAlgo_t;

    typedef struct {
        int chargingAppAddr;
        int swapNumber;
        std::list<nodeAlgo_t> nodeList;
    } groupInfo_t;

protected:
  virtual int numInitStages() const override { return NUM_INIT_STAGES; }
  virtual void initialize(int stage) override;
  virtual void handleMessageWhenUp(cMessage *msg) override;
  virtual void finish(void) override;

  virtual void sendPacket();
  //virtual void processStart();
  virtual void processPacket(cPacket *msg);

  virtual void make1secStats(void);
  virtual void make5secStats(void);

  virtual void updateNeighbourhood(void);

  virtual double calculateInterDistance(double radious);
  virtual void updateVirtualForces(void);

  double calculateRechargeStimuliEnergyFactor(void);
  double calculateRechargeStimuliTimeFactor(void);

  double calculateRechargeStimuli(void);
  double calculateRechargeThreshold(void);

  virtual void checkRecharge(void);
  double calculateSendBackoff(void);

  virtual double calculateDischargeProb(void);
  virtual void checkDischarge(void);

  virtual double calculateSwapPenalitiesEstimationCount(double estimatedSteps);
  virtual double calculateChargeDiff (double myChoice);
  virtual double reinforceTimeVal(double val);
  virtual double calculateRechargeTime(bool log);

  virtual bool checkRechargingStationFree(void);

  virtual void checkAliveDistributed(void);

  virtual int calculateNodeDegree(void);

  void getFilteredNeigh(std::map<int, nodeInfo_t> &filteredNeigh);

  virtual void checkCentralizedRecharge(void);
  virtual void checkCentralizedRechargeGroup(groupInfo_t *actGI);
  virtual void initCentralizedRecharge(void);

  virtual void checkAliveGroup(groupInfo_t *actGI);

  virtual bool checkScheduleFeasibilityGroup(groupInfo_t *actGI);
  virtual bool decideRechargeSceduling(void);
  virtual bool decideRechargeScedulingGroup(groupInfo_t *actGI);
  virtual bool decideRechargeScedulingGroupRR(groupInfo_t *actGI);
  virtual void decideRechargeScedulingGroupLast(groupInfo_t *actGI);

  int getNodeWithMaxEnergy(groupInfo_t *gi, double &battVal);
  int getNodeWithMinEnergy(groupInfo_t *gi, double &battVal);

  void updateBatteryVals(std::list<nodeAlgo_t> *list);

  void printChargingInfo(void);
  void printChargingInfo(std::ostream &ss, const char *str);
  void printChargingInfo(const char *str);

  void printDistributedChargingInfo(std::ostream &ss, const char *str);

  void putNodeInCharging(int addr);
  void putNodeInDischarging(int addr);

  double getFullCoverage(void);
  double getMyCoverageMax(void);
  double getMyCoverageActual(void);
  void printMatrix(std::vector< std::vector<bool> > &matrix);


  virtual void sendRechargeMessage(void);

  double getGameTheoryC_Sigmoid(void);
  double getGameTheoryC_Linear1(void);
  double getGameTheoryC_Linear2(void);
  double getGameTheoryC_Linear3(void);
  double getGameTheoryC_LinearDiscount(void);
  double getGameTheoryC_SigmoidDiscount(void);
  double getGameTheoryC_LinearIncrease(void);
  double getGameTheoryC_SigmoidIncrease(void);
  double getGameTheoryC_LinearIncreaseConsistent1(void);


public:
    UDPBasicRecharge() {}
    virtual ~UDPBasicRecharge();

    double getGameTheoryC(void);

    double getTheta(void);
    double getGamma(void);
    double getAlpha(void);
    double getBeta(void);
    double getP(void);

    double getEavg(bool activeOnly, GameTheoryKnowledge_Type scope);
    double getEmax(bool activeOnly, GameTheoryKnowledge_Type scope);
    double getEmin(bool activeOnly, GameTheoryKnowledge_Type scope);


    double getRechargeProbMax(bool useDishargeProbIfTheCase);

    virtual double calculateNodeDischargeProb(void);
    virtual double calculateRechargeProb(bool useDishargeProbIfTheCase);

    //bool compare_energy (const nodeAlgo_t& first, const nodeAlgo_t& second);

private:
    L3Address myAddr;
    int myAppAddr;

    cMessage *autoMsgRecharge = nullptr;
    cMessage *autoMsgCentralizedRecharge = nullptr;
    cMessage *stat1sec = nullptr;
    cMessage *stat5sec = nullptr;
    cMessage *dischargeTimer = nullptr;
    cMessage *goToCharge = nullptr;

    //cOutVector personalUniqueCoverageVector;
    cOutVector totalCoverageVector;
    cOutVector activeNodesVector;
    cOutVector rechargingNodesVector;
    cOutVector stimulusVector;
    cOutVector thresholdVector;
    cOutVector responseVector;
    cOutVector degreeVector;
    cOutVector timeFactorVector;
    cOutVector energyFactorVector;
    cOutVector energyVector;
    cOutVector failedAttemptVector;
    cOutVector dischargeProbVector;

    VirtualSpringMobility *mob = nullptr;
    power::SimpleBattery *sb = nullptr;

    std::map<int, nodeInfo_t> neigh;
    std::list<groupInfo_t> groupList;

    std::list<int> rechargingNodesStats;

    Coord lastPosBeforeCharge;

    simtime_t lastRechargeTimestamp;

    double reinforcementVal;

    double inRechargingTime;
    simtime_t startRecharge;

    int failedAttemptCount;

    //parameters
    double checkRechargeTimer;
    double sensorRadious;
    Coord rebornPos;

    bool isCentralized;
    Scheduling_Type st;
    RechargeLength_Type rlt;
    int chargingStationNumber;
    int numRechargeSlotsStimulusZeroNeigh;
    DischargeProbEnergyToUse_Type dischargeProbEnergyToUse;
    GameTheoryKnowledge_Type gameTheoryKnowledgeType;

    int roundrobinRechargeSize;

    int numRechargeSlotsProbabilistic;

    double exponential_dischargeProb_decay;

    double stimulusExponent;

    double const_c_dischargeProb;
    bool makeLowEnergyFactorCurves;
    double timeFactorMultiplier;
    bool godCheckIfRechargeStationFree;
    bool firstRecharge;

    bool returnBackAfterRecharge;

    bool stationANDnodeKNOWN;

    char logFile[256];
    bool printAnalticalLog;

    int rechargeLostAccess;

    double dicountminLINEAR4;
    double temp_factorProbDischarge;

    bool reinforcementRechargeTime;
    double reinforcementRechargeAlpha;
    double reinforcementRechargeAlphaFinal;

    double chargeTimeOthersNodeFactor;

    bool makeCoverageLog;

    //bool developingStimuli;
    double constantTheta;
    Stimulus_Type stim_type;
    VarConstant_Type constant_type;
    VarProbability_Type probability_type;
};

} /* namespace inet */

#endif /* INET_APPLICATIONS_UDPAPP_UDPBASICRECHARGE_H_ */
