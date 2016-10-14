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

package inet.applications.udpapp;


import inet.applications.udpapp.UDPBasicApp;

//
// TODO auto-generated type
//
simple UDPBasicRecharge extends UDPBasicApp
{
    parameters:
        double checkRechargeTimer @unit(s) = default(30s);
        double sensorRadious @unit(m) = default(50m);
        string schedulingType @enum("ANALYTICAL", "ROUNDROBIN", "STIMULUS", "PROBABILISTIC") = default("STIMULUS");
        //bool isCentralized = default(false);
        
        int chargingStationNumber = default(3);
        
        int roundrobinRechargeSize = default(1);
        
        double stimulusExponent = default(1);
        
        int numRechargeSlotsProbabilistic = default(10);
        int numRechargeSlotsStimulusZeroNeigh = default(3);
        
        bool makeLowEnergyFactorCurves = default(true);
        double timeFactorMultiplier = default(2);
        
        bool godCheckIfRechargeStationFree = default(true);
        
        string analticalLogFile = default("/dev/null");
        bool printAnalticalLog = default(false);
        
        bool returnBackAfterRecharge = default(true);
        
        @display("i=block/app");
        @class(UDPBasicRecharge);
}
