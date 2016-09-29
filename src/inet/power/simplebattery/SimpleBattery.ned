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

package inet.power.simplebattery;

//
// TODO auto-generated type
//
simple SimpleBattery
{
    double nominalCapacity @unit(J) = default(1000J);
    double initialCapacity @unit(J) = default(nominalCapacity);
    double updateInterval @unit(s) = default(0.2s);
    
    double chargingFactor = default(2);
    double dischargingFactor = default(4);
    
    @class(SimpleBattery);
}
