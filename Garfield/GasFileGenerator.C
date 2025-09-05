#include <iostream>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  const double pressure = 3 * AtmosphericPressure;
  const double temperature = 293.15;
 
  // Setup the gas.
  MediumMagboltz gas("Ar", 93., "CO2", 7.);
  gas.SetTemperature(temperature);
  gas.SetPressure(pressure);

  // Set the field range to be covered by the gas table. 
  const size_t nE = 20;
  const double emin =    100.;
  const double emax = 100000.;
  // Flag to request logarithmic spacing.
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog); 

  const int ncoll = 10;
  // Run Magboltz to generate the gas table.
  gas.GenerateGasTable(ncoll);
  // Save the table. 
  gas.WriteGasFile("ar_93_co2_7.gas");

}