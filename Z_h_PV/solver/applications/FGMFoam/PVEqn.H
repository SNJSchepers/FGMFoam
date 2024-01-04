// Solve the progress variable equation
fvScalarMatrix PVEqn
(
    fvm::ddt(rho, PV)
  + mvConvection->fvmDiv(phi, PV)
  - fvc::div( thermo.alpha()       * fvc::grad(PV)  
            + thermo.prefDPV_H()   * fvc::grad(thermo.H())
            + thermo.prefDPV_H2()  * fvc::grad(thermo.H2())
            + thermo.prefDPV_H2O() * fvc::grad(thermo.H2O()) )
 ==
    thermo.sourcePV()
  + fvOptions(rho, PV)
);

PVEqn.relax();
fvOptions.constrain(PVEqn);
PVEqn.solve();

