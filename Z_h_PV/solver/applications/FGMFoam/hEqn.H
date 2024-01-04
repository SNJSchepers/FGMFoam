// Solve the enthalpy equation
fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + mvConvection->fvmDiv(phi, h)
  - fvc::div( thermo.alpha()       * fvc::grad(h)  
            + thermo.prefDh_H()   * fvc::grad(thermo.H())
            + thermo.prefDh_H2()  * fvc::grad(thermo.H2())
            + thermo.prefDh_H2O() * fvc::grad(thermo.H2O()) )
 ==
    fvOptions(rho, h)
);

hEqn.relax();
fvOptions.constrain(hEqn);
hEqn.solve();

