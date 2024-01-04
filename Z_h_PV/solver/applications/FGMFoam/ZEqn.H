// Convection scheme
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Z)")
    )
);

// Solve the enthalpy equation
fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Z)
  + mvConvection->fvmDiv(phi, Z)
  - fvc::div( thermo.alpha()       * fvc::grad(Z)  
            + thermo.prefDZ_H()   * fvc::grad(thermo.H())
            + thermo.prefDZ_H2()  * fvc::grad(thermo.H2())
            + thermo.prefDZ_H2O() * fvc::grad(thermo.H2O()) )
 ==
    fvOptions(rho, Z)
);

ZEqn.relax();
fvOptions.constrain(ZEqn);
ZEqn.solve();

