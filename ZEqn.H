{
    fvScalarMatrix ZEqn
    (
        fvm::ddt(rho, Z)
      + mvConvection->fvmDiv(phi, Z)
      - fvm::laplacian(turbulence->muEff(), Z)
        ==
        parcels.Srho(Z)
    );

    ZEqn.relax();
    ZEqn.solve("Yi");
    Z.max(0.0);

    Chi = 2*turbulence->epsilon()/turbulence->k()*varZ;

    if (mesh.objectRegistry::foundObject<compressible::LESModel>(turbulenceModel::propertiesName))
    {
        const compressible::LESModel& lesModel
            = mesh.objectRegistry::lookupObject<compressible::LESModel>(turbulenceModel::propertiesName);
        Chi = 2*turbulence->nu()*magSqr(fvc::grad(Z)) + 2*turbulence->nut()*varZ/Foam::sqr(lesModel.delta());
        varZ_algebraic = 2*Foam::sqr(lesModel.delta())*magSqr(fvc::grad(Z));
    }

    fvScalarMatrix varZEqn
    (
        fvm::ddt(rho, varZ)
      + mvConvection->fvmDiv(phi, varZ)
      - fvm::laplacian(turbulence->muEff(), varZ)
        ==
      + 2*turbulence->muEff()*magSqr(fvc::grad(Z))
      - rho*Chi
    );

    varZEqn.relax();
    varZEqn.solve("Yi");
    varZ.max(0.0);
}