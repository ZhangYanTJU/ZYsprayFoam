{
	fvScalarMatrix ZEqn
	(
		fvm::ddt(rho, Z)
      + mvConvection->fvmDiv(phi, Z)
      - fvm::laplacian(turbulence->alphaEff(), Z)
    	==
	    parcels.Srho(Z)
	);

    ZEqn.relax();

    ZEqn.solve("Yi");

    Z.max(0.0);
}