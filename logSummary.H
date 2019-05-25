{
    const dictionary& subModelDict = parcels.subModelProperties();
    vector position = subModelDict.subDict
    (
        "injectionModels"
    ).subDict("model1").lookup("position");

    vector direction = subModelDict.subDict
    (
        "injectionModels"
    ).subDict("model1").lookup("direction");
    direction /= mag(direction);

    scalar vaporPenetration = 0.;

    const wordList fuelList(thermo.subDict("liquids").keys());

    scalarField fuelMass(T.size(), 0.);
    forAll(Y, i)
    {
        forAll(fuelList, fuelID)
        {
            if ( Y[i].name() == fuelList[fuelID])
            {
                //Info<<"Y[i].name()==="<<Y[i].name()<<endl;
                fuelMass += Y[i].primitiveField();
            }
        }
    }

    forAll (fuelMass,cellI)
    {
        if (fuelMass[cellI] >= 0.001)
        {
            vector raw = position - mesh.C()[cellI];
            if (mag(raw&direction) > vaporPenetration)
            {
                vaporPenetration = mag(raw&direction);
            }
        }
    }

    reduce(vaporPenetration, maxOp<scalar>());
    
    label OH_index = 0;
    if (composition.species().contains("OH")) 
    {
        OH_index = composition.species()["OH"];
    }
    scalar Y_OH_max = gMax(Y[OH_index]);
    scalar liftOffLength(0);

    forAll (Y[0],cellI)
    {
        if (Y[OH_index][cellI] >= 0.02*Y_OH_max)
        {
            vector raw = position - mesh.C()[cellI];
            liftOffLength = mag(raw&direction);
            if (mag(raw&direction) < liftOffLength)
            {
                liftOffLength = mag(raw&direction);
            }          
        }
    }

    reduce(liftOffLength, maxOp<scalar>());    


    scalar T_mean = T.weightedAverage(mesh.V()).value();
    scalar T_max = gMax(T);
    scalar T_min = gMin(T);
    scalar T_max_old = gMax(T.oldTime());
    scalar dTdt = (T_max-T_max_old)/runTime.deltaTValue();
    scalar D10 = 1.0e+6*parcels.Dij(1, 0);
    scalar D32 = 1.0e+6*parcels.Dij(3, 2);
    scalar Dmax = 1.0e+6*parcels.Dmax();
    scalar pene95 = 1.0e+3*parcels.penetration(0.95);
    scalar pene97 = 1.0e+3*parcels.penetration(0.97);
    scalar pene98 = 1.0e+3*parcels.penetration(0.98);
    scalar pene99 = 1.0e+3*parcels.penetration(0.99);
    scalar pene999 = 1.0e+3*parcels.penetration(0.999);
    
    scalar HeatRelease = fvc::domainIntegrate(combustion->Qdot()).value(); // J/s    Qdot--dimEnergy/dimVolume/dimTime
    HRtotal += HeatRelease*runTime.deltaTValue(); // J    
    
    logSummaryFile()
        << 1000*runTime.value() << tab //ms
        << p.weightedAverage(mesh.V()).value()/1.e6 << tab //MPa
        << p.weightedAverage(mesh.V()).value()/1.e6 - pInitial << tab //MPa //pressure rise        
        << HeatRelease/1000. << tab //J/ms
        << HRtotal << tab //J       
        << T_mean << tab //K
        << T_max<< tab //K
        << T_min<< tab //K
        << D10 << tab //mu
        << D32 << tab //mu
        << Dmax << tab //mu
        << pene95 << tab //mm
        << pene97 << tab //mm
        << pene98 << tab //mm
        << pene99 << tab //mm
        << pene999 << tab //mm
		<< 1.0e+3*vaporPenetration << tab //mm
		<< Y_OH_max << tab //-
		<< 0.02*Y_OH_max << tab //-
		<< 1.0e+3*liftOffLength << tab //mm
        << dTdt << tab
        << endl;
}
