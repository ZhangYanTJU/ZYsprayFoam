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


    scalar vaporPenetration(0);

    const word fuel(thermo.lookup("fuel"));
    Info<<"fuel==="<<fuel;
    const label fuelIndex(composition.species()[fuel]);
    Info<<"=====fuelIndex==="<<fuelIndex;
    forAll (Y[fuelIndex],cellI)
    {
        if (Y[fuelIndex][cellI] >= 0.001)
        {
            vector raw = position - mesh.C()[cellI];
            if (mag(raw&direction) > vaporPenetration)
            {
                vaporPenetration = mag(raw&direction);
            }    
        }
    }

    reduce(vaporPenetration, maxOp<scalar>());


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
  
    logSummaryFile()
        << 1000*runTime.value() << tab
        << T_mean << tab
        << T_max<< tab
        << T_min<< tab
        << D10 << tab
        << D32 << tab
        << Dmax << tab
        << pene95 << tab
        << pene97 << tab
        << pene98 << tab
        << pene99 << tab	
        << pene999 << tab
        << 1.0e+3*vaporPenetration << tab
        << dTdt << tab
        << endl;
}
