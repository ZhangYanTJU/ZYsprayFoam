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

scalarField Y_fuel(T.size(), 0.);
forAll(fuelList, fuelI)
{
    Y_fuel += Y[composition.species()[fuelList[fuelI]]].primitiveField();
}

forAll (Y_fuel,cellI)
{
    if (Y_fuel[cellI] >= 0.001)
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
if (composition.species().found("OH"))
{
    OH_index = composition.species()["OH"];
}

scalar Y_OH_max = gMax(Y[OH_index]);

if (Y_OH_max > Y_OH_max_max) Y_OH_max_max = Y_OH_max;
Info << "max(Y_OH) until now === " << Y_OH_max_max << endl;

scalar liftOffLength = great;

forAll (Y[0],cellI)
{
    if (Y[OH_index][cellI] >= 0.02*Y_OH_max)
    {
        vector raw = position - mesh.C()[cellI];
        if (mag(raw&direction) < liftOffLength)
        {
            liftOffLength = mag(raw&direction);
        }
    }
}
reduce(liftOffLength, minOp<scalar>());


scalar liftOffLength14 = great;

forAll (Y[0],cellI)
{
    if (Y[OH_index][cellI] >= 0.14*Y_OH_max)
    {
        vector raw = position - mesh.C()[cellI];
        if (mag(raw&direction) < liftOffLength14)
        {
            liftOffLength14 = mag(raw&direction);
        }
    }
}
reduce(liftOffLength14, minOp<scalar>());


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
scalar p_weightedAverage = p.weightedAverage(mesh.V()).value();

scalar parcelsMassInSystem = parcels.massInSystem();
reduce(parcelsMassInSystem, sumOp<scalar>());

if (Pstream::master())
{
    logSummaryFile()
        << scientific << 1000*runTime.value() << tab //ms
        << scientific << p_weightedAverage/1.e6 << tab //MPa
        << scientific << p_weightedAverage/1.e6 - pInitial << tab //MPa //pressure rise
        << scientific << HeatRelease/1000. << tab //J/ms
        << scientific << HRtotal << tab //J
        << scientific << T_mean << tab //K
        << scientific << T_max<< tab //K
        << scientific << T_min<< tab //K
        << scientific << D10 << tab //mu
        << scientific << D32 << tab //mu
        << scientific << Dmax << tab //mu
        << scientific << pene95 << tab //mm
        << scientific << pene97 << tab //mm
        << scientific << pene98 << tab //mm
        << scientific << pene99 << tab //mm
        << scientific << pene999 << tab //mm
        << scientific << 1.0e+3*vaporPenetration << tab //mm
        << scientific << Y_OH_max << tab //-
        << scientific << 1.0e+3*liftOffLength << tab //mm
        << scientific << 1.0e+3*liftOffLength14 << tab //mm
        << scientific << parcelsMassInSystem << tab //kg
        << scientific << dTdt << endl;
}


// output for parcels if needed, you can uncomment this block to see more information
// parcels (cloud) > parcel > particle
/*
Pout << "total nParcels === " << parcels.nParcels() << endl;
Info << "averageParcelMass (sometimes not true when massTotal in SprayCloudProperties is not true) === " << parcels.averageParcelMass() << endl;
Pout << "total parcel mass (sometimes not true, neither)=== " << parcels.nParcels()*parcels.averageParcelMass() << endl;
Pout << "parcelsMassInSystem === " << parcels.massInSystem() << endl;

//Info << "InjectionModelList-averageParcelMass === " << parcels.injectors().averageParcelMass() << endl; // may cause error
Info << "InjectionModelList-size === " << parcels.injectors().size() << endl;


forAll(parcels.injectors(), injectorI)
{
    Info << "InjectionModel-nInjections (cumulative) === " << parcels.injectors()[injectorI].nInjections() << endl;
    Info << "InjectionModel-parcelsAddedTotal (cumulative) === " << parcels.injectors()[injectorI].parcelsAddedTotal() << endl;
    //Info << "InjectionModel-averageParcelMass === " << parcels.injectors()[injectorI].averageParcelMass() << endl; // may cause error
    Info << "InjectionModel-massTotal (massTotal you specified in SprayCloudProperties) === " << parcels.injectors()[injectorI].massTotal() << endl;
    Info << "InjectionModel-mass injected (cumulative) === " << parcels.injectors()[injectorI].massInjected() << endl;
    Info << "InjectionModel-volumeTotal_ ( flowRateProfile.integrate(0, duration) ) === " << parcels.injectors()[injectorI].volumeTotal() << endl;
    Info << "InjectionModel-Volume of parcels to inject = volumeToInject(deltaT)/volumeTotal_ )" << endl;
}
*/


// if needed, you can uncomment this block to see more information
/*
forAllConstIter(basicSprayCloud, parcels, pIter)
{
    const basicReactingParcel& p = pIter();
    const label celli = p.cell();

    //sumParticles[celli] += p.nParticle();
    Pout << "in this parcel, nParticle === " << p.nParticle() << endl;
    Pout << "parcel position (x y z)=== " << p.position() << endl;
    Pout << "particle rho === " << p.rho() << endl;
    Pout << "particle T === " << p.T() << endl;
    Pout << "particle d === " << p.d() << endl;
    Pout << "particle volume === " << p.volume() << endl;
    Pout << "particle volume-ZY === " << Foam::constant::mathematical::pi/6.0*pow3(p.d()) << endl;
    Pout << "particle mass === " << p.mass() << endl;
    Pout << "particle mass-ZY === " << p.rho()*Foam::constant::mathematical::pi*pow3(p.d())/6. << endl;
}
*/
