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
/*
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
*/
forAll(fuelList, fuelI)
{
    fuelMass += Y[composition.species()[fuelList[fuelI]]].primitiveField();
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
if (composition.species().found("OH"))
{
    OH_index = composition.species()["OH"];
}
scalar Y_OH_max = gMax(Y[OH_index]);
scalar liftOffLength = great;

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

reduce(liftOffLength, minOp<scalar>());


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
        << scientific << 0.02*Y_OH_max << tab //-
        << scientific << 1.0e+3*liftOffLength << tab //mm
        << scientific << dTdt << endl;
}
