volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField varZ
(
    IOobject
    (
        "varZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Chi
(
    IOobject
    (
        "Chi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField varZ_algebraic
(
    IOobject
    (
        "varZ_algebraic",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField hc
(
    IOobject
    (
        "hc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.hc()
);

volScalarField ha
(
    IOobject
    (
        "ha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.he()+thermo.hc()
);