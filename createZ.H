// use this constructor:
// GeometricField(const IOobject&, const Internal&, const PtrList<PatchField<Type>>&);

volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0*Y[0],
    Y[0].boundaryField()
);

volScalarField varZ
(
    IOobject
    (
        "varZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0*Y[0],
    Y[0].boundaryField()
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
    dimensionedScalar(dimless/dimTime, 0)
    // default calculatedType
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
    dimensionedScalar(dimless, 0)
    // default calculatedType
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