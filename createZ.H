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