OFstream logSummaryFile
(
    //runTime.path()/("logSummary." + runTime.timeName() + ".dat")
    runTime.time().rootPath()/runTime.time().globalCaseName()/("logSummary." + runTime.timeName() + ".dat")
);

scalar pInitial = p.weightedAverage(mesh.V()).value()/1.e6;
scalar HRtotal = 0.;

logSummaryFile
    << "Time" << tab
    << "p" << tab
    << "pressure rise" << tab
    << "HeatRelease" << tab //J/ms
    << "total heat release" << tab // J
    << "Tmean" << tab
    << "Tmax" << tab
    << "Tmin" << tab
    << "D10" << tab
    << "SMD(D32)" << tab
    << "Dmax" << tab
    << "L95" << tab
    << "L97" << tab
    << "L98" << tab
    << "L99" << tab
    << "L999" << tab
    << "Vapor" << tab
    << "Y_OH_max" << tab
    << "2%_Y_OH_max" << tab
    << "liftOffLength" << tab
    << "dTdt" << tab
    << endl;

logSummaryFile
    << "ms" << tab
    << "MPa" << tab
    << "MPa" << tab
    << "J/ms" << tab
    << "J" << tab
    << "K" << tab
    << "K" << tab
    << "K" << tab
    << "mu" << tab
    << "mu" << tab
    << "mu" << tab
    << "mm" << tab
    << "mm" << tab
    << "mm" << tab
    << "mm" << tab
    << "mm" << tab
    << "mm" << tab
    << "-" << tab
    << "-" << tab
    << "mm" << tab
    << "K/s" << endl;
