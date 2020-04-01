OFstream logSummaryFile
(
    //runTime.path()/("logSummary." + runTime.timeName() + ".dat")
    runTime.time().rootPath()/runTime.time().globalCaseName()/("logSummary." + runTime.timeName() + ".dat")
);

scalar pInitial = p.weightedAverage(mesh.V()).value()/1.e6;
scalar HRtotal = 0.;
scalar Y_OH_max_max = small;

logSummaryFile
    << "Time (ms)" << tab
    << "p (MPa)" << tab
    << "pressure rise (MPa)" << tab
    << "HeatRelease (J/ms)" << tab
    << "total heat release (J)" << tab
    << "Tmean (K)" << tab
    << "Tmax (K)" << tab
    << "Tmin (K)" << tab
    << "D10 (mu)" << tab
    << "SMD(D32) (mu)" << tab
    << "Dmax (mu)" << tab
    << "L95 (mm)" << tab
    << "L97 (mm)" << tab
    << "L98 (mm)" << tab
    << "L99 (mm)" << tab
    << "L999 (mm)" << tab
    << "Vapor (mm)" << tab
    << "Y_OH_max (-)" << tab
    << "2%_Y_OH_max (-)" << tab
    << "liftOffLength (mm)" << tab
    << "parcelsMassInSystem (kg)" << tab
    << "dTdt (K/s)" << endl;
