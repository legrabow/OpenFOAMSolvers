while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix CEqn
        (
            fvm::ddt(rho, C)
            + fvm::div(phi, C)
             - fvm::laplacian(gammaC, C)
            ==
            fvOptions(rho, C)
        );

        CEqn.relax();
        fvOptions.constrain(CEqn);
        CEqn.solve();
        fvOptions.correct(C);
    }
