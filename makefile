EXEC = show
CC = mpicc
OBJS = main.o parameterSetting.o findparam.o boundary.o saveFile.o loadLaser.o fieldSolve.o loadPlasma_crystal.o loadPlasma.o rearrangeParticles.o removeEdge.o particlePush.o updateCurrent.o interpolation.o #fieldShareY.o fieldShareZ.o particleShareY.o particleShareZ.o movingDomain.o loadMovingPlasma.o clean.o pml.o dumpData.o densityShare.o track.o saveFieldHDF.o saveDensityHDF.o saveRamanHDF.o saveDumpHDF.o restoreDumpHDF.o saveParticleHDF.o

INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm -lhdf5
$(EXEC):$(OBJS)
#	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
