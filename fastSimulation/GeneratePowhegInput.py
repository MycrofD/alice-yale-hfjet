#!/usr/bin/env python

import shutil
import argparse
import random
import yaml


def GetParallelInputFileName(powheg_stage, x_grid_iter=1):
    if powheg_stage == 1:
        fname = "powheg_Stage_{}_XGrid_{}.input".format(powheg_stage, x_grid_iter)
    else:
        fname = "powheg_Stage_{}.input".format(powheg_stage)
    return fname


def GenerateParallelPowhegInput(outputdir, powheg_stage, x_grid_iter, powhegEvents, gen, powheg_proc, qmass, facscfact, renscfact, lhans, beamType, ebeam1, ebeam2, bornktmin, nPDFset, nPDFerrSet):
    fname = "{}/{}".format(outputdir, GetParallelInputFileName(powheg_stage, x_grid_iter))
    shutil.copy("{}-powheg.input".format(powheg_proc), fname)

    with open(fname, "a") as myfile:
        myfile.write("numevts {0}\n".format(powhegEvents))
        myfile.write("manyseeds 1\n")
        myfile.write("maxseeds 500\n")
        myfile.write("parallelstage {}\n".format(powheg_stage))
        if powheg_proc == "beauty":
            myfile.write("qmass {0}\n".format(qmass))
            myfile.write("facscfact {0}\n".format(facscfact))
            myfile.write("renscfact {0}\n".format(renscfact))
            myfile.write("ncall1 2000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 5000\n")
            myfile.write("itmx2 5\n")
        elif powheg_proc == "charm":
            myfile.write("qmass {0}\n".format(qmass))
            myfile.write("facscfact {0}\n".format(facscfact))
            myfile.write("renscfact {0}\n".format(renscfact))
            myfile.write("ncall1 10000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 5000\n")
            myfile.write("itmx2 5\n")
        elif powheg_proc == "dijet":
            myfile.write("bornktmin {0}\n".format(bornktmin))
            myfile.write("ncall1 5000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 1000\n")
            myfile.write("itmx2 5\n")

        if powheg_stage == 1: myfile.write("xgriditeration {}\n".format(x_grid_iter))
        myfile.write("lhans1 {0}\n".format(lhans))
        myfile.write("lhans2 {0}\n".format(lhans))
        myfile.write("ebeam1 {0}\n".format(ebeam1))
        myfile.write("ebeam2 {0}\n".format(ebeam2))

        if beamType == "pPb":
            myfile.write("nPDFset {0}        ! (0:EKS98, 1:EPS08, 2:EPS09LO, 3:EPS09NLO)\n".format(nPDFset))
            myfile.write("nPDFerrSet {0}     ! (1:central, 2:+1, 3:-1..., 30:+15, 31:-15)\n".format(nPDFerrSet))
            myfile.write("AA1 208            ! (Atomic number of hadron 1)\n")
            myfile.write("AA2 1              ! (Atomic number of hadron 2)\n")


def GenerateSinglePowhegInput(outputdir, powhegEvents, gen, powheg_proc, qmass, facscfact, renscfact, lhans, beamType, ebeam1, ebeam2, bornktmin, nPDFset, nPDFerrSet):
    fname = "{}/powheg.input".format(outputdir)
    shutil.copy("{}-powheg.input".format(powheg_proc), fname)

    with open(fname, "a") as myfile:
        myfile.write("numevts {0}\n".format(powhegEvents))
        if powheg_proc == "beauty":
            myfile.write("qmass {0}\n".format(qmass))
            myfile.write("facscfact {0}\n".format(facscfact))
            myfile.write("renscfact {0}\n".format(renscfact))
            myfile.write("ncall1 10000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 100000\n")
            myfile.write("itmx2 5\n")
        elif powheg_proc == "charm":
            myfile.write("qmass {0}\n".format(qmass))
            myfile.write("facscfact {0}\n".format(facscfact))
            myfile.write("renscfact {0}\n".format(renscfact))
            myfile.write("ncall1 50000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 100000\n")
            myfile.write("itmx2 5\n")
        elif powheg_proc == "dijet":
            myfile.write("bornktmin {0}\n".format(bornktmin))
            myfile.write("ncall1 20000\n")
            myfile.write("itmx1 5\n")
            myfile.write("ncall2 20000\n")
            myfile.write("itmx2 5\n")
        myfile.write("lhans1 {0}\n".format(lhans))
        myfile.write("lhans2 {0}\n".format(lhans))
        myfile.write("ebeam1 {0}\n".format(ebeam1))
        myfile.write("ebeam2 {0}\n".format(ebeam2))
        if beamType == "pPb":
            myfile.write("nPDFset {0}        ! (0:EKS98, 1:EPS08, 2:EPS09LO, 3:EPS09NLO)\n".format(nPDFset))
            myfile.write("nPDFerrSet {0}     ! (1:central, 2:+1, 3:-1..., 30:+15, 31:-15)\n".format(nPDFerrSet))
            myfile.write("AA1 208            ! (Atomic number of hadron 1)\n")
            myfile.write("AA2 1              ! (Atomic number of hadron 2)\n")


def main(yamlConfigFile, outputdir, powhegEvents, powheg_stage, x_grid_iter=1):
    f = open(yamlConfigFile, 'r')
    config = yaml.load(f)
    f.close()

    gen = config["gen"]
    proc = config["proc"]
    if "qmass" in config: qmass = config["qmass"]
    else: qmass = None
    if "facscfact" in config: facscfact = config["facscfact"]
    else: facscfact = None
    if "renscfact" in config: renscfact = config["renscfact"]
    else: renscfact = None
    lhans = config["lhans"]
    beamType = config["beam_type"]
    ebeam1 = config["ebeam1"]
    ebeam2 = config["ebeam2"]
    if "bornktmin" in config:
        bornktmin = config["bornktmin"]
    else:
        bornktmin = 10
    nPDFset = config["nPDFset"]
    nPDFerrSet = config["nPDFerrSet"]

    if proc == "charm_jets" or proc == "beauty_jets":
        powheg_proc = "dijet"
        powhegEvents *= 5
    else:
        powheg_proc = proc

    if qmass < 0:
        if powheg_proc == "charm":
            qmass = 1.5
        elif powheg_proc == "beauty":
            qmass = 4.75

    shutil.copy("{}-powheg.input".format(powheg_proc), "{}/powheg.input".format(outputdir))

    if powheg_stage > 0 and powheg_stage <= 4:
        GenerateParallelPowhegInput(outputdir, powheg_stage, x_grid_iter, powhegEvents, gen, powheg_proc, qmass, facscfact, renscfact, lhans, beamType, ebeam1, ebeam2, bornktmin, nPDFset, nPDFerrSet)
    else:
        GenerateSinglePowhegInput(outputdir, powhegEvents, gen, powheg_proc, qmass, facscfact, renscfact, lhans, beamType, ebeam1, ebeam2, bornktmin, nPDFset, nPDFerrSet)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate POWHEG input file.')
    parser.add_argument('config', metavar='config.yaml',
                        default="default.yaml", help='YAML configuration file')
    parser.add_argument('--numevents', metavar='NEVT',
                        default=50000, type=int)
    parser.add_argument('-o', metavar='job',
                        default='./')
    parser.add_argument('--powheg-stage', metavar='STAGE',
                        default=0, type=int)
    args = parser.parse_args()

    main(args.config, args.o, args.numevents, args.powheg_stage)
