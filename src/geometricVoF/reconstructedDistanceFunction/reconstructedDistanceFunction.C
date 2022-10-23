/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "emptyPolyPatch.H"
#include "reconstructedDistanceFunction.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "wedgePolyPatch.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::reconstructedDistanceFunction::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches) // loop patch in all patches
    {
        if (isA<coupledPolyPatch>(pp)) // if the patch is coupled with another patch
        {
            nCoupled += pp.size(); // sum the number of faces in all coupled patches
        }
    }
    labelList nCoupledFaces(nCoupled); // a empty label list with the list-size of nCoupled
    nCoupled = 0; // reassigned the nCoupled

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            label facei = pp.start(); // start face index of the patch pp

            forAll(pp, i) // loop all faces in the patch pp
            {
                nCoupledFaces[nCoupled++] = facei++; // fill all faces-indexes(absolute) of all coupled patches in the nCoupledFaces list. 
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New  // new autoPtr to indirect addressing of points
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nCoupledFaces
        ),
        mesh_.points()
    );
}


void Foam::reconstructedDistanceFunction::markCellsNearSurf
(
    const boolList& interfaceCells,
    const label neiRingLevel
)
{
    // performance might be improved by increasing the saving last iterations
    // cells in a Map and loop over the map
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (nextToInterface_.size() != mesh_.nCells())
        {
            nextToInterface_.resize(mesh_.nCells());
        }
        coupledBoundaryPoints_ = coupledFacesPatch()().meshPoints();
    }

    const labelListList& pCells = mesh_.cellPoints(); // cell to point, the points of a cell?? correct
    const labelListList& cPoints = mesh_.pointCells(); // point to cell, point shared cells?? correct

    boolList alreadyMarkedPoint(mesh_.nPoints(), false);
    nextToInterface_ = false;

    // do coupled face first
    Map<bool> syncMap;

    for (int level=0;level<=neiRingLevel;level++)
    {
        // parallel
        if (level > 0)
        {
            forAll(coupledBoundaryPoints_, i)  // all points index on coupledBoundaryPatches.
            {
                const label pi = coupledBoundaryPoints_[i]; // global point index 
                forAll(mesh_.pointCells()[pi], j)
                {
                    const label celli = cPoints[pi][j];
                    if (cellDistLevel_[celli] == level-1) // at level = 1, cellDistLevel_: 0 --> interfacecells; -1 --> all other 
                    {                                     // at level = 1, level -1 = 0, if (isAInterfaceCell)
                        syncMap.insert(pi, true); // map : key: a point on coupled boundary patch, value: true / false
                        break;                    // at level = 1: insert map (coupledPatchPoint, true), 
                    }                             // in which coupledPatchPoint has at least one point interfaceCell
                }
            }

            syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>()); //  Synchronize values on selected points. 
                                     //--> cyclic patches considered, when one of the coupled points is inserted, then the another is also inserted (with value true).

            // mark parallel points first
            forAllConstIters(syncMap, iter)
            {
                const label pi = iter.key();

                if (!alreadyMarkedPoint[pi]) // initialised to all false
                {
                    // loop over all cells attached to the point
                    forAll(cPoints[pi], j)
                    {
                        const label pCelli = cPoints[pi][j];
                        if (cellDistLevel_[pCelli] == -1) // at level = 1, cellDistLevel_: 0-> interfaceCell; 
                        {                                 // 1 --> cell contains no interface, but its vertex-shared neigbor cell contain interface
                            cellDistLevel_[pCelli] = level; // -1 --> all oter non-interface cells
                            nextToInterface_[pCelli] = true; // the cell with cellDistLevel_ = 1
                        }
                    }
                }
                alreadyMarkedPoint[pi] = true;
            }
        }


        forAll(cellDistLevel_, celli)  // cellDistLevel_ --> volScalarField, initialised with -1
        {
            if (level == 0) // for level = 0, all interfaceCells are marked: 
            {               // 1. cellDistLevel_ (initial value -1) is set to 0; 
                if (interfaceCells[celli])  // 2. nextToInterface_ (initial value false) is set to true
                {
                    cellDistLevel_[celli] = 0;
                    nextToInterface_[celli] = true; 
                }
                else // non interfacecells --> unchanged
                {
                    cellDistLevel_[celli] = -1;
                }
            }
            else
            {
                if (cellDistLevel_[celli] == level-1) // at level = 1, this means if(celli isAInterfaceCell)
                {
                    forAll(pCells[celli], i) // loop all labels of vertices belong to interfaceCell celli 
                    {
                        const label pI = pCells[celli][i]; // global index of vertex belonging to celli

                        if (!alreadyMarkedPoint[pI])
                        {
                            forAll(cPoints[pI], j) // loop all cells sharing the vertex
                            {
                                const label pCelli = cPoints[pI][j];
                                if (cellDistLevel_[pCelli] == -1) // if (the vertex-shared cell is non-interface cell)
                                {
                                    cellDistLevel_[pCelli] = level; // the same: 0 interfaceCell; 1 vertex-shared cell; -1 oter
                                    nextToInterface_[pCelli] = true; // cellDistLevel_ 1
                                }
                            }
                        }
                        alreadyMarkedPoint[pI] = true;
                    }
                }
            }
        }
    }
} // after this func: 1. nextToInterface_: a) true for all interfaceCells and all non-interfaceCells sharing at least one vertex with interface cells
  //                                       b) false for all other non-interfaceCell
  //                  2. cellDistLevel_ : -1, 0, 1; when 0,1 --> nextToInterface_ == true
  // Note: coupled BC (cyclic, processor) are considered.


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructedDistanceFunction::reconstructedDistanceFunction
(
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            "RDF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // --> the reason, the RDF shows in every non-zero time folder
        ),
        mesh,
        dimensionedScalar(dimLength, Zero)
    ),
    mesh_(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints()), // label list, size = ....meshPoints().size(), value = 0
    cellDistLevel_
    (
        IOobject
        (
            "cellDistLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("cellDistLevel", dimless, -1) // initialised with value -1
    ),
    nextToInterface_(mesh.nCells(), false) // boolean list, value false, size = the number of cells--> if cell is nextToInterface
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField&  Foam::reconstructedDistanceFunction::constructRDF
(
    const boolList& nextToInterface,
    const volVectorField& centre,
    const volVectorField& normal,
    zoneDistribute& distribute, // !!!!!
    bool updateStencil
)
{
    volScalarField& reconDistFunc = *this;

    if (nextToInterface.size() != centre.size())
    {
        FatalErrorInFunction
            << "size of nextToInterface: " << nextToInterface.size()
            << "size of centre:" <<  centre.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
        return reconDistFunc;
    }


    distribute.setUpCommforZone(nextToInterface, updateStencil); // updatestencil --> boolean value, if stencil is updated

    Map<vector> mapCentres =
        distribute.getDatafromOtherProc(nextToInterface, centre);   // !!!  only considering the data from other proc
    Map<vector> mapNormal =
        distribute.getDatafromOtherProc(nextToInterface, normal);   // !! we need also data from other paired cyclic patch internal cell

    const labelListList& stencil = distribute.getStencil();  // cell point cell stencil


    forAll(nextToInterface,celli)
    {
        if (nextToInterface[celli])
        {
            if (mag(normal[celli]) != 0) // interface cell--> level = 0
            {
                vector n = -normal[celli]/mag(normal[celli]);
                scalar dist = (centre[celli] - mesh_.C()[celli]) & n; // distance between cell center and interface center
                reconDistFunc[celli] = dist;
            }
            else // level == 1 and interface cell with normal == 0
            {
                scalar averageDist = 0;
                scalar avgWeight = 0;
                const point p = mesh_.C()[celli]; // cell center, celli here is also the global index

                forAll(stencil[celli], i)
                {
                    const label gblIdx = stencil[celli][i];                    
                    vector n = -distribute.getValue(normal, mapNormal, gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);
                        vector c = distribute.getValue(centre,mapCentres,gblIdx); // the problem is here, distance should be modified. 
                        distribute.centerTransformation(celli, gblIdx, c);
                        vector distanceToIntSeg = (c - p);
                        scalar distToSurf = distanceToIntSeg & (n);
                        scalar weight = 0;

                        if (mag(distanceToIntSeg) != 0)
                        {
                            distanceToIntSeg /= mag(distanceToIntSeg);
                            weight = sqr(mag(distanceToIntSeg & n));
                        }
                        else // exactly on the center
                        {
                            weight = 1;
                        }
                        averageDist += distToSurf * weight;
                        avgWeight += weight;
                    }
                }

                if (avgWeight != 0)
                {
                    reconDistFunc[celli] = averageDist / avgWeight;
                }
            }
        }
        else // nextToInterface[celli] == false: cells outer narrow band
        {
            reconDistFunc[celli] = 0;
        }
    }

    forAll(reconDistFunc.boundaryField(), patchI)
    {
        fvPatchScalarField& pRDF = reconDistFunc.boundaryFieldRef()[patchI];
        if (isA<calculatedFvPatchScalarField>(pRDF))
        {
            
            const polyPatch& pp = pRDF.patch().patch();
//            Info << "## show calculated patch name: " << pRDF.patch().name() << endl;
            forAll(pRDF, i)
            {                
                const label pCellI = pp.faceCells()[i];

                if (nextToInterface_[pCellI])
                {
                    scalar averageDist = 0;
                    scalar avgWeight = 0;
                    const point p = mesh_.C().boundaryField()[patchI][i];

                    forAll(stencil[pCellI], j)
                    {
                        const label gblIdx = stencil[pCellI][j];
                        vector n = -distribute.getValue(normal, mapNormal, gblIdx);
                        if (mag(n) != 0)
                        {
                            n /= mag(n);
                            vector c =
                                distribute.getValue(centre, mapCentres, gblIdx);
                            vector distanceToIntSeg = (c - p);
                            scalar distToSurf = distanceToIntSeg & (n);
                            scalar weight = 0;

                            if (mag(distanceToIntSeg) != 0)
                            {
                                distanceToIntSeg /= mag(distanceToIntSeg);
                                weight = sqr(mag(distanceToIntSeg & n));
                            }
                            else // exactly on the center
                            {
                                weight = 1;
                            }
                            averageDist += distToSurf * weight;
                            avgWeight += weight;
                        }
                    }

                    if (avgWeight != 0)
                    {
                        pRDF[i] = averageDist / avgWeight;
                    }
                    else
                    {
                        pRDF[i] = 0;
                    }
                }
                else
                {
                    pRDF[i] = 0;
                }
            }
        }
    }

    reconDistFunc.correctBoundaryConditions();

    return reconDistFunc;
}


void Foam::reconstructedDistanceFunction::updateContactAngle
(
    const volScalarField& alpha,
    const volVectorField& U,
    surfaceVectorField::Boundary& nHatb
)
{
    const fvMesh& mesh = alpha.mesh();
    const volScalarField::Boundary& abf = alpha.boundaryField();
    volScalarField::Boundary& RDFbf = this->boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad()*acap.theta(U.boundaryField()[patchi], nHatp)
            );

            RDFbf[patchi] =
                1/acap.patch().deltaCoeffs()*cos(theta)
              + RDFbf[patchi].patchInternalField();
        }
    }
}


// ************************************************************************* //
