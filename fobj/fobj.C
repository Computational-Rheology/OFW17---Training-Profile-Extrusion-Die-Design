/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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

#include "fobj.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fobj, 0);
    addToRunTimeSelectionTable(functionObject, fobj, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fobj::fobj
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    //boolData_(dict.getOrDefault<bool>("boolData", true)),
    nES_(dict.get<label>("nES"))
    //,
    //wordData_(dict.getOrDefault<word>("wordData", "defaultWord")),
    //scalarData_(dict.getOrDefault<scalar>("scalarData", 1.0))
{
    read(dict);
        
    // Open the file
    outputFilePtr_.reset(new OFstream(time_.path()+"/Fobj"));

    // Write
    outputFilePtr_()<<  "nES (";
    
    for (label i = 1; i <= nES_ ; i++)
    {
        outputFilePtr_()<<  "Fobj_ES0" << i;
        if (i < nES_) outputFilePtr_()<<  ", ";    
    }

    outputFilePtr_()<<  ") Fobj_Tot , avVeL" << endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fobj::read(const dictionary& dict)
{
    dict.readEntry("nES", nES_);
    return true;
}


bool Foam::functionObjects::fobj::execute()
{
    return true;
}


bool Foam::functionObjects::fobj::end()
{
    return true;
}


bool Foam::functionObjects::fobj::write()
{
    const surfaceScalarField& phi_ = mesh_.lookupObject<surfaceScalarField>("phi");
 
    DynamicList<scalar> flowES;
    DynamicList<scalar> areaES;
    DynamicList<scalar> fobjES;
    
    scalar flowTot=0.0;
    scalar areaTot=0.0;

    for (label i = 1; i <= nES_ ; i++)
    {
        word lPatch = "es";
        if (i <= 9) lPatch += "0";
        lPatch += Foam::name(i);

        label patchID = mesh_.boundaryMesh().findPatchID(lPatch); // Look for faceZone instead....

        if ( patchID < 0 )
        {
            FatalErrorInFunction
                << "Unable to find patch " << lPatch 
                << abort(FatalError);
        }
         
        scalar flowPatch = mag(gSum(phi_.boundaryField()[patchID])) ;
        scalar areaPatch = mag(gSum(mesh_.Sf().boundaryField()[patchID]));

        flowES.append(flowPatch);
        flowTot += flowPatch;
        areaES.append(areaPatch);
        areaTot += areaPatch;
    }

    scalar aveVel = flowTot / areaTot;

    scalar fobjTot = 0.0;
    for (label i = 0; i < nES_ ; i++)
    {
        scalar fobjESi = ( flowES[i] / areaES[i] / aveVel - 1.0)  / max ( flowES[i] / areaES[i] , 1.0);
        fobjES.append( fobjESi ); 
        fobjTot += mag(fobjESi)*areaES[i];
    }
    
    fobjTot /= areaTot;
    
    outputFilePtr_() << nES_ << " (";

    for (label i = 0; i < nES_ ; i++)
    {
        outputFilePtr_() <<  fobjES[i] << i ;
        if (i < (nES_ - 1)) outputFilePtr_() <<  ", " ;    
    }

    outputFilePtr_() <<  " ) " << fobjTot << ", " << aveVel << endl;
       
    return true;
}


// ************************************************************************* //
