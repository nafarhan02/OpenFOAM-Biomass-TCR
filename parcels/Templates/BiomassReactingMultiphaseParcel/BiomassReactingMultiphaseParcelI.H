/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    TDevol_(this->dict_, 0.0),
    LDevol_(this->dict_, 0.0),
    hRetentionCoeff_(this->dict_, 0.0),
    TCRA1_(this->dict_, 0.0),
    TCRE1_(this->dict_, 0.0),
    TCRA2_(this->dict_, 0.0),
    TCRE2_(this->dict_, 0.0),
    TCRV1_(this->dict_, 0.0),
    TCRV2_(this->dict_, 0.0),
    TCRResidualCoeff_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    TDevol_(cp.TDevol_),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_),
    TCRA1_(cp.TCRA1_),
    TCRE1_(cp.TCRE1_),
    TCRA2_(cp.TCRA2_),
    TCRE2_(cp.TCRE2_),
    TCRV1_(cp.TCRV1_),
    TCRV2_(cp.TCRV2_),
    TCRResidualCoeff_(cp.TCRResidualCoeff_)
{}


template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    TDevol_(this->dict_, "TDevol"),
    LDevol_(this->dict_, "LDevol"),
    hRetentionCoeff_(this->dict_, "hRetentionCoeff"),
    TCRA1_(this->dict_, "TCRA1"),
    TCRE1_(this->dict_, "TCRE1"),
    TCRA2_(this->dict_, "TCRA2"),
    TCRE2_(this->dict_, "TCRE2"),
    TCRV1_(this->dict_, "TCRV1"),
    TCRV2_(this->dict_, "TCRV2"),
    TCRResidualCoeff_(this->dict_, "TCRResidualCoeff")
{}


template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::BiomassReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0),
    zeta_(1.0)
{}


template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::BiomassReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0),
    zeta_(1.0)
{}


template<class ParcelType>
inline Foam::BiomassReactingMultiphaseParcel<ParcelType>::BiomassReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    ),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    canCombust_(0),
    zeta_(1.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TDevol() const
{
    return TDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::LDevol() const
{
    return LDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::
hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRA1() const
{
    return TCRA1_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRE1() const
{
    return TCRE1_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRA2() const
{
    return TCRA2_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRE2() const
{
    return TCRE2_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRV1() const
{
    return TCRV1_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::TCRV2() const
{
    return TCRV2_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::constantProperties::
TCRResidualCoeff() const
{
    const scalar value = TCRResidualCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "TCRResidualCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label
Foam::BiomassReactingMultiphaseParcel<ParcelType>::canCombust() const
{
    return canCombust_;
}

template<class ParcelType>
inline Foam::scalar
Foam::BiomassReactingMultiphaseParcel<ParcelType>::zeta() const
{
    return zeta_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::BiomassReactingMultiphaseParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label& Foam::BiomassReactingMultiphaseParcel<ParcelType>::canCombust()
{
    return canCombust_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::BiomassReactingMultiphaseParcel<ParcelType>::zeta()
{
    return zeta_;
}

// ************************************************************************* //
