/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Class
    Foam::BiomassReactingMultiphaseParcel

Group
    grpLagrangianIntermediateParcels

Description
    Multiphase variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    BiomassReactingMultiphaseParcelI.H
    BiomassReactingMultiphaseParcel.C
    BiomassReactingMultiphaseParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BiomassReactingMultiphaseParcel_H
#define BiomassReactingMultiphaseParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType> class BiomassReactingMultiphaseParcel;

template<class ParcelType>
Ostream& operator<<(Ostream&, const BiomassReactingMultiphaseParcel<ParcelType>&);

/*---------------------------------------------------------------------------*\
                  Class BiomassReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class BiomassReactingMultiphaseParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    // IDs of phases in ReactingParcel phase list (Y)

        static const label GAS;
        static const label LIQ;
        static const label SLD;


    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Devolatilisation activation temperature [K]
            demandDrivenEntry<scalar> TDevol_;

            //- Latent heat of devolatilisation [J/kg]
            demandDrivenEntry<scalar> LDevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            demandDrivenEntry<scalar> hRetentionCoeff_;

            //- TCR first-reaction pre-exponential factor [1/s]
            demandDrivenEntry<scalar> TCRA1_;

            //- TCR first-reaction activation energy [J/kmol]
            demandDrivenEntry<scalar> TCRE1_;

            //- TCR second-reaction pre-exponential factor [1/s]
            demandDrivenEntry<scalar> TCRA2_;

            //- TCR second-reaction activation energy [J/kmol]
            demandDrivenEntry<scalar> TCRE2_;

            //- TCR volatile yield from first reaction [-]
            demandDrivenEntry<scalar> TCRV1_;

            //- TCR volatile yield from second reaction [-]
            demandDrivenEntry<scalar> TCRV2_;

            //- Remaining volatile fraction below which char may combust [-]
            demandDrivenEntry<scalar> TCRResidualCoeff_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the devolatilisation temperature
            inline scalar TDevol() const;

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;

            //- Return first TCR pre-exponential factor
            inline scalar TCRA1() const;

            //- Return first TCR activation energy
            inline scalar TCRE1() const;

            //- Return second TCR pre-exponential factor
            inline scalar TCRA2() const;

            //- Return second TCR activation energy
            inline scalar TCRE2() const;

            //- Return first TCR volatile yield
            inline scalar TCRV1() const;

            //- Return second TCR volatile yield
            inline scalar TCRV2() const;

            //- Return TCR volatile residual coefficient
            inline scalar TCRResidualCoeff() const;
    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


private:

    // Private Member Functions

        //- Return the mixture effective specific heat capacity
        template<class TrackCloudType>
        scalar CpEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackCloudType>
        scalar HsEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackCloudType>
        scalar LEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractions
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );


protected:

    // Protected data

        // Parcel properties

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to identify if the particle can devolatilise and combust
            //  Combustion possible only after volatile content falls below
            //  threshold value.  States include:
            //  0 = can devolatilise, cannot combust but can change
            //  1 = can devolatilise, can combust
            // -1 = cannot devolatilise or combust, and cannot change
            label canCombust_;

	    //- Remaining unreacted raw-fuel fraction for TCR devolatilisation
	    //  1 = fresh raw fuel, 0 = completely reacted raw fuel
	    scalar zeta_;


    // Protected Member Functions


        //- Return change of volume due to mass exchange
        template<class TrackCloudType>
        scalar updatedDeltaVolume
        (
            TrackCloudType& cloud,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid,
            const label idG,
            const label idL,
            const label idS,
            const scalar p,
            const scalar T
        );

        //- Calculate Devolatilisation
        template<class TrackCloudType>
        void calcDevolatilisation
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar age,          // age
            const scalar Ts,           // surface temperature
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const scalar mass0,        // mass (initial on injection)
            const scalarField& YGasEff,// gas component mass fractions
            const scalarField& YLiquidEff,// liquid component mass fractions
            const scalarField& YSolidEff,// solid component mass fractions
            label& canCombust,          // 'can combust' flag
            scalarField& dMassDV,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Calculate surface reactions
        template<class TrackCloudType>
        void calcSurfaceReactions
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar d,            // diameter
            const scalar Re,           // Re
            const scalar nu,           // nu
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassSRGas,   // gas-phase mass transfer - local
            scalarField& dMassSRLiquid,// liquid-phase mass transfer - local
            scalarField& dMassSRSolid, // solid-phase mass transfer - local
            scalarField& dMassSRCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans           // sensible enthalpy transfer to carrier
        ) const;


public:

    // Static data members

        //- Runtime type information
        TypeName("BiomassReactingMultiphaseParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " nGas(Y1..YN)"
          + " nLiquid(Y1..YN)"
          + " nSolid(Y1..YN)"
	  + " zeta"
        );


    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline BiomassReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline BiomassReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline BiomassReactingMultiphaseParcel
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
        );

        //- Construct from Istream
        BiomassReactingMultiphaseParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        BiomassReactingMultiphaseParcel(const BiomassReactingMultiphaseParcel& p);

        //- Construct as a copy
        BiomassReactingMultiphaseParcel
        (
            const BiomassReactingMultiphaseParcel& p,
            const polyMesh& mesh
        );

        //- Return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return particle::Clone(*this);
        }

        //- Return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return particle::Clone(*this, mesh);
        }

        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<BiomassReactingMultiphaseParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<BiomassReactingMultiphaseParcel<ParcelType>>
                (
                    new BiomassReactingMultiphaseParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to the canCombust flag
            inline label canCombust() const;

            //- Return remaining unreacted raw-fuel fraction
            inline scalar zeta() const;

        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to the canCombust flag
            inline label& canCombust();

	    //- Return access to remaining unreacted raw-fuel fraction
	    inline scalar& zeta();

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read - composition supplied
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write - composition supplied
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
            ) const;

            //- Read particle fields as objects from the obr registry
            //  - no composition
            template<class CloudType>
            static void readObjects
            (
                CloudType& c,
                const objectRegistry& obr
            );

            //- Read particle fields as objects from the obr registry
            template<class CloudType, class CompositionType>
            static void readObjects
            (
                CloudType& c,
                const CompositionType& compModel,
                const objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            //  - no composition
            template<class CloudType>
            static void writeObjects
            (
                const CloudType& c,
                objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            template<class CloudType, class CompositionType>
            static void writeObjects
            (
                const CloudType& c,
                const CompositionType& compModel,
                objectRegistry& obr
            );


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const BiomassReactingMultiphaseParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BiomassReactingMultiphaseParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BiomassReactingMultiphaseParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
