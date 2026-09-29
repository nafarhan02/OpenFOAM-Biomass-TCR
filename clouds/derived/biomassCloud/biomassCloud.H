#ifndef biomassCloud_H
#define biomassCloud_H

#include "Cloud.H"
#include "KinematicCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "ReactingMultiphaseCloud.H"
#include "basicBiomassReactingMultiphaseParcel.H"

namespace Foam
{

typedef ReactingMultiphaseCloud
<
    ReactingCloud
    <
        ThermoCloud
        <
            KinematicCloud
            <
                Cloud
                <
                    basicBiomassReactingMultiphaseParcel
                >
            >
        >
    >
> biomassCloud;

}

#endif
