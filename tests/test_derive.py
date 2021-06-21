import pytest

import numpy as np
import xarray as xr

import twinotter
import twinotter.derive


def test_calculate(testdata):
    ds = twinotter.load_flight(testdata["flight_data_path"])

    theta = twinotter.derive.calculate("air_potential_temperature", ds)

    # This number is correct for the testdata using r004
    assert np.isclose(theta[0], 298.76865)


def test_calculate_nonsense(testdata):
    ds = twinotter.load_flight(testdata["flight_data_path"])

    with pytest.raises(ValueError):
        twinotter.derive.calculate("nonsense", ds)


@pytest.mark.parametrize(
    "variable,function,arguments",
    [
        (
            "air_temperature",
            twinotter.derive.combine_temperatures,
            ["TAT_ND_R", "TAT_DI_R"],
        )
    ],
)
def test_calculate_equivalent(testdata, variable, function, arguments):
    # Check that "calculate" returns an equivalent variable to calling the specific
    # function Currently only setup for "first-order" functions. i.e. functions that
    # calculate a new variable using only variables that are already in the basic
    # MASIN dataset
    ds = twinotter.load_flight(testdata["flight_data_path"])

    result1 = twinotter.derive.calculate(variable, ds)
    result2 = function(*[ds[arg] for arg in arguments])

    result1 = _filter_nans(result1)
    result2 = _filter_nans(result2)

    assert isinstance(result1, xr.core.dataarray.DataArray)
    assert isinstance(result2, xr.core.dataarray.DataArray)

    assert (result1 == result2).all()


def _filter_nans(array):
    # NaNs mess up array-wise equality checks
    return array.where(~np.isnan(array), drop=True)


def test_liquid_water_content(testdata):
    ds = xr.open_dataset(testdata["flight_microphysics_file"])
    lwc = twinotter.derive.liquid_water_content(
        ds.ambient_particle_number_per_channel, ds.ambient_particle_diameter / 2
    )

    assert isinstance(lwc, xr.core.dataarray.DataArray)
