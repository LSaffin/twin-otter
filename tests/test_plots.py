from unittest.mock import patch
import pytest
from pathlib import Path

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import twinotter.plots.basic_flight_track
import twinotter.plots.flight_track_frames
import twinotter.plots.vertical_profile
import twinotter.plots.heights_and_legs
import twinotter.quicklook
import twinotter.external.goes


@patch("matplotlib.pyplot.savefig")
def test_basic_flight_path(mock_savefig, testdata):
    twinotter.plots.basic_flight_track.generate(
        flight_data_path=testdata["flight_data_path"]
    )
    mock_savefig.assert_called_once()


def test_flight_track_frame(testdata):
    ds = twinotter.external.goes.load_nc(
        path=testdata["goes_path"],
        time=testdata["goes_time"],
    )
    fig, ax = twinotter.plots.flight_track_frames.make_frame(ds)

    plt.close(fig)


def test_overlay_flight_path_segment(testdata):
    twinotter.plots.flight_track_frames.overlay_flight_path_segment(
        ax=plt.axes(projection=ccrs.PlateCarree()),
        flight_data=twinotter.load_flight(testdata["flight_data_path"]),
        time=testdata["goes_time"],
    )

    plt.close()


@patch("matplotlib.pyplot.show")
def test_vertical_profile_plot(mock_showfig, testdata):
    twinotter.plots.vertical_profile.main(flight_data_path=testdata["flight_data_path"])
    mock_showfig.assert_called_once()


@patch("matplotlib.pyplot.savefig")
def test_heights_and_legs_plot(mock_savefig, testdata):
    twinotter.plots.heights_and_legs.generate(
        flight_data_path=testdata["flight_data_path"],
        flight_segments_file=testdata["flight_segments_file"],
    )
    mock_savefig.assert_called_once()


@patch("matplotlib.figure.Figure.savefig")
def test_quicklook_plot(mock_savefig, testdata):
    twinotter.quicklook.generate(
        flight_data_path=testdata["flight_data_path"],
        flight_segments_file=testdata["flight_segments_file"],
    )

    flight_segments = twinotter.load_segments(testdata["flight_segments_file"])

    for seg in flight_segments["segments"]:
        seg_id = seg["segment_id"]
        if "level" in seg["kinds"]:
            for figname in ["quicklook", "paluch"]:
                target_filename = "{}_{}.png".format(seg_id, figname)
                mock_savefig.assert_any_call(target_filename)
        if "profile" in seg["kinds"]:
            for figname in ["skewt", "rh_0-4km", "theta_0-4km"]:
                target_filename = "{}_{}.png".format(seg_id, figname)
                mock_savefig.assert_any_call(target_filename)

    for figname in ["quicklook", "paluch"]:
        target_filename = "{}_{}.png".format("TO-330_overview", figname)
        mock_savefig.assert_any_call(target_filename)

    for figname in ["skewt", "rh_0-4km", "theta_0-4km"]:
        target_filename = "{}_{}.png".format("TO-330_profiles", figname)
        mock_savefig.assert_any_call(target_filename)


# this doesn't run properly because of the gui needing to be opened
# @patch('matplotlib.pyplot.show')
# def test_interactive_flight_path(mock_show):
# twinotter.plots.interactive_flight_track.start_gui(flight_data_path="obs/flight330")
# mock_show.assert_called_once()
