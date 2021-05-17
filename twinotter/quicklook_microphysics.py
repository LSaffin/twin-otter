"""Quicklook plots for each segment over a single flight.

Use the flight-segments .yaml produced from
:mod:`twinotter.plots.interactive_flight_track`

Usage::

    $ python -m twinotter.quicklook_microphysics
        <flight_data_path> <flight_segments_file>

"""

from pathlib import Path

import parse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from . import MICROPHYSICS_FORMAT, load_segments, matching_segments
from .quicklook import savefigs


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("flight_data_path")
    argparser.add_argument("flight_segments_file")

    args = argparser.parse_args()

    generate(
        flight_data_path=args.flight_data_path,
        flight_segments_file=args.flight_segments_file,
    )

    return


def generate(flight_data_path, flight_segments_file):
    filenames = list(
        Path(flight_data_path).glob(
            MICROPHYSICS_FORMAT.format(flight_num="*", instrument="*", revision=1)
        )
    )

    datasets = dict()
    for filename in filenames:
        meta = parse.parse(MICROPHYSICS_FORMAT, str(filename.name)).named
        datasets[meta["instrument"]] = xr.open_dataset(filename)

    flight_segments = load_segments(flight_segments_file)

    # Quicklook plots for the full flight
    figures = plot_microphysics(datasets)
    for fig, figname in figures:
        fig.suptitle("{}: {}".format(flight_segments["flight_id"], figname))
    savefigs(figures, "{}_overview".format(flight_segments["flight_id"]))

    plot_individual_phases(datasets, flight_segments, "level", plot_microphysics)
    plot_individual_phases(datasets, flight_segments, "profile", plot_microphysics)


def plot_individual_phases(datasets, flight_segments, segment_type, plot_func):
    segments = matching_segments(flight_segments, segment_type)
    for seg in segments:
        ds_sections = {}
        for instrument in datasets:
            ds = datasets[instrument]
            ds_sub = ds.sel(time=slice(seg["start"], seg["end"]))

            if len(ds_sub.time) > 0:
                ds_sections[instrument] = ds_sub

        figures = plot_func(ds_sections)
        for fig, figname in figures:
            fig.suptitle(
                "{}: {} {}".format(flight_segments["flight_id"], seg["name"], figname)
            )
        savefigs(figures, seg["segment_id"])


def plot_microphysics(datasets):
    figures = []

    n_micro = len(datasets)

    fig, axes = plt.subplots(
        nrows=n_micro + 1, ncols=1, sharex="all", figsize=[8, 2.5 * (n_micro + 1)]
    )

    ds = datasets[list(datasets.keys())[0]]
    axes[0].plot(ds.time, ds.altitude)
    axes[0].set_ylabel("Altitude (m)")

    for n, instrument in enumerate(datasets):
        ds = datasets[instrument]
        plt.axes(axes[n + 1])
        plt.pcolormesh(
            ds.time,
            ds.ambient_particle_diameter,
            np.log(ds.ambient_particle_number_per_channel.transpose()),
            cmap="cubehelix_r",
        )

        axes[n + 1].set_ylabel(r"Ambient particle diameter ($\mu$m)")
        axes[n + 1].set_title(instrument + " ambient particle number")

        plt.colorbar(orientation="horizontal")

    axes[-1].set_xlabel("Time")

    figures.append([fig, "microphysics"])

    return figures


if __name__ == "__main__":
    main()
