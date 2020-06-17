"""Interactive flight track plot

Plot a variable of interest and the flight track then click on either figure
to mark the corresponding points on both figures.

> python -m twinotter.plots.interactive_flight_track /path/to/data
"""
import datetime
import tkinter
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import cartopy.crs as ccrs
import pandas as pd

from .. import load_flight
from . import plot_flight_path


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("flight_data_path")

    args = argparser.parse_args()

    start_gui(flight_data_path=args.flight_data_path)

    return


def start_gui(flight_data_path):
    ds = load_flight(flight_data_path)

    # Use pandas datetime functionality as xarray makes this difficult
    time = pd.to_datetime(ds.Time.data)
    # Flight leg times will be recorded as time since the start of the flight day
    flight_day_start = pd.to_datetime(ds.Time[0].dt.floor("D").data)

    root = tkinter.Tk()
    root.wm_title(
        "Interactive Flight Track: Flight {}".format(ds.attrs["flight_number"])
    )

    # Plot the main variable of interest
    # Change this to whatever variable you want or add additional figures here
    fig1, ax1a = plt.subplots()
    ax1a.plot(ds.Time, ds.ROLL_OXTS, linestyle="--", alpha=0.5)
    ax1a.set_label("Roll Angle")
    ax1b = ax1a.twinx()
    ax1b.plot(ds.Time, ds.ALT_OXTS / 1000)
    ax1b.set_ylabel("Altitude (km)")

    # Plot flight path with colours for altitude
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()),)
    ax2.gridlines(draw_labels=True)
    ax2.coastlines()
    plot_flight_path(ax=ax2, ds=ds)

    fig1.tight_layout()
    fig2.tight_layout()

    # Save flight leg start and end points
    leg_info = pd.DataFrame(columns=["Label", "Start", "End"])

    # Add the figures to as TK window
    figure_area = tkinter.Frame()
    figure_area.grid(row=0, column=0, columnspan=2)

    canvas = FigureCanvasTkAgg(fig1, master=figure_area)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0)

    canvas = FigureCanvasTkAgg(fig2, master=figure_area)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=1)

    # Add an area for buttons beneath the figures
    button_area = tkinter.Canvas(root)
    button_area.grid(row=1, column=1)

    def _save():
        filename = filedialog.asksaveasfilename(
            initialfile="flight{}-legs.csv".format(ds.attrs["flight_number"])
        )
        leg_info.to_csv(filename)

    save_button = tkinter.Button(master=button_area, text="Save", command=_save)
    save_button.grid(row=0, column=0)

    def _quit():
        root.quit()  # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent
        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

    quit_button = tkinter.Button(master=button_area, text="Quit", command=_quit)
    quit_button.grid(row=0, column=1)

    # Use an Entry textbox to label the legs
    textbox = tkinter.Entry(master=root)
    textbox.grid(row=1, column=0)

    # Add a span selector to the time-height plot to highlight legs
    # Drag mouse from the start to the end of a leg and save the corresponding
    # times
    def highlight_leg(start, end):
        nonlocal leg_info

        start = _convert_wacky_date_format(start)
        end = _convert_wacky_date_format(end)

        label = textbox.get()
        idx_start = find_nearest_point(start, time)
        idx_end = find_nearest_point(end, time)

        leg_info = leg_info.append(
            {
                "Label": label,
                "Start": str(time[idx_start] - flight_day_start),
                "End": str(time[idx_end] - flight_day_start),
            },
            ignore_index=True,
        )

        return

    selector = SpanSelector(ax1b, highlight_leg, direction="horizontal")

    tkinter.mainloop()


def find_nearest_point(value, points):
    return int(np.argmin(np.abs(value - points)))


t0 = datetime.datetime(1, 1, 1)


def _convert_wacky_date_format(wacky_time):
    # The twinotter MASIN data is loaded in with a datetime coordinate but when this is
    # used on the interactive plot the value returned from the click is in days from the
    # "zeroth" datetime. Use this zeroth datetime (t0) to get the date again.
    # The zeroth datetime is also for day=0 but this can't be handled by datetime so
    # we also have to subtract 1 day from the result
    return t0 + datetime.timedelta(days=wacky_time) - datetime.timedelta(days=1)


if __name__ == "__main__":
    main()
