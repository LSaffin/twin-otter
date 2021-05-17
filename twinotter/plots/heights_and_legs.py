from pathlib import Path

from tqdm import tqdm
from adjustText import adjust_text
import xarray as xr
import matplotlib.pyplot as plt

from .. import load_flight, load_segments


colors = {
    "level": "cyan",
    "profile": "magenta",
}


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("flight_data_path")
    argparser.add_argument("flight_segments_file")
    argparser.add_argument("--show-gui", default=False, action="store_true")
    argparser.add_argument("--output_path", default=None)
    args = argparser.parse_args()

    generate(
        args.flight_data_path,
        args.flight_segments_file,
        show_gui=args.show_gui,
        output_path=args.output_path,
    )


def generate(flight_data_path, flight_segments_file, show_gui=False, output_path=None):
    ds = load_flight(flight_data_path)
    segments = load_segments(flight_segments_file)

    # Produce the basic time-height plot
    fig, ax1 = plt.subplots()
    ax1.plot(ds.Time, ds.ALT_OXTS / 1000, color="k")
    ax1.set_ylabel("Altitude (km)")

    # For each segment overlay a coloured line onto the time-height plot
    texts = []
    for segment in tqdm(segments["segments"]):
        ds_section = ds.sel(Time=slice(segment["start"], segment["end"]))
        color, linestyle = segment_linstyle(segment)

        ax1.plot(
            ds_section.Time,
            ds_section.ALT_OXTS / 1000,
            linestyle=linestyle,
            color=color,
            linewidth=2,
            alpha=0.75,
        )

        texts.append(
            ax1.text(
                ds_section.Time[0],
                ds_section.ALT_OXTS[0] / 1000,
                segment["segment_id"][7:],
                color=color,
            )
        )

    if hasattr(ax1, "secondary_yaxis"):
        # `ax.secondary_yaxis` was added in matplotlib v3.1
        ax1_fl = ax1.secondary_yaxis(
            location=-0.15, functions=(lambda y: (y * 1000 * 3.281) / 100, lambda x: x)
        )
        ax1_fl.set_ylabel(r"Flight level [100ft]")

    for label in ax1.get_xmajorticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment("right")

    adjust_text(texts)

    if show_gui:
        plt.show()
    else:
        if output_path is None:
            output_path = (
                Path(flight_data_path) / "figures" / "height-time-with-legs.png"
            )
        else:
            output_path = Path(output_path)

        output_path.parent.mkdir(exist_ok=True)
        plt.savefig(str(output_path), bbox_inches="tight")


def segment_linstyle(seg):
    if "calibration" in seg["kinds"]:
        return "red", "--"
    elif "bco_flyby" in seg["kinds"]:
        return "red", "-"

    elif "above_cloud" in seg["kinds"]:
        return "orange", "--"
    elif "cloud_base" in seg["kinds"]:
        return "grey", "--"
    elif "cloud" in seg["kinds"]:
        return "grey", "-"

    elif "cold_pool" in seg["kinds"]:
        return "blue", "-"
    elif "rain_shafts" in seg["kinds"]:
        return "blue", "--"
    elif "boundary_layer" in seg["kinds"]:
        return "green", "-"

    elif "sawtooth" in seg["kinds"]:
        return "cyan", "--"

    elif "transit" in seg["kinds"]:
        return "red", ":"
    elif "profile" in seg["kinds"]:
        return "magenta", "--"
    else:
        return "brown", "-"


if __name__ == "__main__":
    main()
