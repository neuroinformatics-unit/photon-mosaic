from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from snakemake.script import snakemake

print("Plotting data...")


# datasets = snakemake.params.datasets
f_path = Path(snakemake.input[0])

print(f"Processing file: {f_path}")

dataset_path = f_path.parent.parent.parent.parent
dataset = dataset_path.name
f_neu_path = f_path.parent / "Fneu.npy"
derotated_full_csv_path = (
    dataset_path / "ses-0" / "funcimg" / "derotation" / "derotated_full.csv"
)
saving_path = Path(dataset_path) / "ses-0" / "traces"
saving_path.mkdir(exist_ok=True)

print(f"Dataset path: {dataset_path}")
print(f"Dataset: {dataset}")
print(f"Derotated full csv path: {derotated_full_csv_path}")
print(f"Fneu path: {f_neu_path}")
print(f"Saving path: {saving_path}")


f = np.load(f_path)
fneu = np.load(f_neu_path)
rotated_frames = pd.read_csv(derotated_full_csv_path)
f_corrected = f - 0.7 * fneu


F_df = pd.DataFrame(f_corrected).T

print(f"Shape of F_df: {F_df.shape}")
print(F_df.head())

full_dataframe = pd.concat([F_df, rotated_frames], axis=1)

# --------------------------------------------------------
# Prepare the dataset

# find where do rotations start
rotation_on = np.diff(full_dataframe["rotation_count"])


def find_zero_chunks(arr):
    zero_chunks = []
    start = None

    for i in range(len(arr)):
        if arr[i] == 0 and start is None:
            start = i
        elif arr[i] != 0 and start is not None:
            zero_chunks.append((start, i - 1))
            start = None

    # Check if the array ends with a chunk of zeros
    if start is not None:
        zero_chunks.append((start, len(arr) - 1))

    return zero_chunks


starts_ends = find_zero_chunks(rotation_on)

frames_before_rotation = 15
# frames_after_rotation = 10

total_len = 100

full_dataframe["rotation_frames"] = np.zeros(len(full_dataframe))
for i, (start, end) in enumerate(starts_ends):
    frame_array = np.arange(total_len)
    column_index_of_rotation_frames = full_dataframe.columns.get_loc(
        "rotation_frames"
    )
    full_dataframe.iloc[
        start - frames_before_rotation : total_len
        + start
        - frames_before_rotation,
        column_index_of_rotation_frames,
    ] = frame_array

    #  extend this value of speed and direction to all this range
    this_speed = full_dataframe.loc[start, "speed"]
    this_direction = full_dataframe.loc[start, "direction"]

    full_dataframe.iloc[
        start - frames_before_rotation : total_len
        + start
        - frames_before_rotation,
        full_dataframe.columns.get_loc("speed"),
    ] = this_speed
    full_dataframe.iloc[
        start - frames_before_rotation : total_len
        + start
        - frames_before_rotation,
        full_dataframe.columns.get_loc("direction"),
    ] = this_direction


#  directtion, change -1 to CCW and 1 to CW
full_dataframe["direction"] = np.where(
    full_dataframe["direction"] == -1, "CCW", "CW"
)

# print(f"Full dataframe shape: {full_dataframe.shape}")
# print(full_dataframe.head())

# #  angle based calculation of ΔF/F
# #  first calculate F0, as the 20th quantile for each angle.
# #  consider angles every 5 degrees, from 0 to 355
# full_dataframe["aproximated_rotation_angle"] = (
#     full_dataframe["rotation_angle"] // 5 * 5
# )

# print("Unique angles:")
# print(full_dataframe["aproximated_rotation_angle"].unique())

# f0_as_20th_quantile_per_angle = np.zeros((360, f_corrected.shape[0]))
# for angle in range(360):
#     for roi in range(f_corrected.shape[0]):
#         angle_indices = full_dataframe["aproximated_rotation_angle"] == angle
#         print(f"Angle: {angle}, ROI: {roi}")
#         print(f"Angle indices: {angle_indices}")
#         #  check for nans / missing values in angle_indices
#         if angle_indices.isnull().values.any():
#             f0_as_20th_quantile_per_angle[angle, roi] = np.nan
#         else:
#             f0_as_20th_quantile_per_angle[angle, roi] = np.quantile(
#                 f_corrected[roi][angle_indices], 0.2
#             )
# print("Shape of f0_as_20th_quantile_per_angle:")
# print(f0_as_20th_quantile_per_angle.shape)
# print(f0_as_20th_quantile_per_angle)

# #  calculate ΔF/F
# for roi in range(f_corrected.T.shape[0]):
#     full_dataframe[roi] = (
#         f_corrected.T[roi] - f0_as_20th_quantile_per_angle[
#             full_dataframe["rotation_angle"], roi
#         ]
#     ) / f0_as_20th_quantile_per_angle[
#         full_dataframe["rotation_angle"], roi
#     ]

# print("Full dataframe with ΔF/F:")
# print(full_dataframe.head())

rois_selection = range(F_df.shape[1])

# --------------------------------------------------------
# Plot single traces

# %%
selected_range = (400, 2000)

for roi in rois_selection:
    roi_selected = full_dataframe.loc[
        :, [roi, "rotation_count", "speed", "direction"]
    ]

    fig, ax = plt.subplots(figsize=(27, 5))
    ax.plot(roi_selected.loc[selected_range[0] : selected_range[1], roi])
    ax.set(xlabel="Frames", ylabel="Neuropil corrected (a.u.)")  # "ΔF/F")

    rotation_on = (
        np.diff(
            roi_selected.loc[
                selected_range[0] : selected_range[1], "rotation_count"
            ]
        )
        == 0
    )

    # add label at the beginning of every block of rotations
    #  if the previous was true, do not write the label
    for i, rotation in enumerate(rotation_on):
        if rotation and not rotation_on[i - 1]:
            ax.text(
                i + selected_range[0] + 3,
                -1100,
                f"{int(roi_selected.loc[i + 5 + selected_range[0], 'speed'])}º/s\n{roi_selected.loc[i + 5 + selected_range[0], 'direction']}",
                fontsize=10,
            )

    #  add gray squares when the rotation is happening using the starst_ends
    for start, end in starts_ends:
        if start > selected_range[0] and end < selected_range[1]:
            ax.axvspan(start, end, color="gray", alpha=0.2)

    fps = 6.74
    # change xticks to seconds
    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    ax.set_xticklabels((xticks / fps).astype(int))
    #  change x label
    ax.set(xlabel="Seconds", ylabel="Neuropil corrected (a.u.)")  # "ΔF/F")

    ax.set_xlim(selected_range)
    # ax.set_ylim(-10, 10)

    # leave some gap between the axis and the plot
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    # remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.savefig(saving_path / f"dff_example_{roi}.pdf")
    plt.savefig(saving_path / f"dff_example_{roi}.png")
    plt.close()


# --------------------------------------------------------
# Plot averages

custom_palette = sns.color_palette("dark:#5A9_r", 4)

for roi in rois_selection:
    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    for i, direction in enumerate(["CW", "CCW"]):
        sns.lineplot(
            x="rotation_frames",
            y=roi,
            data=full_dataframe[(full_dataframe["direction"] == direction)],
            hue="speed",
            palette=custom_palette,
            ax=ax[i],
        )
        ax[i].set_title(f"Direction: {direction}")
        ax[i].legend(title="Speed")

        #  remove top and right spines
        ax[i].spines["top"].set_visible(False)
        ax[i].spines["right"].set_visible(False)

        # add vertical lines to show the start of the rotation
        #  start is always at 11, end at total len - 10
        ax[i].axvline(x=frames_before_rotation, color="gray", linestyle="--")

        #  change x axis to seconds
        fps = 6.74
        xticks = ax[i].get_xticks()
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(np.round(xticks / fps, 1))
        #  change x label
        ax[i].set(
            xlabel="Seconds", ylabel="Neuropil corrected (a.u.)"
        )  # "ΔF/F")

    plt.savefig(saving_path / f"roi_{roi}_direction_speed.pdf")
    plt.savefig(saving_path / f"roi_{roi}_direction_speed.png")
    plt.close()

    #  make also another plot showing all traces (not averaged - no std)

    fig, ax = plt.subplots(figsize=(20, 10))
    for i, direction in enumerate(["CW", "CCW"]):
        # sns.relplot(
        #     x="rotation_frames",
        #     y=roi,
        #     data=full_dataframe[(full_dataframe["direction"] == direction)],
        #     hue="speed",
        #     palette=custom_palette,
        #     kind="line",
        #     estimator=None,
        #     style="direction",
        #     ax=ax,
        # )
        #  plot single traces using matplotlib
        for speed in full_dataframe["speed"].unique():
            ax.plot(
                full_dataframe[
                    (full_dataframe["direction"] == direction)
                    & (full_dataframe["speed"] == speed)
                ]["rotation_frames"],
                full_dataframe[
                    (full_dataframe["direction"] == direction)
                    & (full_dataframe["speed"] == speed)
                ][roi],
                label=f"{speed}º/s",
                # color=custom_palette[speed],
            )

        ax.set_title(f"Direction: {direction}")
        ax.legend(title="Speed")

        #  remove top and right spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # add vertical lines to show the start of the rotation
        #  start is always at 11, end at total len - 10
        ax.axvline(x=frames_before_rotation, color="gray", linestyle="--")

        #  change x axis to seconds
        fps = 6.74
        xticks = ax.get_xticks()
        ax.set_xticks(xticks)
        ax.set_xticklabels(np.round(xticks / fps, 1))
        #  change x label
        ax.set(xlabel="Seconds", ylabel="Neuropil corrected (a.u.)")  # "ΔF/F")

    plt.savefig(saving_path / f"roi_{roi}_direction_speed_all.pdf")
    plt.savefig(saving_path / f"roi_{roi}_direction_speed_all.png")

    plt.close()
