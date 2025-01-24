import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from derotation.analysis.full_derotation_pipeline import FullPipeline
from derotation.analysis.metrics import stability_of_most_detected_blob
from matplotlib import pyplot as plt
from snakemake.script import snakemake

datasets = snakemake.params.datasets
processed_data_base = snakemake.params.base_path

csv_path = Path(snakemake.output[0]).with_suffix(".csv")
img_path = Path(snakemake.output[0])

if not img_path.exists():
    datasets_paths = []
    for idx, dataset in enumerate(datasets):
        datasets_paths.append(
            Path(
                f"{snakemake.params.base_path}/sub-{idx}_{dataset}/ses-0/funcimg"
            )
        )

    movie_bin_paths = []
    for dataset in datasets_paths:
        movie_bin_paths.extend(list(Path(dataset).rglob("*.bin")))

    is_cell_paths = []
    for dataset in datasets_paths:
        is_cell_paths.extend(list(Path(dataset).rglob("iscell.npy")))

    metric_paths = []
    for path_to_bin_file in movie_bin_paths:
        metric_path = (
            path_to_bin_file.parent.parent.parent / "derotation_metrics.csv"
        )
        metric_paths.append(metric_path)

    derotated_full_csv_paths = []
    for path_to_bin_file in movie_bin_paths:
        derotated_full_csv_path = (
            path_to_bin_file.parent.parent.parent
            / "derotation/derotated_full.csv"
        )
        derotated_full_csv_paths.append(derotated_full_csv_path)

    all_metrics_df = pd.DataFrame(
        columns=["dataset", "analysis_type", "metric", "value"]
    )
    analysis_types = [
        "no_adj",
        "adj_track",
        "adj_largest",
        "adj_track_shifted",
    ]

    for path_to_bin_file, metric_path, derotated_full_csv_path in zip(
        movie_bin_paths, metric_paths, derotated_full_csv_paths
    ):
        print(
            f"Processing dataset: {path_to_bin_file.parent.parent.parent.parent.parent.name}..."
        )
        try:
            metric = pd.read_csv(metric_path)

            path_to_bin_file = Path(path_to_bin_file)

            rotation_df = pd.read_csv(derotated_full_csv_path)
            num_frames = len(rotation_df)

            shape_image = (num_frames, 256, 256)
            registered = np.memmap(
                path_to_bin_file, shape=shape_image, dtype="int16"
            )

            #  plot first frame as an image of registered as a way to test if the loading was correct
            plt.imshow(registered[0])
            plt.savefig(path_to_bin_file.parent / "first_frame_registered.png")
            plt.close()

            derotator = FullPipeline.__new__(FullPipeline)

            angles = rotation_df["rotation_angle"].values
            if len(angles) > len(registered):
                angles = angles[: len(registered)]
            elif len(angles) < len(registered):
                angles = np.pad(angles, (0, len(registered) - len(angles)))

            derotator.rot_deg_frame = angles
            mean_images = derotator.calculate_mean_images(
                registered, round_decimals=0
            )

            #  show first mean image
            plt.imshow(mean_images[0])
            plt.savefig(path_to_bin_file.parent / "first_mean_image.png")
            plt.close()

            path_plots = path_to_bin_file.parent
            try:
                ptd, std = stability_of_most_detected_blob(
                    (mean_images, path_plots),
                    # blob_log_kwargs={"min_sigma": 0, "max_sigma": 20, "threshold": 0.5, "overlap": 0},
                    # clip=False
                )
                print(f"ptd: {ptd}, std: {std}")
            except Exception as e:
                print(e)
                print(traceback.format_exc())
                ptd = np.nan
                std = np.nan

            for i, analysis_type in enumerate(analysis_types):
                row_ptd = {
                    "dataset": path_to_bin_file.parent.parent.parent.parent.parent.name,
                    "analysis_type": analysis_type,
                    "metric": "ptd",
                    "value": metric["ptd"][i],
                }
                row_std = {
                    "dataset": path_to_bin_file.parent.parent.parent.parent.parent.name,
                    "analysis_type": analysis_type,
                    "metric": "std",
                    "value": metric["std"][i],
                }
                all_metrics_df = pd.concat(
                    [all_metrics_df, pd.DataFrame([row_ptd, row_std])],
                    ignore_index=True,
                )
            #  add post_suite2p metrics
            row_ptd = {
                "dataset": path_to_bin_file.parent.parent.parent.parent.parent.name,
                "analysis_type": "post_suite2p",
                "metric": "ptd",
                "value": ptd,
            }
            row_std = {
                "dataset": path_to_bin_file.parent.parent.parent.parent.parent.name,
                "analysis_type": "post_suite2p",
                "metric": "std",
                "value": std,
            }
            all_metrics_df = pd.concat(
                [all_metrics_df, pd.DataFrame([row_ptd, row_std])],
                ignore_index=True,
            )
        except Exception as e:
            print(e)
            print("Error in dataset")
            continue

    #  save the dataframe to a csv file (change the file extension from png to csv)
    all_metrics_df.to_csv(csv_path, index=False)

else:
    all_metrics_df = pd.read_csv(csv_path)

sns.set_theme(style="whitegrid")
sns.set_context("paper")
sns.set_palette("pastel")

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

sns.pointplot(
    x="analysis_type",
    y="value",
    hue="dataset",
    data=all_metrics_df[all_metrics_df["metric"] == "ptd"],
    ax=axs[0],
)

sns.pointplot(
    x="analysis_type",
    y="value",
    hue="dataset",
    data=all_metrics_df[all_metrics_df["metric"] == "std"],
    ax=axs[1],
)

axs[0].set_title("PTD")
axs[1].set_title("STD")

plt.tight_layout()
plt.savefig(img_path)
plt.close()

# make another similar plot with these analysis types:
# 1. no_adj
# 2. the min between "no_adj", "adj_track", "adj_largest", "adj_track_shifted" (to be calculated)
# 3. post_suite2p

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

data = pd.DataFrame(columns=["dataset", "analysis_type", "metric", "value"])
for dataset in all_metrics_df["dataset"].unique():
    dataset_df = all_metrics_df[all_metrics_df["dataset"] == dataset]
    #  no_adj
    no_adj_value = dataset_df[
        (dataset_df["analysis_type"] == "no_adj")
        & (dataset_df["metric"] == "ptd")
    ]["value"].values[0]
    row = {
        "dataset": dataset,
        "analysis_type": "no_adj",
        "metric": "ptd",
        "value": no_adj_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)
    #  min but not for post_suite2p
    min_value = dataset_df[
        (dataset_df["analysis_type"] != "post_suite2p")
        & (dataset_df["metric"] == "ptd")
    ]["value"].min()
    row = {
        "dataset": dataset,
        "analysis_type": "min",
        "metric": "ptd",
        "value": min_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)
    #  post_suite2p
    post_suite2p_value = dataset_df[
        (dataset_df["analysis_type"] == "post_suite2p")
        & (dataset_df["metric"] == "ptd")
    ]["value"].values[0]
    row = {
        "dataset": dataset,
        "analysis_type": "post_suite2p",
        "metric": "ptd",
        "value": post_suite2p_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)

#  save dataset
data.to_csv(csv_path.with_name("min_analysis_types_min_ptd.csv"), index=False)

sns.pointplot(
    x="analysis_type",
    y="value",
    hue="dataset",
    data=data[data["metric"] == "ptd"],
    ax=axs[0],
)


axs[0].set_ylabel("Point to point distance (r)")
axs[0].set_xlabel("Derotation adjustment")
axs[0].set_xticklabels(["No", "Yes", "Post Suite2p"])

data = pd.DataFrame(columns=["dataset", "analysis_type", "metric", "value"])

for dataset in all_metrics_df["dataset"].unique():
    dataset_df = all_metrics_df[all_metrics_df["dataset"] == dataset]
    #  no_adj
    no_adj_value = dataset_df[
        (dataset_df["analysis_type"] == "no_adj")
        & (dataset_df["metric"] == "std")
    ]["value"].values[0]
    row = {
        "dataset": dataset,
        "analysis_type": "no_adj",
        "metric": "std",
        "value": no_adj_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)
    #  min but not for post_suite2p
    min_value = dataset_df[
        (dataset_df["analysis_type"] != "post_suite2p")
        & (dataset_df["metric"] == "std")
    ]["value"].min()
    row = {
        "dataset": dataset,
        "analysis_type": "min",
        "metric": "std",
        "value": min_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)
    #  post_suite2p
    post_suite2p_value = dataset_df[
        (dataset_df["analysis_type"] == "post_suite2p")
        & (dataset_df["metric"] == "std")
    ]["value"].values[0]
    row = {
        "dataset": dataset,
        "analysis_type": "post_suite2p",
        "metric": "std",
        "value": post_suite2p_value,
    }
    data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)

#  save dataset
data.to_csv(csv_path.with_name("min_analysis_types_min_std.csv"), index=False)

sns.pointplot(
    x="analysis_type",
    y="value",
    hue="dataset",
    data=data[data["metric"] == "std"],
    ax=axs[1],
)

axs[1].set_ylabel("XY standard deviation (s)")
axs[1].set_xlabel("Derotation adjustment")
axs[1].set_xticklabels(["No", "Yes", "Post Suite2p"])

#  remove legend
axs[0].get_legend().remove()
axs[1].get_legend().remove()

axs[0].set_title("PTD")
axs[1].set_title("STD")

#  despine
sns.despine()

plt.tight_layout()
plt.savefig(img_path.with_name("min_analysis_types.png"))
plt.savefig(img_path.with_name("min_analysis_types.pdf"))
