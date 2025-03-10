import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from derotation.analysis.full_derotation_pipeline import FullPipeline
from derotation.analysis.metrics import ptd_of_most_detected_blob
from matplotlib import pyplot as plt
from snakemake.script import snakemake

from derotation.analysis.mean_images import calculate_mean_images


datasets = snakemake.params.datasets
processed_data_base = snakemake.params.base_path

csv_path = Path(snakemake.output[0]).with_suffix(".csv")
img_path = Path(snakemake.output[0])

if not csv_path.exists():
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
        columns=["dataset", "ptd_pre_suite2p", "ptd_post_suite2p", "#_is_cell"]
    )

    for path_to_bin_file, metric_path, derotated_full_csv_path, is_cell_path in zip(
        movie_bin_paths, metric_paths, derotated_full_csv_paths, is_cell_paths
    ):
        print(
            f"Processing dataset: {path_to_bin_file.parent.parent.parent.parent.parent.name}..."
        )
        try:
            with open(metric_path, "r") as f:
                line = f.readlines()
                ptd_pre_suite2p = float(line[0].split(":")[1].strip())

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

            angles = rotation_df["rotation_angle"].values
            if len(angles) > len(registered):
                angles = angles[: len(registered)]
            elif len(angles) < len(registered):
                angles = np.pad(angles, (0, len(registered) - len(angles)))

            mean_images = calculate_mean_images(registered, angles, round_decimals=0)

            #  show first mean image
            plt.imshow(mean_images[0])
            plt.savefig(path_to_bin_file.parent / "first_mean_image.png")
            plt.close()

            path_plots = path_to_bin_file.parent
            try:
                ptd = ptd_of_most_detected_blob(mean_images, plot=True, debug_plots_folder=path_plots) 
                print(f"ptd: {ptd}")
            except Exception as e:
                print(e)
                print(traceback.format_exc())
                ptd = np.nan

            is_cell = np.load(is_cell_path)[:, 0]
            n_is_cell = int(sum(is_cell))

            row = {
                "dataset": path_to_bin_file.parent.parent.parent.parent.parent.name,
                "ptd_pre_suite2p": ptd_pre_suite2p,
                "ptd_post_suite2p": ptd,
                "#_is_cell": n_is_cell,
            }

            print(f"This row: {row}")
            
            all_metrics_df = pd.concat(
                [all_metrics_df, pd.DataFrame([row])], ignore_index=True
            )


        except Exception as e:
            print(e)
            print("Error in dataset")
            continue

    #  save the dataframe to a csv file (change the file extension from png to csv)
    all_metrics_df.to_csv(csv_path, index=False)

else:
    all_metrics_df = pd.read_csv(csv_path)

# Melt the dataframe to long format
df_long = all_metrics_df.melt(id_vars=["dataset"], value_vars=["ptd_pre_suite2p", "ptd_post_suite2p"],
                  var_name="Condition", value_name="PTD Value")

# Rename conditions for better readability
df_long["Condition"] = df_long["Condition"].replace({"ptd_pre_suite2p": "Pre", "ptd_post_suite2p": "Post"})

# Plot the data
plt.figure(figsize=(12, 6))
for dataset in all_metrics_df["dataset"]:
    subset = df_long[df_long["dataset"] == dataset]
    plt.plot(subset["Condition"], subset["PTD Value"], marker="o", linestyle="-", label=dataset)

plt.xlabel("Condition")
plt.ylabel("PTD Value")
plt.title("Pre vs Post PTD Values per Dataset")
plt.xticks(["Pre", "Post"])
plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Dataset", fontsize=8)

plt.savefig(img_path)
print(f'Image saved at {img_path}')


cleaned_data = [
    {'Subject': '230801CAA1120181', 'Category': 1},
    {'Subject': '230802CAA1120182', 'Category': 1},
    {'Subject': '230803CAA1119915', 'Category': 2},
    {'Subject': '230803CAA1120181', 'Category': 1},
    {'Subject': '230804CAA1119917', 'Category': 2},
    {'Subject': '230818CAA1120210', 'Category': 2},
    {'Subject': '230822CAA1120509', 'Category': 1},
    {'Subject': '230823CAA1120181', 'Category': 1},
    {'Subject': '230824CAA1119915', 'Category': 2},
    {'Subject': '230825CAA1120182', 'Category': 1},
    {'Subject': '230907CAA1120210', 'Category': 2},
    {'Subject': '230907CAA1120509', 'Category': 1},
    # {'Subject': '230907CAA1120509', 'Category': 1},
    {'Subject': '230912CAA1119915', 'Category': 2},
    {'Subject': '230912CAA1120051', 'Category': 3},
    {'Subject': '230913CAA1120182', 'Category': 1},
    {'Subject': '230913CAA1120395', 'Category': 4},
    {'Subject': '230914CAA1120181', 'Category': 1},
    {'Subject': '230914CAA1120210', 'Category': 2},
    {'Subject': '230915CAA1120509', 'Category': 2},
    # {'Subject': '230915CAA1120509', 'Category': 2},
]

# 1 = callosal projecting neurons
# 2 = corticothalamic projecting neurons
# 3 = PV interneurons
# 4 = SST interneurons

# remove sub_## from the all_metrics_df['dataset'] column
all_metrics_df['dataset'] = all_metrics_df['dataset'].apply(lambda x: x.split('_')[1])

# now create a new column with the category by joining the two dataframes
all_metrics_df['Category'] = all_metrics_df['dataset'].map(
    pd.DataFrame(cleaned_data).set_index('Subject')['Category']
)

grouped = all_metrics_df.groupby('Category').apply(lambda x: pd.Series({
    'n_datasets': len(x),
    'n_is_cell': x['#_is_cell'].sum()
}))

#  print by specifying neuron type
print(grouped)
print(f"Callasal projecting neurons: {grouped.loc[0]}")
print(f"Corticothalamic projecting neurons: {grouped.loc[1]}")
print(f"SST interneurons: {grouped.loc[2]}")
