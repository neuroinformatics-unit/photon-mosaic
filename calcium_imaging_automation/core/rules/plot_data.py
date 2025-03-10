from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from snakemake.script import snakemake
from scipy.stats import ttest_rel
import itertools


print("Plotting data...")



window_size_stats = 15
window_size_plot = 60



# datasets = snakemake.params.datasets
f_path = Path(snakemake.input[0])

print(f"Processing file: {f_path}")

dataset_path = f_path.parent.parent.parent.parent
dataset = dataset_path.name
f_neu_path = f_path.parent / "Fneu.npy"
is_cell_path = f_path.parent / "iscell.npy"
derotated_full_csv_path = (
    dataset_path / "ses-0" / "funcimg" / "derotation" / "derotated_full.csv"
)
saving_path = Path(dataset_path) / "ses-0" / "traces"



#  delete pre-existing files
for item in saving_path.glob("*"):
    if item.is_file():
        item.unlink()
    if item.is_dir():
        for file in item.glob("*"):
            file.unlink()
        item.rmdir()

saving_path.mkdir(exist_ok=True)
print(f"Dataset path: {dataset_path}")
print(f"Dataset: {dataset}")
print(f"Derotated full csv path: {derotated_full_csv_path}")
print(f"Fneu path: {f_neu_path}")
print(f"Saving path: {saving_path}")


f = np.load(f_path)
fneu = np.load(f_neu_path)
f_corrected = f - 0.7 * fneu

#  filter by is_cell
classifications = np.load(is_cell_path)
is_cell = classifications[:, 0]
rois_selection = np.where(is_cell)[0] 
f_corrected = f_corrected[rois_selection]
f_corrected = f_corrected.T


rotation_info = pd.read_csv(derotated_full_csv_path)
print("Pivot table of rotation info")
print(rotation_info.pivot_table(values='rotation_count', index='speed', columns='direction', aggfunc=lambda x: len(x.unique())))

df_fluorescence = pd.DataFrame(f_corrected, columns=[i for i in range(f_corrected.shape[1])])

roi_n = df_fluorescence.shape[1]

# join df_fluorescence with rotation_info
df_fluorescence = df_fluorescence.join(rotation_info)

#  discard column called "frame" and "clock"
df_fluorescence = df_fluorescence.loc[:, ~df_fluorescence.columns.isin(["clock"])]

#  melt columns that have a number in the name
df_fluorescence = df_fluorescence.melt(
    id_vars=["rotation_angle", "direction", "speed", "rotation_count", "frame"],
    var_name="roi_n", value_name="roi_fluorescence")


frame_rate = 6.74
time = np.linspace(0, len(f_corrected) / frame_rate, len(f_corrected))

#%%
data = pd.DataFrame(
    columns=["time", "frame", "roi_n", "roi_fluorescence", "speed", "direction", "repetition", "is_start", "is_end"]
)

#  repeat time as many times as there are rois
data["time"] = np.repeat(time, roi_n)
data["frame"] = df_fluorescence["frame"]
data["roi_n"] = [num for num in rois_selection for _ in range(len(time))] 
data["roi_fluorescence"] = df_fluorescence["roi_fluorescence"]

#  set start and end to all false
data["is_start"] = False
data["is_end"] = False

#%%

#  where is rotation happening?
filter = (rotation_info["rotation_count"] > 0) | (rotation_info["rotation_count"] < 0)
rotation_info = rotation_info[filter]
starts = rotation_info.groupby("rotation_count").first()
ends = rotation_info.groupby("rotation_count").last()

data.loc[data["frame"].isin(starts["frame"]), "is_start"] = True
data.loc[data["frame"].isin(ends["frame"]), "is_end"] = True
data["speed"] = df_fluorescence["speed"]
data["direction"] = df_fluorescence["direction"]
data["repetition"] = df_fluorescence["rotation_count"]

# %%


# Initialize the new column with NaNs
data["mean"] = np.nan

# Compute means for 'is_start' and 'is_end' in one line each
data.loc[data["is_start"] == 1, "mean"] = data.index[
    data["is_start"] == 1
    ].map(
        lambda i: data.loc[max(0, i-window_size_plot):i-1, "roi_fluorescence"].mean()
    )
data.loc[data["is_end"] == 1, "mean"] = data.index[
    data["is_end"] == 1
    ].map(
        lambda i: data.loc[i+1:i+window_size_plot, "roi_fluorescence"].mean()
    )


# Function to extract fluorescence before 'is_start' and after 'is_end'
def extract_windows(df):
    before = df.loc[df["is_start"] == 1, "mean"]
    after = df.loc[df["is_end"] == 1, "mean"]
    return before, after

# Perform t-test within each group
def compute_ttest(group):
    before, after = extract_windows(group)

    if len(before) > 1 and len(after) > 1:  # Ensure there are enough samples
        t_stat, p_value = ttest_rel(after, before, alternative="greater")
    else:
        t_stat, p_value = np.nan, np.nan  # Not enough samples for test
    return pd.Series({"t_stat": t_stat, "p_value": p_value, "significant": p_value < 0.05})

#%%
# Group by 'roi_n', 'speed', 'direction' and apply t-test
results_all_groups = data.groupby(["roi_n", "speed", "direction"]).apply(compute_ttest)
results_cw_vs_ccw = data.groupby(["roi_n", "direction"]).apply(compute_ttest)

print(results_all_groups[results_all_groups["significant"] == True])
print(results_cw_vs_ccw[results_cw_vs_ccw["significant"] == True])

results_all_groups.to_csv(saving_path / "paired_t_test.csv")
results_cw_vs_ccw.to_csv(saving_path / "paired_t_test_cw_vs_ccw.csv")

# %%
# Now let's plot the results
for roi in data.roi_n.unique():
    print(f"Plotting ROI {roi}")
    fig, ax = plt.subplots(2, 4, figsize=(20, 10))
    
    product = itertools.product(data.direction.unique(), sorted(data.speed.unique()))
    product = [i for i in product if not np.isnan(i[0]) and not np.isnan(i[1])]
    for i, (direction, speed) in enumerate(product):
        
        ax[i // 4, i % 4].set_title(f"Speed: {speed}, Direction: {'CW' if direction == 1 else 'CCW'}")
        ax[i // 4, i % 4].set_xlabel("Frame")
        ax[i // 4, i % 4].set_ylabel("Fluorescence")

        subset = data[
            (data["roi_n"] == roi) &
            (data["speed"] == speed) &
            (data["direction"] == direction)
        ]
        single_reps = []
        for rep in subset[subset["is_start"] == 1].repetition.unique():
            
            idx_start = subset.loc[
                (subset["repetition"] == rep) & (subset["is_start"] == 1), "frame"
            ].index.item() - window_size_plot
            idx_end = subset.loc[
                (subset["repetition"] == rep) & (subset["is_end"] == 1), "frame"
            ].index.item() + window_size_plot
            
            ax[i // 4, i % 4].plot(
                data.iloc[idx_start:idx_end]["roi_fluorescence"].values,
                color="gray",
            )
            single_reps.append(data.iloc[idx_start:idx_end]["roi_fluorescence"].values.tolist())

        #  single_reps have different lengths, so we need to pad them with NaNs
        max_len = max(len(i) for i in single_reps)
        single_reps = [i + [np.nan] * (max_len - len(i)) for i in single_reps]
        single_reps = np.array(single_reps)

        #  also plot the mean
        filter = (data["roi_n"] == roi) & (data["speed"] == speed) & (data["direction"] == direction)
        ax[i // 4, i % 4].plot(
            np.mean(single_reps, axis=0),
            label="Mean",
            color="red" if results_all_groups.loc[(roi, speed, direction), "significant"] else "blue"
        )

        #  draw vertical lines at start and end based on window_size
        ax[i // 4, i % 4].axvline(window_size_plot, color="black", linestyle="--")
        ax[i // 4, i % 4].axvline(max_len - window_size_plot, color="black", linestyle="--")

        #  remove top and right spines
        ax[i // 4, i % 4].spines["top"].set_visible(False)
        ax[i // 4, i % 4].spines["right"].set_visible(False)

        #  write p-value if significant
        p_value = results_all_groups.loc[(roi, speed, direction), "p_value"]
        if p_value < 0.05:
            ax[i // 4, i % 4].text(0.8, 0.9, f"p-value: {p_value:.3f}", transform
                =ax[i // 4, i % 4].transAxes, color="red")
            

        
    plt.suptitle(f"ROI {roi}")
    plt.tight_layout()

    plt.savefig(saving_path / f"roi_{roi}.png")
    plt.close(fig)

