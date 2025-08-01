import os
import re

base_dir = "/data/a.saricaoglu/repo/COMPAS/Runs/Default_WD_Enabled_Detailed"
start_num = 52  # The next available number

# Find all folders matching COMPAS_Output_X(1)
pattern = re.compile(r"COMPAS_Output_(\d+) \(1\)$")
print(f"Searching in directory: {base_dir}")
print(f"Using pattern: {pattern.pattern}")
folders = [f for f in os.listdir(base_dir) if pattern.match(f)]
print(f"Found folders: {folders}")
# Sort for reproducibility (optional)
folders.sort(key=lambda x: int(pattern.match(x).group(1)))
print(f"Found folders: {folders}")
for i, folder in enumerate(folders):
    new_num = start_num + i
    new_name = f"COMPAS_Output_{new_num}"
    src = os.path.join(base_dir, folder)
    dst = os.path.join(base_dir, new_name)
    print(f"Renaming {src} -> {dst}")
    os.rename(src, dst)