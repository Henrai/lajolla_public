import subprocess
import os

# # Define the range of frames
# i = 0 # Start frame
# j = 528  # End frame

# # Define the paths
# scene_dir = "../scenes/animation/scene/"
# output_dir = "./out/"
# os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

# # Run the command for each frame
# for frame in range(i, j + 1):
#     frame_str = f"{frame:04d}"  # Format as 4-digit number
#     scene_file = f"{scene_dir}scene_{frame_str}.xml"
#     output_file = f"{output_dir}{frame_str}.exr"

#     command = ["./lajolla", scene_file, "-o", output_file]
    
#     print(f"Running: {' '.join(command)}")
#     subprocess.run(command, check=True)  # Run the command

# print("Rendering completed!")


scene_dir = "../scenes/volpath_test/"
output_dir = "./res/"

file_name = [
    # "volpath_test1",
    # "volpath_test2",
    # "volpath_test3",
    # "volpath_test4",
    # "volpath_test4_2",
    # "volpath_test5",
    # "volpath_test5_2",
    # "volpath_test6",
    # "vol_cbox",
    # "vol_cbox_teapot",
    "hetvol",
    "hetvol_colored",
]

for file in file_name:
    scene_file = f"{scene_dir}{file}.xml"
    output_file = f"{output_dir}{file}.exr"
    command = ["./lajolla", scene_file, "-o", output_file]
    print(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True)  # Run the command

