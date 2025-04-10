#!/bin/bash

source /home/dot/miniforge3/etc/profile.d/conda.sh
conda activate bandage || { echo "Failed to activate conda environment 'bandage'"; exit 1; }

# Define the parent directory containing the assembler directories as the current directory containin
current_dir=$(pwd)
parent_dir=${current_dir%/bin}
assembly_dir="$parent_dir/assemblies"


# Output directory for concatenated images
output_dir="$parent_dir/bin/summaries/bandage_plots"
mkdir -p "$output_dir"

# Ensure font exists
font_path="/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
if [ ! -f "$font_path" ]; then
    echo "Error: Font not found at $font_path"
    exit 1
fi


# Step 1: Generate Bandage plots for .gfa files directly inside `assemblies/*/`
find "$assembly_dir" -maxdepth 2 -mindepth 2 -type f -name "*.gfa" | while read -r gfa_file; do
    base_name=$(basename "$gfa_file" .gfa)
    output_png="$(dirname "$gfa_file")/${base_name}.png"

    # Generate Bandage image
    Bandage image "$gfa_file" "$output_png" --height 800 --width 800
    echo "Generated: $output_png"
done

# Step 2: Generate Bandage plots for plasmid GFA files in hybracter_long_01/plasmids and hybracter_hybrid/plasmids
for plasmid_dir in hybracter_long_01/plasmids hybracter_hybrid/plasmids; do
    full_path="$assembly_dir/$plasmid_dir"
    if [ -d "$full_path" ]; then
        find "$full_path" -type f -name "*.gfa" | while read -r gfa_file; do
            base_name=$(basename "$gfa_file" .gfa)
            output_png="$(dirname "$gfa_file")/${base_name}.png"

            Bandage image "$gfa_file" "$output_png" --height 800 --width 800
            echo "Generated: $output_png"
        done
    fi
done



# Step 3: Define Lables for assembly images
declare -A folder_labels=(
    ["hybracter_long_01"]="Hybracter (long)"
    ["flye_01"]="Flye"
    ["autocycler_consensus"]="Autocycler"
    ["hybracter_hybrid"]="Hybracter (hybrid)"
    ["unicycler"]="Unicycler"
    ["unicycler_bold"]="Unicycler (bold)"
)


# Step 4: Define labels for plasmid images
declare -A plasmid_labels=(
    ["hybracter_long_01"]="Hybracter (long) plasmids"
    ["hybracter_hybrid"]="Hybracter (hybrid) plasmids"
)



# Iterate over sample names from 'flye_01' directory
for flye_img in "$assembly_dir/flye_01/"*.png; do
    sample=$(basename "$flye_img" .png)

    # Collect images from the required directories
    images=()
    for folder in hybracter_long_01 flye_01 autocycler_consensus hybracter_hybrid unicycler unicycler_bold; do
        img_path="$assembly_dir/$folder/${sample}.png"
        if [[ -f "$img_path" ]]; then
            images+=("$img_path")
        else
            echo "Warning: Missing image $img_path"
            images+=("xc:white")  # Placeholder blank image
        fi
    done

    # Generate header labels for each image
    header_images=()
    for folder in hybracter_long_01 flye_01 autocycler_consensus hybracter_hybrid unicycler unicycler_bold; do
        label_img="/tmp/label_${folder}.png"
        magick -size 800x100 xc:white -gravity center -font "$font_path" -pointsize 50 -annotate +0+0 "${folder_labels[$folder]}" "$label_img"
        header_images+=("$label_img")
    done

    # Stack headers above images
    output_image="$output_dir/${sample}_row.png"
    magick \( "${header_images[@]}" +append \) \( "${images[@]}" +append \) -append "$output_image"
    echo "Created: $output_image"

    # Step 6: Collect plasmid images for this sample
    plasmid_images=()
    for folder in hybracter_long_01 hybracter_hybrid; do
        plasmid_img="$assembly_dir/$folder/plasmids/${sample}_plasmids.png"
        if [[ -f "$plasmid_img" ]]; then
            plasmid_images+=("$plasmid_img")
        else
            echo "Warning: Missing plasmid image $plasmid_img"
            plasmid_images+=("xc:white")  # Placeholder blank image
        fi
    done

    # Generate labels for plasmid images
    plasmid_header_images=()
    for folder in hybracter_long_01 hybracter_hybrid; do
        label_img="/tmp/label_${folder}_plasmids.png"
        magick -size 800x100 xc:white -gravity center -font "$font_path" -pointsize 50 -annotate +0+0 "${plasmid_labels[$folder]}" "$label_img"
        plasmid_header_images+=("$label_img")
    done

    # Step 7: Append plasmid images to the right of the existing row
    output_image_with_plasmids="$output_dir/${sample}_with_plasmids.png"
    magick \( "$output_image" \( "${plasmid_header_images[@]}" +append \) \( "${plasmid_images[@]}" +append \) -append \) +append "$output_image_with_plasmids"
    echo "Created: $output_image_with_plasmids"





done



conda deactivate
