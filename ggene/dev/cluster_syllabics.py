#!/usr/bin/env python3
"""
Automatically cluster Canadian Aboriginal Syllabics by visual similarity.
Uses image analysis to group characters by base shape and identify diacritical variants.
"""

from PIL import Image
import numpy as np
from pathlib import Path
from collections import defaultdict
import json


def load_image_as_binary(image_path, threshold=128):
    """Load image and convert to binary array."""
    img = Image.open(image_path).convert('L')  # Convert to grayscale
    arr = np.array(img)

    # Threshold to binary (1 = black/character, 0 = white/background)
    binary = (arr < threshold).astype(int)

    return binary


def get_vertical_extent(binary_img):
    """Calculate the vertical extent of the character (as fraction of image height)."""
    # Find rows that contain any black pixels (character pixels)
    rows_with_char = np.any(binary_img == 1, axis=1)

    if not np.any(rows_with_char):
        return 0.0  # Empty image

    # Find first and last row with character
    char_rows = np.where(rows_with_char)[0]
    top = char_rows[0]
    bottom = char_rows[-1]

    height = bottom - top + 1
    total_height = binary_img.shape[0]

    return height / total_height


def count_whitespace(binary_img):
    """Count the number of white pixels (background) in the image."""
    return np.sum(binary_img == 1)


def normalize_image(binary_img):
    """Normalize binary image for comparison (flatten to 1D)."""
    return binary_img.flatten()


def rotate_90(binary_img):
    """Rotate image 90 degrees clockwise."""
    return np.rot90(binary_img, k=-1)


def rotate_180(binary_img):
    """Rotate image 180 degrees."""
    return np.rot90(binary_img, k=2)


def rotate_270(binary_img):
    """Rotate image 270 degrees clockwise (90 counter-clockwise)."""
    return np.rot90(binary_img, k=1)


def flip_horizontal(binary_img):
    """Flip image horizontally."""
    return np.fliplr(binary_img)


def flip_vertical(binary_img):
    """Flip image vertically."""
    return np.flipud(binary_img)


def get_all_transforms(binary_img):
    """Generate all rotations and flips of the image."""
    transforms = []

    # Original and rotations
    transforms.append(("orig", binary_img))

    # Horizontal flips of each rotation
    transforms.append(("flip_h", flip_horizontal(binary_img)))
    transforms.append(("flip_v", flip_vertical(binary_img)))
    transforms.append(("flip_hv", flip_vertical(flip_horizontal(binary_img))))

    return transforms


def similarity(img1, img2):
    """Calculate similarity between two binary images using dot product.

    Returns value between -1 and 1, where 1 means identical.
    We need to normalize: treat 0 (black) as -1 and 1 (white) as +1.
    """
    # Convert 0/1 to -1/+1 for proper correlation
    img1_normalized = img1 * 2 - 1
    img2_normalized = img2 * 2 - 1

    # Flatten images
    vec1 = img1_normalized.flatten()
    vec2 = img2_normalized.flatten()

    # Compute dot product normalized by length
    dot = np.dot(vec1, vec2)
    norm = len(vec1)

    return dot / norm


def find_best_match(binary_img, cluster_keys, threshold=0.7):
    """Find the best matching cluster key for a character.

    Returns (cluster_key_char, best_similarity, best_transform) or (None, 0, None).
    """
    best_match = None
    best_sim = -1
    best_transform = None

    # Get all transforms of the input image
    transforms = get_all_transforms(binary_img)

    for key_char, key_img in cluster_keys.items():
        for transform_name, transformed_img in transforms:
            sim = similarity(key_img, transformed_img)

            if sim > best_sim:
                best_sim = sim
                best_match = key_char
                best_transform = transform_name

    if best_sim >= threshold:
        return best_match, best_sim, best_transform, transformed_img
    else:
        return None, best_sim, None, None


def cluster_characters(image_dir, similarity_threshold=0.7, diacritic_height_threshold=0.5, diacritic_mass_threshold = 0.05):
    """Main clustering function.

    Returns:
        - full_character_clusters: dict mapping cluster key char -> list of chars in cluster
        - diacritics: list of characters that are small (likely diacritics/finals)
        - cluster_keys_imgs: dict mapping cluster key char -> binary image
    """
    image_path = Path(image_dir)

    # Load all character images
    char_images = {}

    for img_file in sorted(image_path.glob("U+*.png")):
        # Extract hex code and character
        hex_code = img_file.stem  # e.g., "U+1403"
        codepoint = int(hex_code.replace("U+", ""), 16)
        char = chr(codepoint)

        # Load as binary
        binary_img = load_image_as_binary(img_file)
        char_images[char] = binary_img

    print(f"Loaded {len(char_images)} character images")

    # Separate diacritics from full characters
    full_characters = {}
    diacritics = []

    for char, binary_img in char_images.items():
        extent = get_vertical_extent(binary_img)
        mass = np.sum(binary_img) / binary_img.shape[0] / binary_img.shape[1]

        if extent < diacritic_height_threshold or mass < diacritic_mass_threshold:
            diacritics.append(char)
            print(f"  Diacritic/final: {char} (U+{ord(char):04X}) - extent: {extent:.2f}")
        else:
            full_characters[char] = binary_img

    print(f"\nFound {len(full_characters)} full characters")
    print(f"Found {len(diacritics)} diacritics/finals\n")

    # Cluster full characters
    clusters = defaultdict(list)  # cluster_key -> list of chars
    cluster_keys_imgs = {}  # cluster_key -> binary image

    # Sort characters by whitespace (most whitespace first = fewest diacritics)
    chars_by_whitespace = sorted(
        full_characters.items(),
        key=lambda x: count_whitespace(x[1]),
        reverse=True
    )

    for char, binary_img in chars_by_whitespace:
        # Try to find a matching cluster
        match, sim, transform, transformed_img = find_best_match(binary_img, cluster_keys_imgs, similarity_threshold)

        if match is not None:
            # Add to existing cluster
            clusters[match].append(char)
            update_key_img(match, cluster_keys_imgs, clusters, transformed_img)
            print(f"  {char} (U+{ord(char):04X}) -> cluster {match} (sim: {sim:.3f}, transform: {transform})")
        else:
            # Create new cluster with this character as the key
            clusters[char].append(char)
            cluster_keys_imgs[char] = binary_img
            print(f"  {char} (U+{ord(char):04X}) -> NEW CLUSTER (key)")

    return clusters, diacritics, cluster_keys_imgs

def update_key_img(key, cluster_keys_imgs, clusters, new_img):
    
    num_chars = len(clusters.get(key, []))
    key_img = num_chars * cluster_keys_imgs[key]
    new_key_img = (key_img + new_img) / (num_chars + 1)
    
    cluster_keys_imgs[key] = new_key_img
    
    return cluster_keys_imgs
    
def save_key_imgs(cluster_keys_imgs):
    
    outpath = Path("./syllabics_keys")
    outpath.mkdir(exist_ok = True)
    
    for char, img in cluster_keys_imgs.items():
        
        char_path = outpath / f"U+{hex(ord(char))[2:]}_key.png"
        
        immin = np.min(img)
        immax = np.max(img)
        img_norm = (255*(img - immin) / (immax - immin)).astype("uint8")
        
        pil_img = Image.fromarray(img_norm, mode="L")
        
        pil_img.save(char_path)


def save_results(clusters, diacritics, output_file="syllabics_clusters.json"):
    """Save clustering results to JSON file."""

    # Convert to serializable format
    output = {
        "clusters": {},
        "diacritics": diacritics
    }

    for key, members in clusters.items():
        output["clusters"][key] = members

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"\nResults saved to {output_file}")


def print_summary(clusters, diacritics):
    """Print a summary of the clustering results."""
    print("\n" + "="*60)
    print("CLUSTERING SUMMARY")
    print("="*60)

    print(f"\nTotal clusters: {len(clusters)}")
    print(f"Total diacritics/finals: {len(diacritics)}")

    print("\nCluster sizes:")
    cluster_sizes = defaultdict(int)
    for members in clusters.values():
        cluster_sizes[len(members)] += 1

    for size in sorted(cluster_sizes.keys()):
        count = cluster_sizes[size]
        print(f"  {count} clusters with {size} members")

    print("\nLargest clusters:")
    sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (key, members) in enumerate(sorted_clusters[:10]):
        print(f"  {i+1}. {key} (U+{ord(key):04X}): {len(members)} members")
        print(f"     Members: {''.join(members)}")


def main():
    """Main function."""
    print("Canadian Aboriginal Syllabics Clustering")
    print("=" * 60)

    # Configuration
    image_dir = "./ggene/dev/syllabics_rendered"
    similarity_threshold = 0.80  # Adjust this to tune clustering (higher = stricter)
    diacritic_height_threshold = 0.27  # Characters smaller than this are diacritics
    diacritic_mass_threshold = 0.043

    # Cluster characters
    clusters, diacritics, cluster_keys = cluster_characters(
        image_dir,
        similarity_threshold=similarity_threshold,
        diacritic_height_threshold=diacritic_height_threshold,
        diacritic_mass_threshold = diacritic_mass_threshold
    )

    save_key_imgs(cluster_keys)

    # Print summary
    print_summary(clusters, diacritics)

    # Save results
    save_results(clusters, diacritics)

    

    print("\n✓ Clustering complete!")


if __name__ == "__main__":
    main()
