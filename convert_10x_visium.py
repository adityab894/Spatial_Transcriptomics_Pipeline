# "#!/usr/bin/env python3
"""
Convert 10x Genomics Visium data to pipeline format
Converts matrix.mtx + barcodes + features → gene expression CSV
Converts tissue_positions → coordinates CSV with pixel mapping
"""

import pandas as pd
import numpy as np
import scipy.io
import gzip
import json
import os
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt

def read_10x_mtx(mtx_file, barcodes_file, features_file):
    # \"\"\"Read 10x matrix format and convert to DataFrame\"\"\"
    print("📊 Reading 10x matrix files...")
    
    # Read matrix
    matrix = scipy.io.mmread(mtx_file)
    
    # Read barcodes (spot identifiers)
    with gzip.open(barcodes_file, 'rt') as f:
        barcodes = [line.strip() for line in f]
    
    # Read features (gene names)
    features_df = pd.read_csv(features_file, sep='	', header=None, compression='gzip')
    if features_df.shape[1] >= 2:
        gene_names = features_df.iloc[:, 1].tolist()  # Gene symbols (column 1)
    else:
        gene_names = features_df.iloc[:, 0].tolist()  # Gene IDs (column 0)
    
    # Convert to dense DataFrame (genes as columns, barcodes as rows)
    print(f"   Matrix shape: {matrix.shape}")
    print(f"   Barcodes: {len(barcodes)}")
    print(f"   Genes: {len(gene_names)}")
    
    # Matrix is usually genes x cells, we want cells x genes
    if matrix.shape[0] == len(gene_names):
        matrix = matrix.T  # Transpose
    
    df = pd.DataFrame(matrix.toarray(), index=barcodes, columns=gene_names)
    
    print(f"✅ Created expression matrix: {df.shape[0]} spots, {df.shape[1]} genes")
    return df

def process_tissue_positions(positions_file, scalefactors_file=None):
    # \"\"\"Process tissue positions and add pixel coordinates\"\"\"
    print("📍 Processing tissue positions...")
    
    # Read positions
    positions_df = pd.read_csv(positions_file, compression='gzip', header=None)
    
    # Standard 10x Visium tissue_positions format:
    # Column 0: barcode
    # Column 1: in_tissue (1 if spot is in tissue, 0 if not)
    # Column 2: array_row
    # Column 3: array_col  
    # Column 4: pxl_col_in_fullres (pixel x coordinate)
    # Column 5: pxl_row_in_fullres (pixel y coordinate)
    
    columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
    positions_df.columns = columns[:positions_df.shape[1]]
    
    # Filter for spots in tissue
    if 'in_tissue' in positions_df.columns:
        tissue_spots = positions_df[positions_df['in_tissue'] == 1].copy()
        print(f"   Spots in tissue: {len(tissue_spots)}")
    else:
        tissue_spots = positions_df.copy()
        print(f"   Total spots: {len(tissue_spots)}")
    
    # Create coordinates in the format expected by pipeline
    coords_df = pd.DataFrame({
        'spot': tissue_spots['barcode'],
        'x': tissue_spots['array_col'],
        'y': tissue_spots['array_row'],
        'new_x': tissue_spots['array_col'],  # Keep same for now
        'new_y': tissue_spots['array_row'],
        'pixel_x': tissue_spots['pxl_col_in_fullres'] if 'pxl_col_in_fullres' in tissue_spots.columns else tissue_spots['array_col'] * 100,
        'pixel_y': tissue_spots['pxl_row_in_fullres'] if 'pxl_row_in_fullres' in tissue_spots.columns else tissue_spots['array_row'] * 100
    })
    
    # Handle scale factors if available
    if scalefactors_file and os.path.exists(scalefactors_file):
        try:
            with gzip.open(scalefactors_file, 'rt') as f:
                scalefactors = json.load(f)
            
            # Apply high-resolution scaling if needed
            if 'tissue_hires_scalef' in scalefactors:
                scale_factor = scalefactors['tissue_hires_scalef']
                coords_df['pixel_x'] = coords_df['pixel_x'] * scale_factor
                coords_df['pixel_y'] = coords_df['pixel_y'] * scale_factor
                print(f"   Applied scale factor: {scale_factor}")
        except Exception as e:
            print(f"   Warning: Could not process scale factors: {e}")
    
    print(f"✅ Created coordinates: {len(coords_df)} spots with pixel coordinates")
    return coords_df

def convert_image(image_file, output_path):
    # \"\"\"Convert and extract tissue image\"\"\"
    print("🖼️  Processing tissue image...")
    
    # Open compressed image
    with gzip.open(image_file, 'rb') as f:
        img = Image.open(f)
        
        # Convert to RGB if needed
        if img.mode != 'RGB':
            img = img.convert('RGB')
        
        # Save as JPG
        img.save(output_path, 'JPEG', quality=95)
        
        print(f"✅ Saved tissue image: {img.size[0]}x{img.size[1]} pixels")
        return img.size

def convert_10x_visium_sample(input_dir, output_dir, sample_name):
    # \"\"\"Convert a 10x Visium sample to pipeline format\"\"\"
    
    print(f"🔄 Converting 10x Visium sample: {sample_name}")
    print("=" * 50)
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Expected 10x files
    files = {
        'matrix': f"{sample_name}_matrix.mtx.gz",
        'barcodes': f"{sample_name}_barcodes.tsv.gz", 
        'features': f"{sample_name}_features.tsv.gz",
        'positions': f"{sample_name}_P6_rep2_tissue_positions_list.csv.gz",
        'image': f"{sample_name}_P6_rep2_tissue_hires_image.png.gz",
        'scalefactors': f"{sample_name}_P6_rep2_scalefactors_json.json.gz"
    }
    
    # Check file existence
    missing_files = []
    for file_type, filename in files.items():
        file_path = input_path / filename
        if not file_path.exists():
            missing_files.append(filename)
            print(f"❌ Missing: {filename}")
        else:
            print(f"✅ Found: {filename}")
    
    if missing_files:
        print(f"Missing files: {missing_files}")
        return False
    
    try:
        # 1. Convert matrix to gene expression CSV
        print("1. Converting gene expression matrix...")
        expression_df = read_10x_mtx(
            input_path / files['matrix'],
            input_path / files['barcodes'], 
            input_path / files['features']
        )
        
        # Save as CSV (similar to cancer data format)
        expression_output = output_path / f"{sample_name}_stdata.csv"
        expression_df.to_csv(expression_output)
        print(f"✅ Saved: {expression_output}")
        
        # 2. Process coordinates
        print("2. Converting tissue positions...")
        coords_df = process_tissue_positions(
            input_path / files['positions'],
            input_path / files['scalefactors']
        )
        
        # Save coordinates
        coords_output = output_path / f"{sample_name}_coordinates.csv"
        coords_df.to_csv(coords_output, index=False)
        print(f"✅ Saved: {coords_output}")
        
        # 3. Convert tissue image
        print("3. Converting tissue image...")
        image_output = output_path / f"{sample_name}_tissue_image.jpg"
        image_size = convert_image(input_path / files['image'], image_output)
        print(f"✅ Saved: {image_output}")
        
        # 4. Create summary
        print("=" * 50)
        print("📊 CONVERSION SUMMARY")
        print("=" * 50)
        print(f"Sample: {sample_name}")
        print(f"Expression data: {expression_df.shape[0]} spots, {expression_df.shape[1]} genes")
        print(f"Coordinates: {len(coords_df)} spots with pixel coordinates")
        print(f"Image: {image_size[0]}x{image_size[1]} pixels")
        print(f"Output directory: {output_path}")
        
        print("📁 Files created: ")
        print(f"  - {sample_name}_stdata.csv")
        print(f"  - {sample_name}_coordinates.csv") 
        print(f"  - {sample_name}_tissue_image.jpg")
        
        print(f"🚀 Ready for pipeline! Create directory structure:")
        print(f"data/{sample_name}/")
        print(f"├── {sample_name}_stdata.csv")
        print(f"├── {sample_name}_coordinates.csv")
        print(f"└── {sample_name}_tissue_image.jpg")
        
        return True
        
    except Exception as e:
        print(f"❌ Conversion failed: {e}")
        return False

def main():
    # \"\"\"Main conversion function\"\"\"
    print("🧬 10x Genomics Visium to Pipeline Converter")
    print("=" * 60)
    
    # Example usage
    sample_name = "GSM4565826"  # Change this to your sample name
    input_directory = "./10x_data"  # Directory containing your 10x files
    output_directory = "./data/normal_sample"  # Output for pipeline
    
    print(f"Sample: {sample_name}")
    print(f"Input: {input_directory}")
    print(f"Output: {output_directory}")
    
    # Run conversion
    success = convert_10x_visium_sample(input_directory, output_directory, sample_name)
    
    if success:
        print("🎉 Conversion completed successfully!")
        print("You can now run the spatial transcriptomics pipeline:")
        print("python spatial_transcriptomics_pipeline.py --sample ./data/normal_sample")
    else:
        print("❌ Conversion failed. Please check the files and try again.")

if __name__ == "__main__":
    main()