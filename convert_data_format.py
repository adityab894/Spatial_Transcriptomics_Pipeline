#!/usr/bin/env python3
"""
Data Format Converter for Spatial Transcriptomics Pipeline
Converts TSV files with GSM names to CSV format with expected naming convention
"""

import pandas as pd
import os
import glob
from pathlib import Path
import shutil
import argparse
from pathlib import Path
import warnings
from datetime import datetime
import logging
from typing import Dict, List, Tuple, Optional
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import pdist, squareform


def find_data_files(data_directory: str = "./data"):
    """Find all TSV/CSV files in sample directories"""
    print(f"🔍 Scanning directory: {data_directory}")
    
    sample_dirs = [d for d in Path(data_directory).iterdir() if d.is_dir()]
    found_files = {}
    
    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        print(f"\n📁 Checking sample: {sample_name}")
        
        # Look for various file patterns
        file_patterns = [
            "*.tsv", "*.csv", "*.txt", 
            "GSM*.tsv", "GSM*.csv", "GSM*.txt",
            "*_stdata.*", "*stdata.*", "*ST*", "*st*"
        ]
        
        sample_files = {
            'expression': [],
            'coordinates': [],
            'images': [],
            'other': []
        }
        
        for pattern in file_patterns:
            files = list(sample_dir.glob(pattern))
            for file_path in files:
                file_name = file_path.name.lower()
                
                # Categorize files
                if any(keyword in file_name for keyword in ['stdata', 'expression', 'count', 'matrix', 'gsm']):
                    if file_path.suffix.lower() in ['.tsv', '.csv', '.txt']:
                        sample_files['expression'].append(file_path)
                if any(keyword in file_name for keyword in ['coord', 'position', 'spatial', 'tissue_positions', 'spot_data']):
                    if file_path.suffix.lower() in ['.tsv', '.csv', '.txt']:
                        sample_files['coordinates'].append(file_path)
                elif any(keyword in file_name for keyword in ['stdata', 'expression', 'count', 'matrix']):
                    if file_path.suffix.lower() in ['.tsv', '.csv', '.txt']:
                        sample_files['expression'].append(file_path)
                elif file_path.suffix.lower() in ['.jpg', '.jpeg', '.png', '.tiff', '.tif']:
                    sample_files['images'].append(file_path)
                elif file_name.startswith('gsm') and file_path.suffix.lower() in ['.tsv', '.csv', '.txt']:
                    # GSM files that don't match specific patterns - likely expression data
                    sample_files['expression'].append(file_path)
                else:
                    sample_files['other'].append(file_path)
        
        # Remove duplicates
        for category in sample_files:
            sample_files[category] = list(set(sample_files[category]))
        
        found_files[sample_name] = sample_files
        
        # Print findings
        for category, files in sample_files.items():
            if files:
                print(f"  {category.capitalize()}: {len(files)} files")
                for file_path in files:
                    print(f"    - {file_path.name}")
    
    return found_files

def detect_file_format(file_path: Path) -> str:
    """Detect if file is TSV or CSV format"""
    try:
        # Read first few lines to detect separator
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            second_line = f.readline().strip()
        
        tab_count = first_line.count('\t')
        comma_count = first_line.count(',')
        
        if tab_count > comma_count and tab_count > 0:
            return 'tsv'
        elif comma_count > tab_count and comma_count > 0:
            return 'csv'
        else:
            # Try to load with pandas to detect
            try:
                df_tab = pd.read_csv(file_path, sep='\t', nrows=5)
                df_comma = pd.read_csv(file_path, sep=',', nrows=5)
                
                if df_tab.shape[1] > df_comma.shape[1]:
                    return 'tsv'
                else:
                    return 'csv'
            except:
                return 'unknown'
    except Exception as e:
        print(f"  ⚠️  Could not detect format for {file_path.name}: {e}")
        return 'unknown'

def convert_file_format(input_path: Path, output_path: Path, file_format: str) -> bool:
    """Convert file from TSV to CSV or verify CSV format"""
    try:
        print(f"  🔄 Converting: {input_path.name} → {output_path.name}")
        
        # Read file with appropriate separator
        if file_format == 'tsv':
            df = pd.read_csv(input_path, sep='\t', index_col=0)
        elif file_format == 'csv':
            df = pd.read_csv(input_path, index_col=0)
        else:
            # Try both separators
            try:
                df = pd.read_csv(input_path, sep='\t', index_col=0)
                if df.shape[1] < 10:  # If too few columns, try comma
                    df = pd.read_csv(input_path, sep=',', index_col=0)
            except:
                df = pd.read_csv(input_path, sep=',', index_col=0)
        
        # Basic validation
        if df.shape[0] < 10 or df.shape[1] < 10:
            print(f"  ⚠️  Warning: Small dataset ({df.shape[0]} spots, {df.shape[1]} genes)")
        
        # Save as CSV
        df.to_csv(output_path)
        
        print(f"  ✅ Converted successfully: {df.shape[0]} spots, {df.shape[1]} genes")
        return True
        
    except Exception as e:
        print(f"  ❌ Error converting {input_path.name}: {e}")
        return False

def standardize_filenames(sample_dir: Path, files_info: dict) -> dict:
    """Standardize filenames to match pipeline expectations"""
    conversions = {}
    sample_name = sample_dir.name
    
    print(f"\n🏷️  Standardizing filenames for: {sample_name}")
    
    # Handle expression files
    if files_info['expression']:
        for i, expr_file in enumerate(files_info['expression']):
            # Generate standard name
            if i == 0:
                new_name = f"{sample_name}_stdata.csv"
            else:
                new_name = f"{sample_name}_stdata_{i+1}.csv"
            
            new_path = sample_dir / new_name
            conversions[expr_file] = new_path
    
    # Handle coordinate files
    if files_info['coordinates']:
        for i, coord_file in enumerate(files_info['coordinates']):
            if i == 0:
                new_name = f"{sample_name}_coordinates.csv"
            else:
                new_name = f"{sample_name}_coordinates_{i+1}.csv"
            
            new_path = sample_dir / new_name
            conversions[coord_file] = new_path
    
    # Handle image files
    if files_info['images']:
        for i, img_file in enumerate(files_info['images']):
            extension = img_file.suffix.lower()
            if i == 0:
                new_name = f"{sample_name}_tissue_image{extension}"
            else:
                new_name = f"{sample_name}_tissue_image_{i+1}{extension}"
            
            new_path = sample_dir / new_name
            
            # For images, just copy/rename (no format conversion needed)
            if img_file != new_path:
                conversions[img_file] = new_path
    
    return conversions

def cleanup_remaining_gsm_files(sample_dir: Path, converted_files: List[str]):
    # \"\"\"Remove any remaining GSM files after successful conversion\"\"\"
    try:
        # Find any remaining GSM files
        remaining_gsm_files = []
        for file_path in sample_dir.iterdir():
            if file_path.is_file() and file_path.name.startswith('GSM'):
                # Don't remove if it's actually one of our converted files
                if file_path.name not in converted_files:
                    remaining_gsm_files.append(file_path)
        
        if remaining_gsm_files:
            print(f"  🧹 Cleaning up {len(remaining_gsm_files)} remaining original files:")
            for file_path in remaining_gsm_files:
                print(f"    🗑️  {file_path.name}")
                file_path.unlink()
    except Exception as e:
        print(f"  ⚠️  Warning: Could not cleanup remaining files: {e}")

def fix_expression_file_naming(sample_dir: Path, sample_name: str):
    # \"\"\"Fix incorrectly named expression files from *_stdata_2.csv to *_stdata.csv\"\"\"
    try:
        wrong_name = sample_dir / f"{sample_name}_stdata_2.csv"
        correct_name = sample_dir / f"{sample_name}_stdata.csv"
        
        if wrong_name.exists():
            if correct_name.exists():
                print(f"  🔧 Removing duplicate incorrect expression file")
                wrong_name.unlink()
            else:
                print(f"  🔧 Fixing expression file name: {wrong_name.name} → {correct_name.name}")
                wrong_name.rename(correct_name)
    except Exception as e:
        print(f"  ⚠️  Warning: Could not fix expression file naming: {e}")

def backup_original_files(data_directory: str):
    """Create backup of original data"""
    backup_dir = Path(data_directory + "_backup")
    
    if backup_dir.exists():
        print(f"⚠️  Backup already exists: {backup_dir}")
        response = input("Do you want to overwrite it? (y/n): ").strip().lower()
        if response != 'y':
            print("📁 Using existing backup")
            return str(backup_dir)
        else:
            shutil.rmtree(backup_dir)
    
    print(f"💾 Creating backup: {backup_dir}")
    shutil.copytree(data_directory, backup_dir)
    print("✅ Backup created successfully")
    
    return str(backup_dir)

def convert_all_files(data_directory: str = "./data", create_backup: bool = True):
    """Main conversion function"""
    print("🔄 Data Format Converter for Spatial Transcriptomics")
    print("=" * 60)
    
    # Check if data directory exists
    if not os.path.exists(data_directory):
        print(f"❌ Data directory not found: {data_directory}")
        return False
    
    # Create backup if requested
    if create_backup:
        backup_path = backup_original_files(data_directory)
        print(f"📁 Original data backed up to: {backup_path}")
    
    # Find all data files
    found_files = find_data_files(data_directory)
    
    if not found_files:
        print("❌ No sample directories found")
        return False
    
    # Process each sample
    conversion_summary = {
        'successful': 0,
        'failed': 0,
        'details': {}
    }
    
    for sample_name, files_info in found_files.items():
        print(f"\n🔬 Processing sample: {sample_name}")
        sample_dir = Path(data_directory) / sample_name
        
        # Check if we have expression data
        if not files_info['expression']:
            print(f"  ⚠️  No expression data files found in {sample_name}")
            conversion_summary['failed'] += 1
            conversion_summary['details'][sample_name] = 'No expression files found'
            continue
        
        # Generate standardized filenames
        conversions = standardize_filenames(sample_dir, files_info)
        
        sample_success = True
        converted_files = []
        
        # Process each file conversion
        for original_path, new_path in conversions.items():
            if original_path.suffix.lower() in ['.jpg', '.jpeg', '.png', '.tiff', '.tif']:
                # Handle image files (just rename/copy)
                if original_path != new_path:
                    try:
                        if new_path.exists():
                            new_path.unlink()  # Remove if exists
                        shutil.copy2(original_path, new_path)
                        print(f"  📷 Renamed image: {original_path.name} → {new_path.name}")
                        converted_files.append(new_path.name)
                        # Remove original image file after successful copy
                        original_path.unlink()
                        print(f"  🗑️  Removed original: {original_path.name}")
                    except Exception as e:
                        print(f"  ❌ Error renaming image {original_path.name}: {e}")
                        sample_success = False
                else:
                    # If names are the same, no conversion needed
                    print(f"  📷 Image already has correct name: {original_path.name}")
                    converted_files.append(original_path.name)
            else:
                # Handle data files (TSV/CSV conversion)
                file_format = detect_file_format(original_path)
                print(f"  📊 File format detected: {file_format}")
                
                success = convert_file_format(original_path, new_path, file_format)
                if success:
                    converted_files.append(new_path.name)
                    # Remove original file if conversion was successful and names are different
                    if original_path != new_path and original_path.exists():
                        original_path.unlink()
                        print(f"  🗑️  Removed original: {original_path.name}")
                else:
                    sample_success = False
        
        # Update summary
        if sample_success:
            cleanup_remaining_gsm_files(sample_dir, converted_files)
            fix_expression_file_naming(sample_dir, sample_name)
            conversion_summary['successful'] += 1
            conversion_summary['details'][sample_name] = f"Converted {len(converted_files)} files"
            print(f"  ✅ Sample {sample_name} processed successfully")
            print(f"     Files created: {', '.join(converted_files)}")
        else:
            conversion_summary['failed'] += 1
            conversion_summary['details'][sample_name] = 'Some conversions failed'
            print(f"  ❌ Sample {sample_name} had conversion errors")
    
    # Print final summary
    print("\n" + "=" * 60)
    print("📊 CONVERSION SUMMARY")
    print("=" * 60)
    print(f"✅ Successful samples: {conversion_summary['successful']}")
    print(f"❌ Failed samples: {conversion_summary['failed']}")
    print(f"📈 Success rate: {conversion_summary['successful']/(conversion_summary['successful']+conversion_summary['failed'])*100:.1f}%")
    
    print("\n📋 Detailed Results:")
    for sample, result in conversion_summary['details'].items():
        status = "✅" if "Converted" in result else "❌"
        print(f"  {status} {sample}: {result}")
    
    if conversion_summary['successful'] > 0:
        print(f"\n🚀 Ready to run pipeline:")
        print(f"python run_analysis.py")
        print(f"python spatial_transcriptomics_pipeline.py --data_dir {data_directory}")
    
    return conversion_summary['successful'] > 0

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(description='Convert spatial transcriptomics data format')
    parser.add_argument('--data_dir', type=str, default='./data', 
                       help='Directory containing sample data (default: ./data)')
    parser.add_argument('--no_backup', action='store_true', 
                       help='Skip creating backup of original data')
    parser.add_argument('--dry_run', action='store_true',
                       help='Show what would be converted without actually converting')
    
    args = parser.parse_args()
    
    if args.dry_run:
        print("🔍 DRY RUN - Showing what would be converted:")
        found_files = find_data_files(args.data_dir)
        
        for sample_name, files_info in found_files.items():
            print(f"\n📁 Sample: {sample_name}")
            sample_dir = Path(args.data_dir) / sample_name
            conversions = standardize_filenames(sample_dir, files_info)
            
            for original, new in conversions.items():
                action = "CONVERT & REMOVE" if original != new else "KEEP"
                print(f"  {original.name} → {new.name} ({action})")
    else:
        convert_all_files(args.data_dir, not args.no_backup)

if __name__ == "__main__":
    main()