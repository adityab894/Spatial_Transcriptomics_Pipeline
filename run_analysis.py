#!/usr/bin/env python3
"""
Simple runner script for the Spatial Transcriptomics Pipeline
"""

import os
import sys
from pathlib import Path
from Spatial_Transcriptomics_Pipeline import SpatialTranscriptomicsAnalyzer

def setup_directories():
    """Create necessary directories"""
    dirs = ['data', 'results']
    for dir_name in dirs:
        os.makedirs(dir_name, exist_ok=True)
        print(f"✓ Created directory: {dir_name}")

def main():
    print("🧬 Spatial Transcriptomics Cancer Analysis Pipeline")
    print("=" * 50)
    
    # Setup directories
    setup_directories()
    
    # Check if data directory has samples
    data_dir = Path('./data')
    sample_dirs = [d for d in data_dir.iterdir() if d.is_dir()]
    
    if not sample_dirs:
        print("⚠️  No sample directories found in ./data/")
        print("Please organize your data as follows:")
        print("data/")
        print("├── sample1/")
        print("│   ├── *_stdata.csv")
        print("│   ├── *.jpg")
        print("│   └── *_coordinates.csv (optional)")
        print("├── sample2/")
        print("│   ├── *_stdata.csv")
        print("│   ├── *.jpg")
        print("│   └── *_coordinates.csv (optional)")
        print("└── ...")
        return
    
    print(f"\n📁 Found {len(sample_dirs)} sample directories:")
    for sample_dir in sample_dirs:
        print(f"   - {sample_dir.name}")
    
    # Initialize and run analyzer
    print(f"\n🚀 Starting analysis...")
    analyzer = SpatialTranscriptomicsAnalyzer(config_path='config.json')
    
    # Run batch analysis
    results = analyzer.run_batch_analysis()
    
    # Print summary
    successful = len([r for r in results.values() if 'error' not in r])
    failed = len(results) - successful
    
    print(f"\n✅ Analysis Complete!")
    print(f"   - Successful: {successful}")
    print(f"   - Failed: {failed}")
    print(f"   - Results saved in: ./results/")
    
    if successful > 0:
        print(f"\n📊 Generated files for each sample:")
        print(f"   - Spatial expression plots (*_spatial_expression.png)") 
        print(f"   - Clustering analysis (*_clustering.png)")
        print(f"   - Differential expression (*_differential_expression.png)") 
        print(f"   - Cancer marker analysis (*_cancer_markers.png)")
        print(f"   - Analysis report (*_analysis_report.md)")
        print(f"   - Data tables (*.csv)")

if __name__== "__main__":
    main()