#!/usr/bin/env python3
"""
Spatial Transcriptomics Cancer Analysis Pipeline
Automated pipeline for analyzing spatial gene expression in cancer tissue samples
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
import umap
import igraph as ig
import leidenalg
from PIL import Image
import os
import json
import argparse
from pathlib import Path
import warnings
from datetime import datetime
import logging
from typing import Dict, List, Tuple, Optional
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import matplotlib.patches as patches

# Disable PIL's decompression bomb check for large tissue images
Image.MAX_IMAGE_PIXELS = None
warnings.filterwarnings('ignore')

class SpatialTranscriptomicsAnalyzer:
    """Main class for spatial transcriptomics analysis"""
    
    def __init__(self, config_path: str = None):
        """Initialize the analyzer with configuration"""
        self.setup_logging()
        self.config = self.load_config(config_path)
        self.results = {}

        # Load comprehensive marker categories first
        self.comprehensive_markers = self.load_comprehensive_markers()
        
        # Load cancer marker genes from config, comprehensive file, or minimal default
        if 'cancer_markers' in self.config and len(self.config['cancer_markers']) > 20:
            # Use config markers if it has a comprehensive list
            self.cancer_markers = self.config['cancer_markers']
            self.logger.info(f"Using {len(self.cancer_markers)} cancer markers from config")
        elif self.comprehensive_markers:
            # Use all markers from comprehensive file
            all_markers = set()
            if 'cancer_markers_by_category' in self.comprehensive_markers:
                for category_markers in self.comprehensive_markers['cancer_markers_by_category'].values():
                    all_markers.update(category_markers)
            if 'cancer_type_specific' in self.comprehensive_markers:
                for type_markers in self.comprehensive_markers['cancer_type_specific'].values():
                    all_markers.update(type_markers)
            if 'pathway_markers' in self.comprehensive_markers:
                for pathway_markers in self.comprehensive_markers['pathway_markers'].values():
                    all_markers.update(pathway_markers)
            
            self.cancer_markers = sorted(list(all_markers))
            self.logger.info(f"Using {len(self.cancer_markers)} cancer markers from comprehensive database")
        else:
            # Fallback to default only if nothing else available
            self.cancer_markers = self.get_default_cancer_markers()
            self.logger.warning(f"Using default {len(self.cancer_markers)} cancer markers - comprehensive database not found")

        panglao_filepath = 'PanglaoDB_markers_27_Mar_2020.tsv' 
        
        # Define the TME cell types you want to score
        desired_tme_cells = self.config.get("desired_tme_cell_types", {
            'T-Cells': ['T cell'],
            'B-Cells': ['B cell'],
            'Macrophages': ['Macrophage'],
            'Fibroblasts': ['Fibroblast'],
            'Endothelial-Cells': ['Endothelial cell'],
            'NK-Cells': ['NK cell', 'Natural killer cell'],
            'Mast-Cells': ['Mast cell']
        })
        self.tme_marker_genes = self.load_panglaodb_markers(
            filepath=panglao_filepath,
            desired_cell_mappings=desired_tme_cells
        )
        
        
        # # Cancer marker genes (can be customized)
        # self.cancer_markers = [
        #     'TP53', 'KRAS', 'PIK3CA', 'APC', 'BRAF', 'PTEN', 'EGFR', 'MYC',
        #     'CDKN2A', 'VHL', 'RB1', 'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'ATM'
        # ]
        
    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('spatial_analysis.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def load_config(self, config_path: str) -> Dict:
        """Load configuration from JSON file"""
        default_config = {
            "data_directory": "./data/",
            "output_directory": "./results/",
            "image_format": "jpg",
            "coordinate_file_suffix": "_coordinates.csv",
            "expression_file_pattern": "*_stdata.csv",
            "image_file_pattern": "*_tissue_image.*",
            "min_genes_per_spot": 200,
            "min_spots_per_gene": 5,
            "n_clusters": 5,
            "resolution": 0.5,
            "n_top_genes": 50,
            "figure_dpi": 300,
            "analysis_types": ["spatial_plots", "clustering", "differential_expression", "pathway_analysis"]
        }
        
        if config_path and os.path.exists(config_path):
            with open(config_path, 'r') as f:
                user_config = json.load(f)
                default_config.update(user_config)
                
        return default_config
    
    def get_default_cancer_markers(self) -> List[str]:
        # \"\"\"Get default cancer markers if not specified in config\"\"\"
        return [
            'TP53', 'KRAS', 'PIK3CA', 'APC', 'BRAF', 'PTEN', 'EGFR', 'MYC',
            'CDKN2A', 'VHL', 'RB1', 'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'ATM'
        ]
    
    def load_comprehensive_markers(self) -> Dict:
        # \"\"\"Load comprehensive cancer markers from JSON file if available\"\"\"
        try:
            if os.path.exists('cancer_markers_comprehensive.json'):
                with open('cancer_markers_comprehensive.json', 'r') as f:
                    return json.load(f)
        except:
            pass
        return {}




    def load_panglaodb_markers(self, filepath: str, desired_cell_mappings: Dict[str, List[str]]) -> Dict[str, List[str]]:
        """
        Loads and parses the PanglaoDB marker file using a flexible keyword search.
        """
        self.logger.info(f"Loading TME markers from PanglaoDB file: {filepath}")
        if not Path(filepath).exists():
            self.logger.error(f"PanglaoDB marker file not found at {filepath}.")
            return {}
            
        try:
            markers_df = pd.read_csv(filepath, sep='\t')
            
            # --- THE FIX: Make the species search case-insensitive ---
            human_markers = markers_df[markers_df['species'].str.contains('Hs', na=False)].copy()
            self.logger.info(f"Found {len(human_markers)} total human marker entries in the file.")

            if human_markers.empty:
                self.logger.warning("No 'Homo sapiens' entries found. Please check the 'species' column in the TSV file.")
                return {}

            human_markers['sensitivity_human'] = pd.to_numeric(human_markers['sensitivity_human'], errors='coerce')
            human_markers.dropna(subset=['sensitivity_human'], inplace=True)

            tme_markers = {}
            for clean_name, keywords in desired_cell_mappings.items():
                pattern = '|'.join(keywords)
                matched_cells_df = human_markers[human_markers['cell type'].str.contains(pattern, case=False, na=False)]
                
                self.logger.info(f"For '{clean_name}', found {len(matched_cells_df)} potential marker entries using keywords: {keywords}")
                
                sensitivity_threshold = 0.40 
                sensitive_markers = matched_cells_df[matched_cells_df['sensitivity_human'] > sensitivity_threshold]
                self.logger.info(f"--> After filtering for sensitivity > {sensitivity_threshold}, {len(sensitive_markers)} markers remain.")
                
                if sensitive_markers.empty: continue

                top_genes = sensitive_markers.loc[sensitive_markers.groupby('official gene symbol')['sensitivity_human'].idxmax()]
                gene_list = top_genes.sort_values(by='sensitivity_human', ascending=False).head(20)['official gene symbol'].tolist()

                if gene_list:
                    tme_markers[clean_name] = gene_list

            if not tme_markers:
                 self.logger.warning("Could not load any TME markers after filtering. Check keywords in config or adjust sensitivity threshold.")
            else:
                 self.logger.info(f"Successfully loaded markers for {len(tme_markers)} cell types.")
            return tme_markers
            
        except Exception as e:
            self.logger.error(f"Error parsing PanglaoDB file: {e}", exc_info=True)
            return {}
        
    
    
    
    def extract_coordinates(self, spot_names: List[str]) -> pd.DataFrame:
        """Extract x,y coordinates from spot names like '10x12'"""
        coords = []
        for spot in spot_names:
            try:
                # Handle different coordinate formats
                if 'x' in spot:
                    x, y = spot.split('x')
                    coords.append({'spot': spot, 'x': int(x), 'y': int(y)})
                else:
                    # If no 'x' separator, try other patterns
                    coords.append({'spot': spot, 'x': 0, 'y': 0})
            except:
                coords.append({'spot': spot, 'x': 0, 'y': 0})
        
        return pd.DataFrame(coords)
    
   
        

    def load_sample_data(self, sample_path: str) -> Dict:
        """Load data for a single sample"""
        self.logger.info(f"Loading sample data from {sample_path}")
        
        sample_data = {}
        
        # Load expression data (st_data)
        expression_files = list(Path(sample_path).glob(self.config["expression_file_pattern"]))
        if not expression_files:
            raise FileNotFoundError(f"No expression data files found in {sample_path}")
        expression_file = expression_files[0]
        st_data = pd.read_csv(expression_file, index_col=0)
        
        # Load coordinate mapping file (coord_mapping)
        coord_mapping = None
        coord_files = list(Path(sample_path).glob(f"*{self.config['coordinate_file_suffix']}"))
        if coord_files:
            coord_mapping = pd.read_csv(coord_files[0])
            sample_data['coordinate_mapping'] = coord_mapping
        else:
            self.logger.error("Coordinate file not found. Cannot proceed with spatial analysis.")
            # Depending on your needs, you might want to raise an error here
            raise FileNotFoundError(f"No coordinate file found in {sample_path}")

        # =====================================================================
        # ========= START: NEW CODE TO ALIGN EXPRESSION AND COORDINATES =========
        # =====================================================================
        self.logger.info("Aligning expression data with available coordinates...")

        # 1. Create the 'spot' ID in the coordinate file from the x and y columns.
        if 'spot' not in coord_mapping.columns and 'x' in coord_mapping.columns and 'y' in coord_mapping.columns:
            coord_mapping['spot'] = coord_mapping['x'].astype(str) + 'x' + coord_mapping['y'].astype(str)

        # 2. Get a list of all spot IDs that have a valid spatial coordinate.
        valid_spatial_spots = set(coord_mapping['spot'])
        
        # 3. Filter the main expression table (st_data) to keep ONLY the spots from that list.
        original_spot_count = len(st_data)
        st_data = st_data[st_data.index.isin(valid_spatial_spots)]
        filtered_spot_count = len(st_data)
        
        self.logger.info(f"Filtered expression data to {filtered_spot_count} spots with valid coordinates (down from {original_spot_count}).")
        
        # =====================================================================
        # ========= END: NEW CODE TO ALIGN EXPRESSION AND COORDINATES ===========
        # =====================================================================

        # Clean column names if needed
        if 'Unnamed: 0' in st_data.columns:
            st_data = st_data.rename(columns={'Unnamed: 0': 'Spots'})
            st_data.index = st_data['Spots']
            st_data = st_data.drop('Spots', axis=1)
        
        sample_data['expression'] = st_data # This is now the FILTERED data
        sample_data['sample_name'] = Path(sample_path).name
        
        # Extract coordinates (this will now be based on the filtered st_data)
        coords_df = self.extract_coordinates(st_data.index.tolist())
        sample_data['coordinates'] = coords_df

        # ... (the rest of the function for loading the image continues as before) ...
        # Load tissue image if available
        image_files = list(Path(sample_path).glob(self.config["image_file_pattern"]))
        image_pattern = self.config["image_file_pattern"]
        self.logger.info(f"Looking for image files with pattern: {image_pattern}")
        image_files = list(Path(sample_path).glob(image_pattern))
        self.logger.info(f"Found {len(image_files)} image files: {[f.name for f in image_files]}")
        
        if image_files:
            self.logger.info(f"Loading image: {image_files[0].name}")
            img = Image.open(image_files[0])
            sample_data['image'] = np.array(img)
            try:
                # Load and resize large images safely
                img, resize_factor = self.load_and_resize_image(image_files[0])
                sample_data['image'] = img
                sample_data['image_path'] = str(image_files[0])
                sample_data['image_resize_factor'] = resize_factor
                if resize_factor < 1.0:
                    self.logger.info(f"Image resized by factor {resize_factor:.3f} to fit memory limits")
            except Exception as e:
                self.logger.error(f"Failed to load image {image_files[0].name}: {e}")
                sample_data['image'] = None
                sample_data['image_resize_factor'] = 1.0
        else:
            # Try a broader search to see what files exist
            all_files = list(Path(sample_path).iterdir())
            self.logger.warning(f"No tissue image found for sample {sample_path}")
            self.logger.info(f"Available files in directory: {[f.name for f in all_files if f.is_file()]}")
            sample_data['image'] = None
        
        # Load coordinate mapping file if available
        coord_files = list(Path(sample_path).glob(f"*{self.config['coordinate_file_suffix']}"))
        if coord_files:
            coord_mapping = pd.read_csv(coord_files[0])
            sample_data['coordinate_mapping'] = coord_mapping
        
        self.logger.info(f"Loaded sample: {st_data.shape[0]} spots, {st_data.shape[1]} genes")
        return sample_data
    
    def load_and_resize_image(self, image_path: Path, max_pixels: int = 100_000_000) -> Tuple[np.ndarray, float]:
        """
        Load and resize image if it's too large, maintaining aspect ratio
        Returns: (image_array, resize_factor)
        """
        
        
        try:
            # Open image and get dimensions
            img = Image.open(image_path)
            width, height = img.size
            total_pixels = width * height
            
            self.logger.info(f"Original image size: {width}x{height} ({total_pixels:,} pixels)")
            
            # Calculate resize factor if needed
            if total_pixels > max_pixels:
                resize_factor = (max_pixels / total_pixels) ** 0.5
                new_width = int(width * resize_factor)
                new_height = int(height * resize_factor)
                
                self.logger.info(f"Resizing to: {new_width}x{new_height} (factor: {resize_factor:.3f})")
                
                # Resize with high-quality resampling
                img = img.resize((new_width, new_height), Image.Resampling.LANCZOS)
            else:
                resize_factor = 1.0
                self.logger.info(f"Image size acceptable, no resizing needed")
            
            # Convert to array
            img_array = np.array(img)
            
            return img_array, resize_factor
            
        except Exception as e:
            self.logger.error(f"Error loading/resizing image {image_path}: {e}")
            raise

    def _scale_coordinates_to_image(self, coords_df: pd.DataFrame, img: np.ndarray, resize_factor: float = 1.0):
        # \"\"\"Scale coordinates to image dimensions when pixel coordinates are not available\"\"\"
        if img is not None:
            # Scale coordinates to image dimensions, accounting for image resizing
            x_scale = img.shape[1] / (coords_df['x'].max() - coords_df['x'].min() + 1)
            y_scale = img.shape[0] / (coords_df['y'].max() - coords_df['y'].min() + 1)
            
            coords_df['pixel_x'] = (coords_df['x'] - coords_df['x'].min()) * x_scale
            coords_df['pixel_y'] = (coords_df['y'] - coords_df['y'].min()) * y_scale
        else:
            coords_df['pixel_x'] = coords_df['x']
            coords_df['pixel_y'] = coords_df['y']

    def quality_control(self, sample_data: Dict) -> Dict:
        """Perform quality control on the data"""
        self.logger.info("Performing quality control...")
        
        st_data = sample_data['expression']
        
        # Calculate QC metrics
        qc_metrics = {
            'n_spots_original': st_data.shape[0],
            'n_genes_original': st_data.shape[1],
            'total_counts': st_data.sum().sum(),
            'spots_per_gene': (st_data > 0).sum(axis=0),
            'genes_per_spot': (st_data > 0).sum(axis=1),
            'mean_counts_per_spot': st_data.sum(axis=1).mean(),
            'median_counts_per_spot': st_data.sum(axis=1).median()
        }
        
        # Filter genes (present in at least min_spots_per_gene spots)
        genes_keep = qc_metrics['spots_per_gene'] >= self.config['min_spots_per_gene']
        st_data_filtered = st_data.loc[:, genes_keep]
        
        # Filter spots (express at least min_genes_per_spot genes)
        spots_keep = (st_data_filtered > 0).sum(axis=1) >= self.config['min_genes_per_spot']
        st_data_filtered = st_data_filtered.loc[spots_keep, :]
        
        qc_metrics.update({
            'n_spots_filtered': st_data_filtered.shape[0],
            'n_genes_filtered': st_data_filtered.shape[1],
            'spots_removed': qc_metrics['n_spots_original'] - st_data_filtered.shape[0],
            'genes_removed': qc_metrics['n_genes_original'] - st_data_filtered.shape[1]
        })
        
        sample_data['expression_filtered'] = st_data_filtered
        sample_data['qc_metrics'] = qc_metrics
        
        self.logger.info(f"QC complete: {st_data_filtered.shape[0]} spots, {st_data_filtered.shape[1]} genes retained")
        return sample_data
    

    def perform_leiden_clustering(self, data: np.ndarray, resolution: float = 0.5, n_neighbors: int = 15) -> Tuple[np.ndarray, int]:
        """
        Perform Leiden clustering on the data
        Returns: (cluster_labels, n_clusters)
        """
        try:
            # Build k-NN graph
            self.logger.info(f"Building k-NN graph with {n_neighbors} neighbors...")
            nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
            nn.fit(data)
            
            # Get k-NN graph
            knn_graph = nn.kneighbors_graph(data, mode='connectivity')
            
            # Convert to igraph
            sources, targets = knn_graph.nonzero()
            edges = list(zip(sources.tolist(), targets.tolist()))
            
            # Remove self-loops
            edges = [(s, t) for s, t in edges if s != t]
            
            # Create igraph object
            g = ig.Graph(n=data.shape[0], edges=edges, directed=False)
            
            # Remove duplicate edges
            g.simplify()
            
            # Perform Leiden clustering
            self.logger.info(f"Performing Leiden clustering with resolution {resolution}...")
            partition = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, 
                                                resolution_parameter=resolution, 
                                                seed=42)
            
            # Extract cluster labels
            clusters = np.array(partition.membership)
            n_clusters = len(set(clusters))
            
            self.logger.info(f"Leiden clustering completed: {n_clusters} clusters found")
            
            return clusters, n_clusters
            
        except Exception as e:
            self.logger.warning(f"Leiden clustering failed: {e}. Falling back to K-means...")
            # Fallback to K-means if Leiden fails
            kmeans = KMeans(n_clusters=self.config.get('n_clusters', 5), random_state=42)
            clusters = kmeans.fit_predict(data)
            n_clusters = self.config.get('n_clusters', 5)
            return clusters, n_clusters


    def perform_clustering(self, sample_data: Dict, output_dir: str) -> Dict:
        """Perform spatial clustering analysis"""
        self.logger.info("Performing clustering analysis...")
        
        st_data = sample_data['expression_filtered']
        coords_df = sample_data['coordinates']
        sample_name = sample_data['sample_name']
        
        # Normalize data
        scaler = StandardScaler()
        st_data_scaled = pd.DataFrame(
            scaler.fit_transform(st_data.T).T,
            index=st_data.index,
            columns=st_data.columns
        )
        
        # PCA
        pca = PCA(n_components=50)
        pca_result = pca.fit_transform(st_data_scaled)
        
        # UMAP
        umap_model = umap.UMAP(n_components=2, random_state=42)
        umap_result = umap_model.fit_transform(pca_result)
        
        # t-SNE
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, st_data.shape[0]-1))
        tsne_result = tsne.fit_transform(pca_result)

        clusters, n_clusters_found = self.perform_leiden_clustering(pca_result, self.config.get('resolution', 0.5))
        
        self.logger.info(f"Leiden clustering found {n_clusters_found} clusters (resolution: {self.config.get('resolution', 0.5)})")
        # Also perform K-means for comparison (using found cluster number)
        kmeans = KMeans(n_clusters=n_clusters_found, random_state=42)
        kmeans_clusters = kmeans.fit_predict(pca_result)

        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        
        # PCA plot with Leiden clusters
        scatter = axes[0, 0].scatter(pca_result[:, 0], pca_result[:, 1], c=clusters, cmap='viridis')
        axes[0, 0].set_title('PCA with Leiden Clusters')
        axes[0, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        axes[0, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.colorbar(scatter, ax=axes[0, 0])
        
        # UMAP plot with Leiden clusters
        scatter = axes[0, 1].scatter(umap_result[:, 0], umap_result[:, 1], c=clusters, cmap='viridis')
        axes[0, 1].set_title('UMAP with Leiden Clusters')
        axes[0, 1].set_xlabel('UMAP1')
        axes[0, 1].set_ylabel('UMAP2')
        plt.colorbar(scatter, ax=axes[0, 1])
        
        # t-SNE plot with Leiden clusters
        scatter = axes[0, 2].scatter(tsne_result[:, 0], tsne_result[:, 1], c=clusters, cmap='viridis')
        axes[0, 2].set_title('t-SNE with Leiden Clusters')
        axes[0, 2].set_xlabel('t-SNE1')
        axes[0, 2].set_ylabel('t-SNE2')
        plt.colorbar(scatter, ax=axes[0, 2])
        
        # K-means comparison plots
        scatter = axes[1, 0].scatter(pca_result[:, 0], pca_result[:, 1], c=kmeans_clusters, cmap='plasma')
        axes[1, 0].set_title('PCA with K-means Clusters')
        axes[1, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        axes[1, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        plt.colorbar(scatter, ax=axes[1, 0])

        scatter = axes[1, 1].scatter(umap_result[:, 0], umap_result[:, 1], c=kmeans_clusters, cmap='plasma')
        axes[1, 1].set_title('UMAP with K-means Clusters')
        axes[1, 1].set_xlabel('UMAP1')
        axes[1, 1].set_ylabel('UMAP2')
        plt.colorbar(scatter, ax=axes[1, 1])
        
        # Spatial clustering with Leiden clusters
        coords_filtered = coords_df[coords_df['spot'].isin(st_data.index)].copy()
        coords_filtered['cluster'] = clusters
        
        scatter = axes[1, 2].scatter(coords_filtered['x'], coords_filtered['y'], 
                                   c=coords_filtered['cluster'], cmap='viridis', s=50)
        axes[1, 2].set_title('Spatial Distribution (Leiden)')
        axes[1, 2].set_xlabel('X Coordinate')
        axes[1, 2].set_ylabel('Y Coordinate')
        plt.colorbar(scatter, ax=axes[1, 2])
        
        plt.tight_layout()
        clustering_path = os.path.join(output_dir, f'{sample_name}_clustering.png')
        plt.savefig(clustering_path, dpi=self.config['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        # Save clustering results
        clustering_results = {
            'clusters': clusters,
            'pca_result': pca_result,
            'umap_result': umap_result,
            'tsne_result': tsne_result,
            'pca_variance_ratio': pca.explained_variance_ratio_,
            'clustering_plot': clustering_path
        }
        
        return clustering_results
    

    
    
    def differential_expression_analysis(self, sample_data: Dict, clustering_results: Dict, output_dir: str) -> Dict:
        """Perform differential expression analysis between clusters"""
        self.logger.info("Performing differential expression analysis...")
        
        st_data = sample_data['expression_filtered']
        clusters = clustering_results['clusters']
        sample_name = sample_data['sample_name']
        
        # Create cluster dataframe
        cluster_df = pd.DataFrame({'spot': st_data.index, 'cluster': clusters})
        
        de_results = {}
        
        # Compare each cluster vs all others
        for cluster_id in np.unique(clusters):
            cluster_spots = cluster_df[cluster_df['cluster'] == cluster_id]['spot']
            other_spots = cluster_df[cluster_df['cluster'] != cluster_id]['spot']
            
            cluster_expr = st_data.loc[cluster_spots]
            other_expr = st_data.loc[other_spots]
            
            # Perform t-test for each gene
            de_genes = []
            for gene in st_data.columns:
                try:
                    stat, pval = stats.ttest_ind(cluster_expr[gene], other_expr[gene])
                    fold_change = cluster_expr[gene].mean() / (other_expr[gene].mean() + 1e-8)
                    
                    de_genes.append({
                        'gene': gene,
                        'cluster': cluster_id,
                        'mean_cluster': cluster_expr[gene].mean(),
                        'mean_other': other_expr[gene].mean(),
                        'fold_change': fold_change,
                        'log2_fold_change': np.log2(fold_change + 1e-8),
                        'pvalue': pval,
                        't_statistic': stat
                    })
                except:
                    continue
            
            de_df = pd.DataFrame(de_genes)
            # Adjust p-values using Benjamini-Hochberg FDR
            if not de_df.empty:
                de_df = de_df.sort_values('pvalue')
                m = len(de_df)
                de_df['pvalue_adj'] = de_df['pvalue'] * m / (np.arange(1, m+1))
                de_df['pvalue_adj'] = np.minimum(de_df['pvalue_adj'], 1.0)
            else:
                de_df['pvalue_adj'] = []
            
            de_results[f'cluster_{cluster_id}'] = de_df
        
        # Create volcano plots
        n_clusters = len(de_results)
        n_cols = min(3, n_clusters)
        n_rows = int(np.ceil(n_clusters / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
        if n_clusters == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        axes = np.array(axes).flatten()
        
        for i, (cluster_name, de_df) in enumerate(de_results.items()):
            ax = axes[i]
            
            # Volcano plot
            if not de_df.empty:
                significant = de_df['pvalue_adj'] < 0.05
                ax.scatter(de_df[~significant]['log2_fold_change'], 
                          -np.log10(de_df[~significant]['pvalue_adj'] + 1e-10), 
                          alpha=0.6, s=30, color='gray')
                ax.scatter(de_df[significant]['log2_fold_change'], 
                          -np.log10(de_df[significant]['pvalue_adj'] + 1e-10), 
                          alpha=0.8, s=30, color='red')
            
                ax.set_xlabel('Log2 Fold Change')
                ax.set_ylabel('-Log10 Adjusted P-value')
                ax.set_title(f'{cluster_name.replace("_", " ").title()} vs Others')
                ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
                ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
            else:
                ax.text(0.5, 0.5, 'No DE genes', ha='center', va='center')
        
        # Hide empty subplots
        for j in range(n_clusters, len(axes)):
            axes[j].axis('off')
        
        plt.tight_layout()
        de_plot_path = os.path.join(output_dir, f'{sample_name}_differential_expression.png')
        plt.savefig(de_plot_path, dpi=self.config['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        # Save DE results
        de_summary_path = os.path.join(output_dir, f'{sample_name}_de_results.csv')
        all_de_results = pd.concat([df.assign(comparison=cluster) for cluster, df in de_results.items()])
        all_de_results.to_csv(de_summary_path, index=False)
        
        return {
            'de_results': de_results,
            'de_plot': de_plot_path,
            'de_summary_file': de_summary_path
        }
    
    def analyze_cancer_markers(self, sample_data: Dict, output_dir: str) -> Dict:

    # Analyzes cancer marker genes by identifying significant markers based on expression
    # and prevalence thresholds, and generates relevant plots and statistics.
   
        self.logger.info("Analyzing cancer markers...")
        
        st_data = sample_data['expression_filtered']
        sample_name = sample_data['sample_name']
        
        # 1. Identify all potential markers present in the dataset
        potential_markers = [gene for gene in self.cancer_markers if gene in st_data.columns]
        self.logger.info(f"{len(potential_markers)} potential markers from the database were found in the dataset.")

        if not potential_markers:
            self.logger.warning("No cancer marker genes from the database were found in this sample's data.")
            return {'significant_markers': [], 'marker_analysis': None}
            
        # 2. Calculate initial expression statistics for all potential markers
        marker_expr_potential = st_data[potential_markers]
        marker_stats_potential = pd.DataFrame({
            'gene': potential_markers,
            'mean_expression': marker_expr_potential.mean(),
            'percentage_expressing': (marker_expr_potential > 0).sum() / len(marker_expr_potential) * 100
        })

        # 3. Apply thresholds to filter for biologically significant markers
        mean_thresh = self.config.get('marker_mean_expression_threshold', 0.5)
        prevalence_thresh = self.config.get('marker_prevalence_threshold_percent', 3)

        significant_markers_df = marker_stats_potential[
            (marker_stats_potential['mean_expression'] >= mean_thresh) &
            (marker_stats_potential['percentage_expressing'] >= prevalence_thresh)
        ]
        
        significant_markers = significant_markers_df['gene'].tolist()

        self.logger.info(f"Found {len(significant_markers)} significant markers based on thresholds "
                        f"(Mean Expr > {mean_thresh}, Prevalence > {prevalence_thresh}%)")

        if not significant_markers:
            self.logger.warning("No potential markers passed the significance thresholds.")
            return {
                'potential_markers': potential_markers,
                'significant_markers': [], 
                'marker_stats': marker_stats_potential
            }

        # 4. Proceed with analysis and plotting using ONLY the significant markers
        self.logger.info(f"Generating analysis for {len(significant_markers)} significant markers.")
        significant_marker_expr = st_data[significant_markers]
        
        
        # Calculate full summary statistics for the significant markers for the report
        final_marker_stats = pd.DataFrame({
            'gene': significant_markers,
            'mean_expression': significant_marker_expr.mean(),
            'std_expression': significant_marker_expr.std(),
            'max_expression': significant_marker_expr.max(),
            'spots_expressing': (significant_marker_expr > 0).sum(),
            'percentage_expressing': (significant_marker_expr > 0).sum() / len(significant_marker_expr) * 100
        }).sort_values(by='mean_expression', ascending=False)
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Heatmap of significant marker correlation
        if len(significant_markers) > 1:
            correlation_matrix = significant_marker_expr.corr()
            sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, ax=axes[0, 0], annot_kws={"size": 8})
            axes[0, 0].set_title('Significant Cancer Marker Gene Correlation')
        else:
            axes[0, 0].text(0.5, 0.5, 'Insufficient markers for correlation plot', 
                        ha='center', va='center', transform=axes[0, 0].transAxes)
            axes[0, 0].set_title('Significant Cancer Marker Gene Correlation')
        
        # Expression distribution of significant markers
        significant_marker_expr.boxplot(ax=axes[0, 1], rot=45)
        axes[0, 1].set_title('Significant Marker Expression Distribution')
        axes[0, 1].set_ylabel('Expression Level')
        
        # Percentage of spots expressing each significant marker
        axes[1, 0].bar(final_marker_stats['gene'], final_marker_stats['percentage_expressing'])
        axes[1, 0].set_title('Percentage of Spots Expressing Each Marker')
        axes[1, 0].set_ylabel('Percentage (%)')
        axes[1, 0].tick_params(axis='x', rotation=45, labelsize=8)
        
        # Mean expression levels of significant markers
        axes[1, 1].bar(final_marker_stats['gene'], final_marker_stats['mean_expression'])
        axes[1, 1].set_title('Mean Expression Levels of Markers')
        axes[1, 1].set_ylabel('Mean Expression')
        axes[1, 1].tick_params(axis='x', rotation=45, labelsize=8)
        
        plt.tight_layout(pad=3.0)
        fig.suptitle(f'Significant Cancer Marker Analysis for {sample_name}', fontsize=16)
        plt.subplots_adjust(top=0.92)

        markers_plot_path = os.path.join(output_dir, f'{sample_name}_cancer_markers.png')
        plt.savefig(markers_plot_path, dpi=self.config['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        # Save statistics for significant markers
        marker_stats_path = os.path.join(output_dir, f'{sample_name}_cancer_marker_stats.csv')
        final_marker_stats.to_csv(marker_stats_path, index=False)
        
        return {
            'potential_markers': potential_markers,
            'significant_markers': significant_markers,
            'marker_stats': marker_stats_potential, # Full stats for all potential markers for inspection
            'markers_plot': markers_plot_path,
            'marker_stats_file': marker_stats_path
        }
    


    def create_spatial_plots(self, sample_data: Dict, output_dir: str) -> Dict:
    # """
    # Create spatial gene expression plots. This function is robustly designed to handle
    # both barcode-based (e.g., Visium) and array-based ('XxY') coordinate formats.
    # It corrects for coordinate system differences and uses semi-transparent spots
    # for a clear visual overlay of expression on the tissue image.
    # """
        self.logger.info("Creating spatial plots with universal coordinate handler and enhanced visualization...")

        st_data = sample_data['expression_filtered']
        coord_mapping = sample_data.get('coordinate_mapping')
        img = sample_data.get('image')
        sample_name = sample_data['sample_name']

        if coord_mapping is None:
            self.logger.error("A *_coordinates.csv file is required but was not found. Skipping spatial plots.")
            return {'plots_created': [], 'genes_plotted': []}

        if 'spot' not in coord_mapping.columns and 'x' in coord_mapping.columns and 'y' in coord_mapping.columns:
            self.logger.info("Detected array-based coordinates ('x', 'y'). Creating 'spot' column for merging.")
            coord_mapping['spot'] = coord_mapping['x'].astype(str) + 'x' + coord_mapping['y'].astype(str)

        st_data_with_spots = st_data.reset_index().rename(columns={'index': 'spot'})
        merged_data = pd.merge(st_data_with_spots, coord_mapping, on='spot', how='inner')

        if merged_data.empty:
            self.logger.error("Merge failed: Spot names/barcodes in expression data do not match any in the coordinates file.")
            return {'plots_created': [], 'genes_plotted': []}
        
        self.logger.info(f"Successfully merged expression and coordinate data for {len(merged_data)} spots.")

        resize_factor = sample_data.get('image_resize_factor', 1.0)
        merged_data['plot_x'] = merged_data['pixel_x'] * resize_factor
        merged_data['plot_y'] = merged_data['pixel_y'] * resize_factor

        plots_created = []
        
        genes_plotted = []
        significant_markers = self.results.get(sample_name, {}).get('cancer_analysis', {}).get('significant_markers', [])
        if significant_markers:
            top_markers = st_data[significant_markers].mean().nlargest(16).index.tolist()
            genes_plotted = top_markers
        else:
            genes_plotted = st_data.sum(axis=0).nlargest(16).index.tolist()

        self.logger.info(f"Plotting top genes: {', '.join(genes_plotted)}")
        
        n_cols = 4
        n_rows = int(np.ceil(len(genes_plotted) / n_cols))
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 5))
        axs = np.array(axs).flatten()

        for i, gene in enumerate(genes_plotted):
            if gene in merged_data.columns:
                if img is not None:
                    axs[i].imshow(img) # Display the image with its natural orientation
                
                    # --- ENHANCEMENT: Spots are now semi-transparent for a true overlay effect ---
                    sc = axs[i].scatter(
                        merged_data['plot_x'],
                        merged_data['plot_y'], # Use direct coordinates, as imshow now handles orientation
                        c=merged_data[gene],
                        cmap='viridis',
                        s=25,               # Spot size
                        alpha=0.65,         # Make spots semi-transparent (65% opaque)
                        edgecolor='black',
                        linewidth=0.2
                    )
                else: # Fallback for no image
                    sc = axs[i].scatter(
                        merged_data['plot_x'], merged_data['plot_y'],
                        c=merged_data[gene], cmap='viridis', s=25
                    )

                axs[i].set_title(f'{gene}', fontsize=12)
                axs[i].axis('off')
                plt.colorbar(sc, ax=axs[i], fraction=0.046, pad=0.04)
                
        for j in range(len(genes_plotted), len(axs)):
            axs[j].axis('off')

        plt.tight_layout(pad=2.0)
        fig.suptitle(f'Spatial Gene Expression for {sample_name}', fontsize=16)
        plt.subplots_adjust(top=0.93)

        plot_path = os.path.join(output_dir, f'{sample_name}_spatial_expression.png')
        plt.savefig(plot_path, dpi=self.config['figure_dpi'], bbox_inches='tight')
        plt.close()
        plots_created.append(plot_path)
        
        self.results.setdefault(sample_name, {})['merged_spatial_data'] = merged_data
        
        return {'plots_created': plots_created, 'genes_plotted': genes_plotted}
    

    

    def analyze_tme(self, sample_data: Dict, output_dir: str) -> Dict:
        """
        Analyze the tumor microenvironment by scoring spots for various cell types.
        """
        self.logger.info("Analyzing tumor microenvironment (TME)...")
        sample_name = sample_data['sample_name']
        st_data = sample_data['expression_filtered']
        img = sample_data.get('image')
        
        # We need the merged data with pixel coordinates, which is in self.results
        merged_data = self.results.get(sample_name, {}).get('merged_spatial_data')
        if merged_data is None:
            self.logger.error("Merged spatial data not found. Skipping TME analysis.")
            return {}

        # # 2. Convert pandas DataFrame to AnnData for scanpy processing
        # adata = sc.AnnData(st_data)
        tme_marker_genes = self.tme_marker_genes
        if not tme_marker_genes:
            self.logger.warning("TME marker gene dictionary is empty. Skipping TME analysis.")
            return {}

        adata = sc.AnnData(st_data)

        for cell_type, markers in tme_marker_genes.items():
            markers_in_data = [m for m in markers if m in adata.var_names]
            if markers_in_data:
                sc.tl.score_genes(adata, gene_list=markers_in_data, score_name=f"{cell_type}_Score", ctrl_size=min(50, len(adata.var_names)-1))
                self.logger.info(f"Calculated enrichment score for {cell_type}")

        cell_types_with_scores = [col for col in adata.obs.columns if col.endswith('_Score')]
        if not cell_types_with_scores:
            self.logger.warning("No TME scores were calculated. Skipping visualization.")
            return {}
            
        n_cols = min(4, len(cell_types_with_scores))
        n_rows = int(np.ceil(len(cell_types_with_scores) / n_cols))
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 5))
        axs = np.array(axs).flatten()

        plot_data = pd.merge(merged_data, adata.obs[cell_types_with_scores], left_on='spot', right_index=True)

        for i, score_name in enumerate(cell_types_with_scores):
            if sample_data['image'] is not None:
                axs[i].imshow(sample_data['image'])
            scatter_plot = axs[i].scatter(
                plot_data['plot_x'],
                plot_data['plot_y'],
                c=plot_data[score_name],
                cmap='viridis',
                s=25,
                alpha=0.7
            )
            axs[i].set_title(score_name.replace('_', ' '))
            axs[i].axis('off')
            plt.colorbar(scatter_plot, ax=axs[i])
        
        for j in range(len(cell_types_with_scores), len(axs)): axs[j].axis('off')

        plt.tight_layout()
        plot_path = Path(output_dir) / f"{sample_name}_tme_spatial_scores.png"
        plt.savefig(plot_path, dpi=self.config['figure_dpi'])
        plt.close()

        scores_path = Path(output_dir) / f"{sample_name}_tme_scores.csv"
        adata.obs.to_csv(scores_path)
        
        self.logger.info(f"TME analysis complete. Plot saved to {plot_path}")
        return {'tme_plot': str(plot_path), 'tme_scores_file': str(scores_path)}



    def generate_report(self, sample_results: Dict, output_dir: str):
        """Generate comprehensive analysis report for a single sample."""
        self.logger.info(f"Generating analysis report for {sample_results['sample_name']}...")
        
        sample_name = sample_results['sample_name']
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report_content = f"""# Spatial Transcriptomics Analysis Report

    **Sample:** {sample_name}  
    **Generated:** {timestamp}

    ## Summary Statistics

    """
        
        # Add QC metrics
        if 'qc_metrics' in sample_results:
            qc = sample_results['qc_metrics']
            report_content += f"""
    ### Quality Control Metrics

    - **Original data:** {qc.get('n_spots_original', 'N/A'):,} spots, {qc.get('n_genes_original', 'N/A'):,} genes
    - **Filtered data:** {qc.get('n_spots_filtered', 'N/A'):,} spots, {qc.get('n_genes_filtered', 'N/A'):,} genes
    - **Spots removed:** {qc.get('spots_removed', 'N/A'):,} ({qc.get('spots_removed', 0)/qc.get('n_spots_original', 1)*100:.1f}%)
    - **Genes removed:** {qc.get('genes_removed', 'N/A'):,} ({qc.get('genes_removed', 0)/qc.get('n_genes_original', 1)*100:.1f}%)
    - **Mean counts per spot:** {qc.get('mean_counts_per_spot', 'N/A'):.1f}
    - **Median counts per spot:** {qc.get('median_counts_per_spot', 'N/A'):.1f}

    """
        
        # Add clustering results
        if 'clustering_results' in sample_results:
            clustering = sample_results['clustering_results']
            pca_var = clustering.get('pca_variance_ratio', [0, 0])[:2].sum()
            n_clusters_found = len(np.unique(clustering.get('clusters', [])))
            report_content += f"""
    ### Clustering Analysis

    - **Clusters found (Leiden):** {n_clusters_found}
    - **PCA variance explained (PC1+PC2):** {pca_var:.1%}
    - **Dimensionality reduction:** PCA → UMAP, t-SNE

    """
        
        # Add cancer marker analysis using the new keys
        if 'cancer_analysis' in sample_results:
            cancer = sample_results.get('cancer_analysis', {})
            potential_markers_count = len(cancer.get('potential_markers', []))
            significant_markers = cancer.get('significant_markers', [])
            
            report_content += f"""
    ### Cancer Marker Analysis

    - **Potential markers found in data:** {potential_markers_count}
    - **Significant markers (passing thresholds):** {len(significant_markers)}
    """
            # Only list the markers if there are any to show
            if significant_markers:
                # Show top 20 significant markers for brevity in the report
                markers_to_show = significant_markers[:20]
                report_content += f"- **Significant markers found:** {', '.join(markers_to_show)}"
                if len(significant_markers) > 10:
                    report_content += f", ... (and {len(significant_markers) - 20} more)"
            report_content += "\n\n"

        # Add differential expression summary
        if 'de_analysis' in sample_results:
            de = sample_results['de_analysis']
            report_content += f"""
    ### Differential Expression Analysis

    - **Comparisons performed:** {len(de.get('de_results', {}))} cluster comparisons (each vs. all others)
    - **Significant genes per cluster:** Variable (see detailed results in the CSV file)

    """
        
        # Add file paths
        report_content += "## Generated Files\n\n"
        
        if 'spatial_plots' in sample_results and sample_results['spatial_plots'].get('plots_created'):
            report_content += f"- **Spatial expression plots:** `{os.path.basename(sample_results['spatial_plots']['plots_created'][0])}`\n"
        
        if 'clustering_results' in sample_results and sample_results['clustering_results'].get('clustering_plot'):
            report_content += f"- **Clustering analysis:** `{os.path.basename(sample_results['clustering_results']['clustering_plot'])}`\n"
        
        if 'de_analysis' in sample_results and sample_results['de_analysis'].get('de_plot'):
            report_content += f"- **Differential expression plots:** `{os.path.basename(sample_results['de_analysis']['de_plot'])}`\n"
            report_content += f"- **DE results table:** `{os.path.basename(sample_results['de_analysis']['de_summary_file'])}`\n"
        
        if 'cancer_analysis' in sample_results and sample_results['cancer_analysis'].get('markers_plot'):
            report_content += f"- **Cancer marker analysis plots:** `{os.path.basename(sample_results['cancer_analysis']['markers_plot'])}`\n"
            report_content += f"- **Marker statistics table:** `{os.path.basename(sample_results['cancer_analysis']['marker_stats_file'])}`\n"
        
        # Save report
        report_path = os.path.join(output_dir, f'{sample_name}_analysis_report.md')
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        self.logger.info(f"Report generated successfully: {report_path}")
        return report_path
    

    


    def process_sample(self, sample_path: str) -> Dict:
        """Process a single sample through the entire pipeline"""
        sample_name = Path(sample_path).name
        self.logger.info(f"Processing sample: {sample_name}")
        
        output_dir = os.path.join(self.config['output_directory'], sample_name)
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Load data and perform QC (no changes here)
            sample_data = self.load_sample_data(sample_path)
            sample_data = self.quality_control(sample_data)
            
            results = {
                'sample_name': sample_name,
                'sample_path': sample_path,
                'output_directory': output_dir,
                'qc_metrics': sample_data['qc_metrics']
            }
            
            # --- START: CORRECTED ORDER OF OPERATIONS ---

            # 1. Run clustering and DE analysis first
            if 'clustering' in self.config['analysis_types']:
                clustering_results = self.perform_clustering(sample_data, output_dir)
                results['clustering_results'] = clustering_results
            
            if 'differential_expression' in self.config['analysis_types'] and 'clustering_results' in results:
                de_results = self.differential_expression_analysis(
                    sample_data, results['clustering_results'], output_dir
                )
                results['de_analysis'] = de_results

            # 2. Run cancer marker analysis to find significant genes
            cancer_results = self.analyze_cancer_markers(sample_data, output_dir)
            results['cancer_analysis'] = cancer_results
            
            # 3. NOW create spatial plots, which can use the cancer_results
            if 'spatial_plots' in self.config['analysis_types']:
                # IMPORTANT: We store the cancer results in self.results so the plotting function can find it
                self.results[sample_name] = results
                spatial_results = self.create_spatial_plots(sample_data, output_dir)
                results['spatial_plots'] = spatial_results

             # 4. Perform the TME analysis using the final dataset
            tme_results = self.analyze_tme(sample_data, output_dir)
            results['tme_analysis'] = tme_results
            
            # --- END: CORRECTED ORDER OF OPERATIONS ---

            # 4. Finally, generate the report
            report_path = self.generate_report(results, output_dir)
            results['report_path'] = report_path
            
            self.logger.info(f"Sample {sample_name} processed successfully")
            return results
            
        except Exception as e:
            self.logger.error(f"Error processing sample {sample_name}: {str(e)}", exc_info=True)
            return {'sample_name': sample_name, 'error': str(e)}
        

    def run_batch_analysis(self, data_directory: str = None) -> Dict:
        """Run analysis on all samples in the data directory"""
        if data_directory is None:
            data_directory = self.config['data_directory']
        
        self.logger.info(f"Starting batch analysis in directory: {data_directory}")
        
        # Find all sample directories
        sample_dirs = [d for d in Path(data_directory).iterdir() if d.is_dir()]
        
        if not sample_dirs:
            self.logger.error(f"No sample directories found in {data_directory}")
            return {}
        
        self.logger.info(f"Found {len(sample_dirs)} sample directories")
        
        batch_results = {}
        
        for sample_dir in sample_dirs:
            sample_results = self.process_sample(str(sample_dir))
            batch_results[sample_results['sample_name']] = sample_results
        
        # Generate batch summary report
        self.generate_batch_report(batch_results)
        
        return batch_results
    
    def generate_batch_report(self, batch_results: Dict):
        """Generate summary report for batch analysis"""
        self.logger.info("Generating batch summary report...")
        
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report_content = f"""# Batch Spatial Transcriptomics Analysis Report

**Generated:** {timestamp}  
**Total Samples:** {len(batch_results)}

## Batch Summary

"""
        
        successful_samples = [name for name, results in batch_results.items() if 'error' not in results]
        failed_samples = [name for name, results in batch_results.items() if 'error' in results]
        
        report_content += f"""
- **Successful analyses:** {len(successful_samples)}
- **Failed analyses:** {len(failed_samples)}
- **Success rate:** {len(successful_samples)/len(batch_results)*100:.1f}%

"""
        
        if failed_samples:
            report_content += f"""
### Failed Samples
{chr(10).join([f"- {sample}" for sample in failed_samples])}

"""
        
        # Summary statistics across samples
        if successful_samples:
            report_content += """
### Cross-Sample Summary

| Sample | Spots | Genes | Cancer Markers | Clusters | Status |
|--------|-------|-------|----------------|----------|--------|
"""
            
            for sample_name in successful_samples:
                results = batch_results[sample_name]
                spots = results.get('qc_metrics', {}).get('n_spots_filtered', 'N/A')
                genes = results.get('qc_metrics', {}).get('n_genes_filtered', 'N/A')
                markers = len(results.get('cancer_analysis', {}).get('available_markers', []))
                clusters = self.config.get('n_clusters', 'N/A')
                
                report_content += f"| {sample_name} | {spots} | {genes} | {markers} | {clusters} | ✓ |\n"
        
        # Save batch report
        report_path = os.path.join(self.config['output_directory'], 'batch_analysis_report.md')
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        self.logger.info(f"Batch report generated: {report_path}")


def main():
    """Main function to run the pipeline"""
    parser = argparse.ArgumentParser(description='Spatial Transcriptomics Cancer Analysis Pipeline')
    parser.add_argument('--data_dir', type=str, help='Directory containing sample data')
    parser.add_argument('--output_dir', type=str, help='Output directory for results')
    parser.add_argument('--config', type=str, help='Configuration file path')
    parser.add_argument('--sample', type=str, help='Process single sample directory')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = SpatialTranscriptomicsAnalyzer(config_path=args.config)
    
    # Update config with command line arguments
    if args.data_dir:
        analyzer.config['data_directory'] = args.data_dir
    if args.output_dir:
        analyzer.config['output_directory'] = args.output_dir
    
    # Create output directory
    os.makedirs(analyzer.config['output_directory'], exist_ok=True)
    
    if args.sample:
        # Process single sample
        results = analyzer.process_sample(args.sample)
        print(f"Single sample analysis complete. Results: {results}")
    else:
        # Run batch analysis
        results = analyzer.run_batch_analysis(args.data_dir)
        print(f"Batch analysis complete. Processed {len(results)} samples.")


if __name__ == "__main__":
    main()
    