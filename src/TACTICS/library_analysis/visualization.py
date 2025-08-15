import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import polars as pl
import altair as alt
from typing import List, Optional, Tuple, Dict, Union, Literal, Any
from .library_analysis import LibraryAnalysis
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display

class LibraryVisualization:
    """
    A class for creating visualizations from multiple LibraryAnalysis instances.
    """
    
    def __init__(self, analysis: LibraryAnalysis):
        """
        Initialize the LibraryVisualization class.
        
        Parameters:
        -----------
        analysis : LibraryAnalysis
            LibraryAnalysis instance containing the data to visualize
        """
        self.analysis = analysis

    def _sort_dataframe(self, df: Any, score_column: str, ascending: bool = False) -> Any:
        """
        Sort a DataFrame by the specified column, handling both pandas and polars DataFrames.
        
        Parameters:
        -----------
        df : Any
            DataFrame to sort (pandas or polars)
        score_column : str
            Column to sort by
        ascending : bool
            Whether to sort in ascending order
            
        Returns:
        --------
        Any
            Sorted DataFrame
        """
        # Check if DataFrame is pandas
        if isinstance(df, pd.DataFrame):
            return df.sort_values(by=score_column, ascending=ascending)
        # Check if DataFrame is polars
        elif hasattr(df, 'sort'):
            return df.sort(score_column, descending=not ascending)
        else:
            raise TypeError(f"Unsupported DataFrame type: {type(df)}")
    
    def _get_top_products(self, df: Any, product_code_column: str, score_column: str, 
                         top_n: int, ascending: bool = False) -> set:
        """
        Get the top N products from a DataFrame.
        
        Parameters:
        -----------
        df : Any
            DataFrame to extract products from
        product_code_column : str
            Column containing product codes
        score_column : str
            Column to sort by
        top_n : int
            Number of top products to extract
        ascending : bool
            Whether to sort in ascending order
            
        Returns:
        --------
        set
            Set of top product codes
        """
        # Sort the DataFrame
        sorted_df = self._sort_dataframe(df, score_column, ascending)
        
        # Extract top N products
        if isinstance(sorted_df, pd.DataFrame):
            top_products = sorted_df[product_code_column].head(top_n).tolist()
        else:  # polars DataFrame
            top_products = sorted_df[product_code_column].head(top_n).to_list()
            
        return set(top_products)
    
    def plot_top_products_comparison(
        self,
        analysis_instances: Union[List[LibraryAnalysis], List[List[LibraryAnalysis]]],
        reference_instance: LibraryAnalysis,
        top_n: int = 100,
        title: str = "Top Products Comparison",
        figsize: Tuple[int, int] = (15, 10),
        analysis_labels: Optional[List[str]] = None,
        save_path: Optional[str] = None
    ) -> None:
        """
        Generate subplots comparing the overlap of top products between analysis instances
        and a reference instance using the compare_analysis_overlap method.
        Each subplot represents a different building block position.
        
        Parameters:
        -----------
        analysis_instances : Union[List[LibraryAnalysis], List[List[LibraryAnalysis]]]
            Either a single list of LibraryAnalysis instances to compare (e.g., Thompson sampling replicates)
            or a list of lists, where each inner list contains instances of the same type
            (e.g., [thompson_standard_instances, thompson_boltzmann_instances])
        reference_instance : LibraryAnalysis
            A single reference LibraryAnalysis instance (e.g., standard docking result)
        top_n : int
            Number of top products to consider for comparison
        title : str
            Title for the plot
        figsize : Tuple[int, int]
            Figure size (width, height)
        analysis_labels : Optional[List[str]]
            Labels for each analysis group (e.g., ["Thompson Standard", "Thompson Boltzmann"])
        save_path : Optional[str]
            Path to save the plot (if None, plot is displayed)
        """
        if reference_instance is None:
            raise ValueError("reference_instance must be provided")
            
        # Handle different input formats for analysis_instances
        if isinstance(analysis_instances[0], LibraryAnalysis):
            # Single list of instances
            analysis_groups = [analysis_instances]
            if analysis_labels is None:
                analysis_labels = ["Analysis"]
        else:
            # List of lists of instances
            analysis_groups = analysis_instances
            if analysis_labels is None:
                analysis_labels = [f"Analysis Group {i+1}" for i in range(len(analysis_groups))]
        
        if len(analysis_labels) != len(analysis_groups):
            raise ValueError("Number of analysis_labels must match number of analysis groups")
        
        # Get the first analysis instance to determine the number of positions
        first_analysis = analysis_groups[0][0]
        overlap_results = first_analysis.compare_analysis_overlap(reference_instance, top_n=top_n)
        num_positions = len(overlap_results)
        
        # Initialize data structures to store results for each position
        position_data = {i: [] for i in range(num_positions)}
        
        # Calculate overlaps for each analysis group
        for group_idx, group in enumerate(analysis_groups):
            group_label = analysis_labels[group_idx]
            
            # Calculate overlaps for each instance in the group
            for i, analysis in enumerate(group):
                # Use the compare_analysis_overlap method from LibraryAnalysis
                overlap_results = analysis.compare_analysis_overlap(reference_instance, top_n=top_n)
                
                # Store results for each position
                for position, num_overlap, _ in overlap_results:
                    overlap_percentage = (num_overlap / top_n) * 100
                    position_data[position].append({
                        'Instance': f"Instance {i+1}",
                        'Group': group_label,
                        'Overlap Percentage': overlap_percentage
                    })
        
        # Create the figure with subplots
        fig, axes = plt.subplots(num_positions, 1, figsize=figsize, sharex=True)
        if num_positions == 1:
            axes = [axes]  # Make axes a list for consistent indexing
        
        # Set the main title
        fig.suptitle(title, fontsize=16)
        
        # Create a barplot for each position
        for position, data in position_data.items():
            ax = axes[position]
            
            # Convert position data to DataFrame
            plot_data = pd.DataFrame(data)
            
            # Create grouped barplot with Set1 palette
            sns.barplot(
                data=plot_data,
                x='Instance',
                y='Overlap Percentage',
                hue='Group',
                palette="Set1",
                ax=ax
            )
            
            # Customize the subplot
            ax.set_title(f"Position {position + 1}")
            ax.set_ylabel('Overlap Percentage (%)')
            ax.set_ylim(0, 100)  # Set y-axis from 0 to 100%
            
            # Add value labels on top of bars
            for i, v in enumerate(plot_data['Overlap Percentage']):
                ax.text(i, v + 1, f'{v:.1f}%', ha='center', fontsize=8)
            
            # Rotate x-axis labels for better readability
            ax.tick_params(axis='x', rotation=45)
        
        # Set x-axis label for the bottom subplot
        axes[-1].set_xlabel('Analysis Instance')
        
        # Adjust layout
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for the suptitle
        
        # Save or show the plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
        
        # Print summary statistics by group and position
        print("\nSummary Statistics by Position:")
        for position in range(num_positions):
            print(f"\nPosition {position + 1}:")
            position_df = pd.DataFrame(position_data[position])
            
            for group in analysis_labels:
                group_data = position_df[position_df['Group'] == group]['Overlap Percentage']
                print(f"\n  {group}:")
                print(f"    Mean overlap: {group_data.mean():.2f}%")
                print(f"    Std overlap: {group_data.std():.2f}%")
                print(f"    Min overlap: {group_data.min():.2f}%")
                print(f"    Max overlap: {group_data.max():.2f}%")
            
            # Print overall statistics for this position
            print(f"\n  Overall:")
            print(f"    Mean overlap: {position_df['Overlap Percentage'].mean():.2f}%")
            print(f"    Std overlap: {position_df['Overlap Percentage'].std():.2f}%")
            print(f"    Min overlap: {position_df['Overlap Percentage'].min():.2f}%")
            print(f"    Max overlap: {position_df['Overlap Percentage'].max():.2f}%")
        
        # Print overall statistics across all positions
        print("\nOverall Statistics (All Positions):")
        all_data = []
        for position_data_list in position_data.values():
            all_data.extend(position_data_list)
        
        all_df = pd.DataFrame(all_data)
        
        for group in analysis_labels:
            group_data = all_df[all_df['Group'] == group]['Overlap Percentage']
            print(f"\n{group}:")
            print(f"  Mean overlap: {group_data.mean():.2f}%")
            print(f"  Std overlap: {group_data.std():.2f}%")
            print(f"  Min overlap: {group_data.min():.2f}%")
            print(f"  Max overlap: {group_data.max():.2f}%")
        
        print("\nOverall:")
        print(f"  Mean overlap: {all_df['Overlap Percentage'].mean():.2f}%")
        print(f"  Std overlap: {all_df['Overlap Percentage'].std():.2f}%")
        print(f"  Min overlap: {all_df['Overlap Percentage'].min():.2f}%")
        print(f"  Max overlap: {all_df['Overlap Percentage'].max():.2f}%")

    def visualize_top_building_blocks(self,
                                    show_overlap: bool = False,
                                    mols_per_row: int = 5,
                                    sub_img_size: Tuple[int, int] = (300, 300),
                                    comparison_analysis: Optional[LibraryAnalysis] = None,
                                    top_n: int = 20) -> None:
        """
        Visualize the top building blocks for each position using RDKit.
        
        Parameters:
        -----------
        show_overlap : bool
            If True, shows overlapping building blocks
            If False, shows top building blocks
        mols_per_row : int
            Number of molecules to display per row
        sub_img_size : Tuple[int, int]
            Size of each molecule image
        comparison_analysis : Optional[LibraryAnalysis]
            Another LibraryAnalysis instance to compare with for overlap visualization
        top_n : int
            Number of top building blocks to consider for overlap visualization
        """
        if show_overlap:
            if comparison_analysis is not None:
                # Get overlapping building blocks from comparison
                overlap_results = self.analysis.compare_analysis_overlap(comparison_analysis, top_n=top_n)
                
                for position, num_overlap, overlap_blocks in overlap_results:
                    if not overlap_blocks:
                        print(f"No overlapping building blocks at position {position}")
                        continue
                        
                    print(f"Overlapping Building Blocks at Position {position} ({num_overlap} overlaps)")
                    # Create RDKit molecules for overlapping blocks
                    mols = []
                    legends = []
                    for block in overlap_blocks:
                        if block in self.analysis.smiles_dict:
                            mol = Chem.MolFromSmiles(self.analysis.smiles_dict[block])
                            if mol is not None:
                                mols.append(mol)
                                legends.append(f"{block}")
                    
                    if not mols:
                        print(f"No valid SMILES found for overlapping blocks at position {position}")
                        continue
                    
                    # Create and display the grid image
                    img = Draw.MolsToGridImage(
                        mols, 
                        legends=legends,
                        molsPerRow=mols_per_row, 
                        subImgSize=sub_img_size
                    )
                    display(img)
                    
            elif self.analysis.overlap_results is not None:
                # Show overlapping building blocks between current and new cutoff
                for position, num_overlap, overlap_blocks in self.analysis.overlap_results:
                    if not overlap_blocks:
                        print(f"No overlapping building blocks at position {position}")
                        continue
                        
                    print(f"Overlapping Building Blocks at Position {position} ({num_overlap} overlaps)")
                    # Create RDKit molecules for overlapping blocks
                    mols = []
                    legends = []
                    for block in overlap_blocks:
                        if block in self.analysis.smiles_dict:
                            mol = Chem.MolFromSmiles(self.analysis.smiles_dict[block])
                            if mol is not None:
                                mols.append(mol)
                                legends.append(f"{block}")
                    
                    if not mols:
                        print(f"No valid SMILES found for overlapping blocks at position {position}")
                        continue
                    
                    # Create and display the grid image
                    img = Draw.MolsToGridImage(
                        mols, 
                        legends=legends,
                        molsPerRow=mols_per_row, 
                        subImgSize=sub_img_size
                    )
                    display(img)
            else:
                raise ValueError("Must either provide a comparison_analysis or run check_overlap before visualizing overlapping blocks")
                
        else:
            # Show top building blocks
            if self.analysis.position_counters is None:
                raise ValueError("Must call find_top_building_blocks before visualization")
                
            for i, counter in enumerate(self.analysis.position_counters):
                print(f"Top 20 Building Blocks for Position {i + 1}")
                # Get the top 20 building blocks
                top_blocks = counter.most_common(20)
                
                # Create RDKit molecules for the top blocks
                mols = []
                legends = []
                for block, freq in top_blocks:
                    if block in self.analysis.smiles_dict:
                        mol = Chem.MolFromSmiles(self.analysis.smiles_dict[block])
                        if mol is not None:
                            mols.append(mol)
                            # Add building block name and frequency on the first line, fraction on the second line
                            legends.append(f"{block} (Freq: {freq})\nFraction: {round((freq / self.analysis.total_molecules) * 100, 2)}%")
                
                if not mols:
                    print(f"No valid SMILES found for top blocks at position {i + 1}")
                    continue
                
                # Create and display the grid image
                img = Draw.MolsToGridImage(
                    mols, 
                    legends=legends,
                    molsPerRow=mols_per_row, 
                    subImgSize=sub_img_size
                )
                display(img)

class TS_Benchmarks:
    """
    A class to generate visualizations of TS results. The goal here is to compare different cycles of TS runs with the same search strategy.
    This is mainly used to compare different search strategies to ground truth values.
    Each search strategy is used a no of times equal to the number of cycles.
    It is recommended that random baseline data and brute-force (exhaustive search) data be included for comparison.
    Reference data is optional, but is required for generating the barplots. A strip plot can be generated without reference data for comparison of methods.
    
    All required data (TS runs data, bar plot data, line plot data, and grouped statistics) is automatically generated during initialization.
    After creating an instance, you can directly call the plotting methods without additional data generation steps.
    """
    def __init__(self, no_of_cycles: int, methods_list: List[str], TS_runs_data: Dict[str, list], 
                 reference_data: Optional[pl.DataFrame] = None, top_n: int = 100, 
                 sort_type: str = "minimize", top_ns: Optional[List[int]] = None):
        """
        Initialize the TS_Benchmarks class and automatically generate all required data.

        Parameters:
        -----------
        no_of_cycles: int
            Number of cycles to run
        methods_list: List[str]
            List of method names (Names of the search strategies used for TS runs)
        TS_runs_data: Dict[str, list]
            Dictionary of TS runs, where the keys are the method names and the values are lists of different TS instances
        reference_data: Optional[pl.DataFrame]
            Reference data to compare against, this is the ground truth data
        top_n: int
            Number of top products to consider for bar plot analysis (default: 100)
        sort_type: str
            Type of sorting to perform ("minimize" or "maximize", default: "minimize")
        top_ns: Optional[List[int]]
            List of top N values for line plot analysis (default: [50, 100, 200, 300, 400, 500])
        """
        self.no_of_cycles = no_of_cycles
        self.methods_list = methods_list
        self.TS_runs_data = TS_runs_data
        self.reference_data = reference_data
        self.top_n = top_n
        self.sort_type = sort_type
        self.top_ns = top_ns if top_ns is not None else [50, 100, 200, 300, 400, 500]
        
        # Automatically generate all required data during initialization
        print("🔄 Initializing TS_Benchmarks and generating data...")
        self._generate_all_data()

    def _generate_all_data(self):
        """
        Generate all required data for plotting during initialization.
        This includes TS runs data, barplot data, and line plot performance data.
        """
        # Generate basic TS runs data
        self.combined_df_top_n, self.combined_df_all = self.gen_TS_runs_data(
            top_n=self.top_n, sort_type=self.sort_type
        )
        
        # Generate barplot data if reference data is available
        if self.reference_data is not None:
            self.bar_plot_df = self.get_barplot_TS_results_data(top_n=self.top_n)
            
            # Generate line plot performance data and grouped statistics
            self.line_plot_df = self.gen_line_plot_performance_data(top_ns=self.top_ns)
            self._generate_line_plot_grouped_stats()
        else:
            print("⚠️  Reference data not provided - barplot and line plot data will not be generated")
            self.bar_plot_df = None
            self.line_plot_df = None
            self.grouped_stats = None
        
        print("✅ Data generation completed! All plotting methods are now ready to use.")

    def _generate_line_plot_grouped_stats(self):
        """
        Generate grouped statistics for line plot with error bars.
        This calculates mean, std, upper, and lower bounds across cycles for each method and top_n.
        """
        if self.line_plot_df is None:
            print("⚠️  Line plot data not available - skipping grouped statistics generation")
            return
        
        # Calculate mean and std across cycles for error bars
        grouped_stats = self.line_plot_df.group_by(["method", "top_n"]).agg([
            pl.col("frac_top_n").mean().alias("mean"),
            pl.col("frac_top_n").std(ddof=1).alias("std"),  # Use sample std with ddof=1
            pl.col("frac_top_n").count().alias("n_cycles")
        ])
        
        # Handle cases where we might have only 1 cycle (std would be null)
        grouped_stats = grouped_stats.with_columns([
            pl.col("std").fill_null(0.0)
        ])
        
        # Add calculated upper and lower bounds for error bars
        grouped_stats = grouped_stats.with_columns([
            (pl.col("mean") + pl.col("std")).alias("upper"),
            (pl.col("mean") - pl.col("std")).alias("lower")
        ])
        
        # Store the grouped statistics and related data in the class
        self.grouped_stats = grouped_stats
        self.unique_top_ns = sorted(grouped_stats["top_n"].unique().to_list())
        self.actual_methods = sorted(grouped_stats["method"].unique().to_list())
        
        # Create cap data for error bars
        cap_width = (max(self.unique_top_ns) - min(self.unique_top_ns)) * 0.015  # 1.5% of x-axis range
        self.grouped_stats_caps = grouped_stats.with_columns([
            (pl.col("top_n") - cap_width).alias("cap_left"),
            (pl.col("top_n") + cap_width).alias("cap_right")
        ])
        self.cap_width = cap_width
        
        print(f"📊 Generated grouped statistics: {grouped_stats.shape} (mean, std, upper, lower)")

    def _get_color_scheme(self, include_ref: bool = True):
        """
        Generate a consistent color scheme for all plots.
        
        Parameters:
        -----------
        include_ref : bool
            Whether to include 'ref' in the domain for reference data
            
        Returns:
        --------
        alt.Scale
            Altair color scale with consistent colors across all plots
        """
        if include_ref and self.reference_data is not None:
            domain = self.methods_list + ["ref"]
        else:
            domain = self.methods_list
            
        return alt.Scale(
            domain=domain,
            range=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"][:len(domain)]
        )

    def gen_TS_runs_data(self, top_n: int = 100, sort_type: str = "minimize"):
        """
        Generates a single dataframe with all the TS runs data. This is a concatenated polars dataframe for each method

        Parameters:
        -----------
        top_n: int
            Number of top products to consider for each method
        sort_type: str
            Type of sorting to perform on the data
            "minimize" - sorts the data in ascending order
            "maximize" - sorts the data in descending order

        """
        # Initialize dictionary of dataframes
        all_ts_runs_data = {}
        all_ts_runs_data_top_n = {}
        for method in self.methods_list:
            ts_data = self.TS_runs_data[method] # Extract list of dataframes for each method
            method_dfs = [] # List to collect dataframes for this method
            for cycle in range(0, self.no_of_cycles):
                # Use cycle-1 as index since Python lists are 0-indexed but cycles are 1-indexed
                ts_data_temp = ts_data[cycle] # Extract dataframe for each cycle
                # Add cycle and method columns
                ts_data_temp = ts_data_temp.with_columns([
                    pl.lit(method).alias("method"),
                    pl.lit(str(cycle+1)).alias("cycle")
                ])
                # Only drop SMILES if it exists
                if "SMILES" in ts_data_temp.columns:
                    ts_data_temp = ts_data_temp.drop("SMILES")
                method_dfs.append(ts_data_temp)
                
            # Apply top_n filtering to each cycle separately - get top_n from EACH cycle
            if sort_type == "minimize":
                # This for the top n compounds found by each method
                method_dfs_filtered_top_n = [df.sort("score", descending=False).head(top_n) for df in method_dfs]
                # This for all the compounds found by each method in each cycle
                method_dfs_filtered_all = [df.sort("score", descending=False) for df in method_dfs]
                # Extract the top n compounds from the reference data
                if self.reference_data is not None:
                    ref_df_top_n = self.reference_data.sort("score", descending=False).head(top_n)
                    ref_df_all = self.reference_data.sort("score", descending=False) # Use all reference data
            elif sort_type == "maximize":
                # This for the top n compounds found by each method
                method_dfs_filtered_top_n = [df.sort("score", descending=True).head(top_n) for df in method_dfs]
                # This for all the compounds found by each method in each cycle
                method_dfs_filtered_all = [df.sort("score", descending=True) for df in method_dfs]
                # Extract the top n compounds from the reference data
                if self.reference_data is not None:
                    ref_df_top_n = self.reference_data.sort("score", descending=True).head(top_n)
                    ref_df_all = self.reference_data.sort("score", descending=True) # Use all reference data
            else:
                raise ValueError(f"Invalid sort_type: {sort_type}")
            
            # Concatenate the filtered cycles for this method
            ts_data_concat_top_n = pl.concat(method_dfs_filtered_top_n, how="vertical")
            ts_data_concat_all = pl.concat(method_dfs_filtered_all, how="vertical")
            all_ts_runs_data[method] = ts_data_concat_all # Add to dictionary
            all_ts_runs_data_top_n[method] = ts_data_concat_top_n # Add to dictionary
            
            # Track method statistics
            if len(ts_data_concat_all) > 0:
                method_min = ts_data_concat_all["score"].min()
                method_max = ts_data_concat_all["score"].max()
                print(f"  {method}: {len(ts_data_concat_all)} compounds, score range: {method_min:.3f} to {method_max:.3f}")
            
        self.all_ts_runs_data = all_ts_runs_data # Contains all the concatenated dataframes (by cycle) for each method
        
        # Combine all the dataframes into a single dataframe
        combined_df_top_n = pl.concat([all_ts_runs_data_top_n[method] for method in self.methods_list], how="vertical")
        combined_df_all = pl.concat([all_ts_runs_data[method] for method in self.methods_list], how="vertical")
        
        # Add reference data if provided
        if self.reference_data is not None:
            # Ensure reference data has proper cycle and method columns
            ref_df_top_n = ref_df_top_n.with_columns([
                pl.lit("ref").alias("method"),
                pl.lit("ref").alias("cycle")  # Reference data gets "ref" as cycle
            ])
            # Only drop SMILES if it exists
            if "SMILES" in ref_df_top_n.columns:
                ref_df_top_n = ref_df_top_n.drop("SMILES")
            
            # Process ref_df_all the same way as ref_df_top_n
            ref_df_all = ref_df_all.with_columns([
                pl.lit("ref").alias("method"),
                pl.lit("ref").alias("cycle")  # Reference data gets "ref" as cycle
            ])
            # Only drop SMILES if it exists
            if "SMILES" in ref_df_all.columns:
                ref_df_all = ref_df_all.drop("SMILES")
            
            combined_df_top_n = pl.concat([combined_df_top_n, ref_df_top_n], how="vertical")
            combined_df_all = pl.concat([combined_df_all, ref_df_all], how="vertical")
        
        # Filter out any invalid cycles (should only be 1-10 or "ref")
        valid_cycles = [str(i) for i in range(1, self.no_of_cycles + 1)] + ["ref"]
        combined_df_top_n = combined_df_top_n.filter(pl.col("cycle").is_in(valid_cycles))
        combined_df_all = combined_df_all.filter(pl.col("cycle").is_in(valid_cycles))
        
        # Convert method to categorical
        combined_df_top_n = combined_df_top_n.with_columns(pl.col("method").cast(pl.Categorical))
        combined_df_all = combined_df_all.with_columns(pl.col("method").cast(pl.Categorical))
        self.combined_df_top_n = combined_df_top_n # Store the combined dataframe in the class
        self.combined_df_all = combined_df_all # Store the combined dataframe in the class
        
        # Print comprehensive summary
        print(f"\n📊 Data Generation Summary:")
        print(f"  - Top-N dataset: {combined_df_top_n.shape} (rows × columns)")
        print(f"  - Complete dataset: {combined_df_all.shape} (rows × columns)")
        print(f"  - Cycles analyzed: {sorted([c for c in combined_df_top_n['cycle'].unique().to_list() if c != 'ref'])}")
        print(f"  - Methods: {sorted([m for m in combined_df_top_n['method'].unique().to_list() if m != 'ref'])}")
        if self.reference_data is not None:
            ref_count = len(combined_df_top_n.filter(pl.col("method") == "ref"))
            print(f"  - Reference compounds included: {ref_count}")
        
        return combined_df_top_n, combined_df_all
    
    def stripplot_TS_results(self, width: Optional[int] = None, height: Optional[int] = None, 
                         save_path: Optional[str] = None, show_plot: bool = True):
        """
        Generate a stripplot for TS results using altair. 
        This visualizes the distribution of scores across cycles and methods.
        Data is automatically generated during class initialization.

        Parameters:
        -----------
        width: Optional[int]
            Width of the plot
        height: Optional[int]
            Height of the plot
        save_path: Optional[str]
            Path to save the plot
        show_plot: bool
            If True, shows the plot in Jupyter
            
        Returns:
        --------
        altair.Chart or None
            The altair chart if show_plot is True, None otherwise
        """

         # Calculate dynamic dimensions if not provided
        if width is None:
            width = 400 + (len(self.methods_list) * 100)
        
        if height is None:
            height = 300 + (len(self.methods_list) * 50)
        
        # Calculate automatic y-axis scale based on data distribution
        score_min = self.combined_df_top_n["score"].min()
        score_max = self.combined_df_top_n["score"].max()
        score_range = score_max - score_min
        
        # Add 10% padding to both ends for better visualization
        y_min = score_min - (score_range * 0.1)
        y_max = score_max + (score_range * 0.1)
        
        # Use the standardized color scheme for consistency across all plots
        color_scheme = self._get_color_scheme(include_ref=True)
        
        # Create a stripplot for TS results using altair
        # Single plot with methods grouped within each cycle using x-axis positioning
        stripplot = alt.Chart(self.combined_df_top_n).mark_circle(
            size=40,
            opacity=0.7
        ).encode(
            x=alt.X("cycle:O", 
                   title="Cycle",
                   axis=alt.Axis(
                       labelAngle=0,
                       labelFontSize=16,
                       titleFontSize=18,
                       titlePadding=10
                   )),
            y=alt.Y("score:Q", 
                   title="Score",
                   scale=alt.Scale(domain=[y_min, y_max]),
                   axis=alt.Axis(
                       labelFontSize=16,
                       titleFontSize=18,
                       titlePadding=10
                   )),
            color=alt.Color("method:N", 
                           title="Method",
                           scale=color_scheme,
                           legend=alt.Legend(
                               orient="right",
                               titleFontSize=20,
                               labelFontSize=18,
                               titlePadding=10,
                               symbolSize=100
                           )),
            # Use xOffset to separate methods horizontally within each cycle
            xOffset=alt.XOffset("method:N", 
                               scale=alt.Scale(
                                   type="band", 
                                   paddingInner=0.3,
                                   paddingOuter=0.1
                               )),
            # Add small vertical jitter to separate overlapping points within each method
            yOffset=alt.YOffset("jitter:Q", scale=alt.Scale(range=[-3, 3]))
        ).transform_calculate(
            # Generate small random jitter for y-axis to separate overlapping points
            jitter="random()"
        ).properties(
            width=width,
            height=height
        )
        # Save the plot if save_path is provided
        if save_path:
            if save_path.endswith('.html'):
                stripplot.save(save_path)
            elif save_path.endswith('.png') or save_path.endswith('.svg'):
                stripplot.save(save_path, scale_factor=2.0)  # Higher resolution for images
            else:
                # Default to HTML if no extension specified
                stripplot.save(save_path + '.html')
        
        # Display in Jupyter if requested
        if show_plot:
            return stripplot
        else:
            return None
        
    def get_barplot_TS_results_data(self, top_n: int = 100):
        """
        Generate data for bar plot, for checking the number of hits recovered by each search strategy.
        To use this plot, you must have the reference compounds that serve as the ground truth.
        This shows what fraction of the top N reference compounds each method finds.

        Parameters:
        -----------
        top_n: int
            Number of top products to consider for each method. Ensure that the top_n is the same for all methods.
        Returns:
        --------
        bar_plot_df: polars DataFrame
            Dataframe with the number of hits found by each method in each cycle compared to the reference method
        """
        if self.reference_data is None:
            raise ValueError("Please ensure that reference_data is provided")
        
        # Get top N reference compounds (sorted by score to get the actual top compounds)
        top_ref_compounds = self.reference_data.sort("score", descending=False).head(top_n)["Name"].unique().to_list()
        
        # Filter combined data to exclude the reference method itself (to avoid duplication)
        # Only look at TS methods to see how many reference compounds they found
        ts_data = self.combined_df_top_n.filter(
            (pl.col("Name").is_in(top_ref_compounds)) & 
            (pl.col("method") != "ref")  # Exclude reference data to avoid duplication
        )
        
        # Count the number of hits found by each TS method in each cycle
        cycle_counts = ts_data.group_by(["cycle", "method"]).agg(
            pl.count().alias("found")
        )
        
        # Count unique hits per method across all cycles (for "concat" bars)
        method_totals = ts_data.group_by("method").agg(
            pl.col("Name").unique().count().alias("found")
        ).with_columns(pl.lit("concat").alias("cycle")).select(["cycle", "method", "found"])
        
        # Add reference baseline showing the total top N compounds as the benchmark
        ref_baseline = pl.DataFrame({
            "cycle": ["ref"], 
            "method": ["ref"], 
            "found": [top_n]  # This represents the total available top compounds
        })
        
        # Combine all data - ensure consistent data types
        cycle_counts_str = cycle_counts.with_columns([
            pl.col("method").cast(pl.String),
            pl.col("found").cast(pl.Int64)
        ])
        method_totals_str = method_totals.with_columns([
            pl.col("method").cast(pl.String),
            pl.col("found").cast(pl.Int64)
        ])
        ref_baseline_typed = ref_baseline.with_columns([
            pl.col("method").cast(pl.String),
            pl.col("found").cast(pl.Int64)
        ])
        
        bar_plot_df = pl.concat([cycle_counts_str, method_totals_str, ref_baseline_typed])
        
        # Convert method to categorical with specific order
        # Note: Using pl.Enum instead of pl.Categorical.cat.set_ordering for polars compatibility
        all_methods = self.methods_list + ["ref"] if "ref" not in self.methods_list else self.methods_list
        bar_plot_df = bar_plot_df.with_columns(
            pl.col("method").cast(pl.Enum(all_methods))
        )
        
        # Print summary of barplot data generation
        print(f"\n📊 Barplot Data Summary:")
        print(f"  - Reference baseline: Top {top_n} compounds from reference data")
        print(f"  - Analysis: How many reference compounds each method found")
        print(f"  - Data shape: {bar_plot_df.shape} (rows × columns)")
        cycle_counts = bar_plot_df.filter(pl.col("cycle") != "concat").filter(pl.col("cycle") != "ref")
        if len(cycle_counts) > 0:
            avg_found = cycle_counts["found"].mean()
            print(f"  - Average compounds found per cycle: {avg_found:.1f}")
        
        self.bar_plot_df = bar_plot_df  # Store the dataframe in the class
        return bar_plot_df
        

    def plot_barplot_TS_results(self, width: Optional[int] = None, height: Optional[int] = None,
                            save_path: Optional[str] = None, show_plot: bool = True):
        """
        Generate a barplot for TS results using altair. 
        This visualizes the number of reference hits recovered by each search strategy.
        Data is automatically generated during class initialization.
        
        Parameters:
        -----------
        width : Optional[int]
            Width of the plot in pixels
        height : Optional[int]
            Height of the plot in pixels
        save_path : Optional[str]
            Path to save the plot
        show_plot : bool
            If True, shows the plot in Jupyter
            
        Returns:
        --------
        altair.Chart or None
            The altair chart if show_plot is True, None otherwise
        """
        if width is None:
            width = max(400, len(self.bar_plot_df["cycle"].unique()) * 120)
        
        if height is None:
            height = 400
        
        # Use the standardized color scheme for consistency across all plots
        color_scheme = self._get_color_scheme(include_ref=True)
        
        # Create grouped barplot (not stacked) for better readability
        barplot = alt.Chart(self.bar_plot_df).mark_bar(
            stroke='white',
            strokeWidth=1
        ).encode(
            x=alt.X("cycle:O", 
                   title="Cycle",
                   axis=alt.Axis(
                       labelAngle=0,
                       labelFontSize=18,
                       titleFontSize=20
                   )),
            y=alt.Y("found:Q", 
                   title="Number of Top Reference Compounds Found",
                   axis=alt.Axis(
                       labelFontSize=18,
                       titleFontSize=20
                   )),
            color=alt.Color("method:N", 
                           title="Method",
                           scale=color_scheme,
                           legend=alt.Legend(
                               orient="right",
                               titleFontSize=20,
                               labelFontSize=18,
                               symbolSize=100,
                               padding=10
                           )),
            xOffset=alt.XOffset("method:N"),
            tooltip=["cycle:O", "method:N", "found:Q"]
        ).properties(
            width=width,
            height=height
        )
        
        # Add text labels positioned correctly on top of each bar
        text = alt.Chart(self.bar_plot_df).mark_text(
            align='center',
            baseline='bottom',
            fontSize=11,
            fontWeight='bold',
            dy=-5,
            color='black'
        ).encode(
            x=alt.X("cycle:O"),
            xOffset=alt.XOffset("method:N"),
            y=alt.Y("found:Q"),
            text=alt.condition(
                alt.datum.found > 0,  # Only show text if found > 0
                alt.Text("found:Q"),
                alt.value("")
            )
        )
        
        # Combine bar chart and text labels
        final_chart = alt.layer(barplot, text).resolve_scale(
            y='shared'
        )
        
        # Save the plot if save_path is provided
        if save_path:
            if save_path.endswith('.html'):
                final_chart.save(save_path)
            elif save_path.endswith('.png') or save_path.endswith('.svg'):
                final_chart.save(save_path, scale_factor=2.0)
            else:
                final_chart.save(save_path + '.html')
        
        # Display in Jupyter if requested
        if show_plot:
            return final_chart
        else:
            return None

    def gen_line_plot_performance_data(self, top_ns: List[int] = None):
        """
        Generate data for line plot, for checking the performance of each method.
        The plot looks at how the performance of the methods changes with comparison to the top N compounds found by the reference method.
        Efficiently calculates fraction of hits found for each top_n cutoff.

        Parameters:
        -----------
        top_ns : List[int], optional
            List of top N values to test (e.g., [50, 100, 200, 300, 400, 500])
            If None, defaults to [50, 100, 200, 300, 400, 500]
            
        Returns:
        --------
        line_plot_df : polars DataFrame
            DataFrame with columns: ['cycle', 'top_n', 'method', 'frac_top_n']
        """
        if self.reference_data is None:
            raise ValueError("Please ensure that reference_data is provided for performance analysis")
        
        if top_ns is None:
            top_ns = [50, 100, 200, 300, 400, 500]
        
        print("\n📈 Generating Line Plot Performance Data...")
        
        # Pre-calculate reference compound sets for each top_n cutoff
        ref_sorted = self.reference_data.sort("score", descending=False)
        ref_sets = {}
        for n in top_ns:
            ref_sets[n] = set(ref_sorted.head(n)["Name"].to_list())
        
        performance_data = []
        
        # Process each method
        for method in self.methods_list:
            print(f"  Processing method: {method}")
            
            # Get all compounds found by this method across all cycles - once per method
            method_data = self.combined_df_all.filter(pl.col("method") == method)
            
            # Process individual cycles
            for cycle in range(1, self.no_of_cycles + 1):
                cycle_id = str(cycle)
                
                # Get compounds found by this method in this specific cycle
                cycle_compounds = set(
                    method_data.filter(pl.col("cycle") == cycle_id)["Name"].to_list()
                )
                
                # Calculate fraction for each top_n value using set intersection
                for n in top_ns:
                    ref_set = ref_sets[n]
                    hits_found = len(cycle_compounds.intersection(ref_set))
                    frac_top_n = hits_found / n
                    
                    performance_data.append({
                        "cycle": cycle_id,
                        "top_n": n,
                        "method": method,
                        "frac_top_n": frac_top_n
                    })
        
        # Convert to polars DataFrame
        line_plot_df = pl.DataFrame(performance_data)
        
        # Ensure consistent data types
        line_plot_df = line_plot_df.with_columns([
            pl.col("cycle").cast(pl.String),
            pl.col("top_n").cast(pl.Int64),
            pl.col("method").cast(pl.String),
            pl.col("frac_top_n").cast(pl.Float64)
        ])
        
        # Print comprehensive summary
        print(f"\n📊 Line Plot Data Summary:")
        print(f"  - Data shape: {line_plot_df.shape} (rows × columns)")
        print(f"  - Top-N values tested: {top_ns}")
        print(f"  - Methods analyzed: {self.methods_list}")
        print(f"  - Cycles per method: {self.no_of_cycles}")
        
        # Show performance range
        if len(line_plot_df) > 0:
            min_frac = line_plot_df["frac_top_n"].min()
            max_frac = line_plot_df["frac_top_n"].max()
            print(f"  - Performance range: {min_frac:.3f} to {max_frac:.3f} (fraction found)")
        print("✅ Line plot data generation completed successfully!")
        
        self.line_plot_df = line_plot_df  # Store in the class
        return line_plot_df
    
    def plot_line_performance_with_error_bars(self, width: Optional[int] = None, height: Optional[int] = None,
                                            save_path: Optional[str] = None, show_plot: bool = True):
        """
        Generate a line plot with error bars for method performance using altair.
        Shows mean fraction of reference compounds found across cycles with standard deviation error bars.
        Data and grouped statistics are automatically generated during class initialization.
        
        Parameters:
        -----------
        width : Optional[int]
            Width of the plot in pixels
        height : Optional[int]
            Height of the plot in pixels
        save_path : Optional[str]
            Path to save the plot
        show_plot : bool
            If True, shows the plot in Jupyter
            
        Returns:
        --------
        altair.Chart or None
            The altair chart if show_plot is True, None otherwise
        """
        if not hasattr(self, 'grouped_stats') or self.grouped_stats is None:
            raise ValueError("Grouped statistics not available. This should have been generated during initialization.")
        
        if width is None:
            width = 800
        
        if height is None:
            height = 500
        
        # Use pre-generated grouped statistics
        grouped_stats = self.grouped_stats
        
        # Use the standardized color scheme for consistency across all plots
        # Line plot doesn't include reference data, so set include_ref=False
        color_scheme = self._get_color_scheme(include_ref=False)
        
        # Create the base chart
        base = alt.Chart(grouped_stats)
        
        # Main line plot with thicker lines and larger points
        line_plot = base.mark_line(
            point=alt.OverlayMarkDef(size=120, filled=True),
            strokeWidth=4
        ).encode(
            x=alt.X("top_n:Q", 
                   title="Top N Compounds",
                   scale=alt.Scale(domain=[min(self.unique_top_ns), max(self.unique_top_ns)]),
                   axis=alt.Axis(
                       labelFontSize=24,
                       titleFontSize=28,
                       labelAngle=0,
                       values=self.unique_top_ns,  # Explicitly set tick values
                       format="d"  # Format as integers
                   )),
            y=alt.Y("mean:Q", 
                   title="Mean Fraction Found",
                   scale=alt.Scale(domain=[0, 1]),
                   axis=alt.Axis(
                       labelFontSize=24,
                       titleFontSize=28,
                       format=".1%"
                   )),
            color=alt.Color("method:N", 
                           title="Method",
                           scale=color_scheme,
                           legend=alt.Legend(
                               orient="right",
                               titleFontSize=20,
                               labelFontSize=18,
                               symbolSize=100
                           )),
            tooltip=["method:N", "top_n:O", "mean:Q", "std:Q", "n_cycles:O"]
        )
        
        # Error bars using mark_rule for vertical lines
        error_bars = base.mark_rule(
            strokeWidth=3,
            opacity=0.8
        ).encode(
            x=alt.X("top_n:Q"),
            y=alt.Y("lower:Q"),
            y2=alt.Y2("upper:Q"),
            color=alt.Color("method:N", scale=color_scheme, legend=None)
        )
        
        # Error bar caps (horizontal lines at top and bottom) using mark_rule with x2 for horizontal lines
        # Use pre-generated cap data
        grouped_stats_caps = self.grouped_stats_caps
        base_caps = alt.Chart(grouped_stats_caps)
        
        error_caps_top = base_caps.mark_rule(
            strokeWidth=3,
            opacity=0.8
        ).encode(
            x=alt.X("cap_left:Q"),
            x2=alt.X2("cap_right:Q"),
            y=alt.Y("upper:Q"),
            color=alt.Color("method:N", scale=color_scheme, legend=None)
        )
        
        error_caps_bottom = base_caps.mark_rule(
            strokeWidth=3,
            opacity=0.8
        ).encode(
            x=alt.X("cap_left:Q"),
            x2=alt.X2("cap_right:Q"),
            y=alt.Y("lower:Q"),
            color=alt.Color("method:N", scale=color_scheme, legend=None)
        )
        
        # Create line plot separately to ensure legend shows
        line_only = line_plot.properties(
            width=width,
            height=height,
            title="Line Plot Only (Testing Legend)"
        )
        
        # Create the full chart with error bars
        final_chart = alt.layer(
            error_bars, 
            error_caps_top, 
            error_caps_bottom, 
            line_plot
        ).resolve_scale(
            color='independent'
        ).properties(
            width=width,
            height=height,
            title=alt.TitleParams(
                text=f"Mean Top N Fraction Found Across {self.no_of_cycles} Cycles",
                fontSize=22,
                anchor="start"
            )
        )
        
        # Store chart components for potential later access
        self.line_only_chart = line_only
        self.final_chart = final_chart
        
        # Save the plot if save_path is provided
        if save_path:
            if save_path.endswith('.html'):
                final_chart.save(save_path)
            elif save_path.endswith('.png') or save_path.endswith('.svg'):
                final_chart.save(save_path, scale_factor=2.0)
            else:
                final_chart.save(save_path + '.html')
        
        # Display in Jupyter if requested
        if show_plot:
            return final_chart
        else:
            return grouped_stats

    def get_performance_summary(self):
        """
        Get a summary of the performance data generated by plot_line_performance_with_error_bars.
        
        Returns:
        --------
        dict
            Dictionary containing all stored performance data and statistics
        """
        if not hasattr(self, 'grouped_stats') or self.grouped_stats is None:
            raise ValueError("Please run plot_line_performance_with_error_bars() first to generate performance data")
        
        return {
            'grouped_stats': self.grouped_stats,
            'grouped_stats_caps': self.grouped_stats_caps if hasattr(self, 'grouped_stats_caps') else None,
            'unique_top_ns': self.unique_top_ns if hasattr(self, 'unique_top_ns') else None,
            'actual_methods': self.actual_methods if hasattr(self, 'actual_methods') else None,
            'cap_width': self.cap_width if hasattr(self, 'cap_width') else None,
            'color_scheme': self.color_scheme if hasattr(self, 'color_scheme') else None,
            'line_plot': self.line_plot if hasattr(self, 'line_plot') else None,
            'error_bars': self.error_bars if hasattr(self, 'error_bars') else None,
            'final_chart': self.final_chart if hasattr(self, 'final_chart') else None
        }



