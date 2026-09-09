"""Altair benchmark visualisations for Thompson Sampling recovery runs.

Requires the optional ``viz`` extra (``pip install chem-tactics[viz]``).
"""

from typing import Dict, List, Optional

import polars as pl
import altair as alt


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
                         save_path: Optional[str] = None, show_plot: bool = True,
                         legend_position: str = "right"):
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
        legend_position: str
            Position of the legend. "right" (default) or "bottom" for horizontal legend below plot.

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

        # Generate proper sort order for cycles (1, 2, 3, ..., 10, ref) instead of lexicographic (1, 10, 2, ...)
        cycle_order = [str(i) for i in range(1, self.no_of_cycles + 1)]
        if self.reference_data is not None:
            cycle_order.append("ref")

        # Build legend config based on position
        if legend_position == "bottom":
            legend_config = alt.Legend(
                orient="bottom",
                direction="horizontal",
                titleFontSize=16,
                labelFontSize=14,
                columns=0  # Auto-wrap
            )
        else:
            legend_config = alt.Legend(
                orient="right",
                titleFontSize=20,
                labelFontSize=18,
                titlePadding=10,
                symbolSize=100
            )

        # Create a stripplot for TS results using altair
        # Single plot with methods grouped within each cycle using x-axis positioning
        stripplot = alt.Chart(self.combined_df_top_n).mark_circle(
            size=40,
            opacity=0.7
        ).encode(
            x=alt.X("cycle:O",
                   title="Cycle",
                   sort=cycle_order,
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
                           legend=legend_config),
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
        # Use sort_type to determine sort direction: "minimize" = ascending, "maximize" = descending
        _descending = self.sort_type == "maximize"
        top_ref_compounds = self.reference_data.sort("score", descending=_descending).head(top_n)["Name"].unique().to_list()
        
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
                            save_path: Optional[str] = None, show_plot: bool = True,
                            legend_position: str = "right", dark_mode: bool = False):
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
        legend_position : str
            Position of the legend. "right" (default) or "bottom" for horizontal legend below plot.
        dark_mode : bool
            If True, uses white text for bar labels (for dark backgrounds). Default is False (black text).

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

        # Generate proper sort order for cycles (1, 2, 3, ..., 10, ref) instead of lexicographic (1, 10, 2, ...)
        cycle_order = [str(i) for i in range(1, self.no_of_cycles + 1)]
        if self.reference_data is not None:
            cycle_order.append("ref")

        # Build legend config based on position
        if legend_position == "bottom":
            legend_config = alt.Legend(
                orient="bottom",
                direction="horizontal",
                titleFontSize=16,
                labelFontSize=14,
                columns=0  # Auto-wrap
            )
        else:
            legend_config = alt.Legend(
                orient="right",
                titleFontSize=20,
                labelFontSize=18,
                symbolSize=100,
                padding=10
            )

        # Create grouped barplot (not stacked) for better readability
        barplot = alt.Chart(self.bar_plot_df).mark_bar(
            stroke='white',
            strokeWidth=1
        ).encode(
            x=alt.X("cycle:O",
                   title="Cycle",
                   sort=cycle_order,
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
                           legend=legend_config),
            xOffset=alt.XOffset("method:N"),
            tooltip=["cycle:O", "method:N", "found:Q"]
        ).properties(
            width=width,
            height=height
        )

        # Text color based on dark mode
        text_color = 'white' if dark_mode else 'black'

        # Add text labels positioned correctly on top of each bar
        text = alt.Chart(self.bar_plot_df).mark_text(
            align='center',
            baseline='bottom',
            fontSize=11,
            fontWeight='bold',
            dy=-5,
            color=text_color
        ).encode(
            x=alt.X("cycle:O", sort=cycle_order),
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
        # Use sort_type to determine sort direction: "minimize" = ascending, "maximize" = descending
        _descending = self.sort_type == "maximize"
        ref_sorted = self.reference_data.sort("score", descending=_descending)
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
                                            save_path: Optional[str] = None, show_plot: bool = True,
                                            legend_position: str = "right"):
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
        legend_position : str
            Position of the legend. "right" (default) or "bottom" for horizontal legend below plot.

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

        # Build legend config based on position
        if legend_position == "bottom":
            legend_config = alt.Legend(
                orient="bottom",
                direction="horizontal",
                titleFontSize=16,
                labelFontSize=14,
                columns=0  # Auto-wrap
            )
        else:
            legend_config = alt.Legend(
                orient="right",
                titleFontSize=20,
                labelFontSize=18,
                symbolSize=100
            )

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
                           legend=legend_config),
            order=alt.Order("top_n:Q"),  # Ensures lines connect in ascending x order
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
