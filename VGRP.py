
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from sklearn import linear_model
from sklearn.metrics import r2_score
from PIL import Image, ImageTk
import os
import sys
import geopandas as gpd
from shapely.geometry import Point
from scipy.interpolate import griddata
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
import shapefile
from PIL import Image
from matplotlib.widgets import SpanSelector
from scipy.signal import savgol_filter
from tkinter import filedialog, simpledialog, messagebox
from sklearn.linear_model import LinearRegression
from scipy.interpolate import griddata
from pykrige.ok import OrdinaryKriging
from matplotlib.backend_bases import NavigationToolbar2
from scipy.interpolate import interp1d, PchipInterpolator
import xml.etree.ElementTree as ET
import serial
import serial.tools.list_ports
from joblib import dump
import datetime
import threading
import threading
import datetime
import json
import serial.tools.list_ports
import pandas as pd
from sklearn.linear_model import LinearRegression
from joblib import dump
from PIL import Image, ImageTk
from tkinter import ttk, messagebox, simpledialog
import serial.tools.list_ports
from PIL import Image, ImageTk
import datetime
import serial
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from tkinter import scrolledtext
import webbrowser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import json
import datetime
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.signal import filtfilt

class VGR_CalibrationToolbox:
    def __init__(self, root):
        self.root = root
        self.root.title("VGR - Multigas Calibration Toolbox (Open Science Edition)")
        self.root.geometry("1100x750")
        
        # Internal database for caching runtime data and calculated parameters
        self.raw_df = None
        self.results = {}
        self.results = {} # Dictionary to temporarily cache calculated calibration coefficients
        
        # Ensure required directories for models and results exist
        os.makedirs("models", exist_ok=True)
        os.makedirs("calibration_results", exist_ok=True)
        
        self.setup_ui()

    # ---------------- User Interface (UI) Setup ----------------
    def setup_ui(self):
        style = ttk.Style()
        style.configure("TButton", font=("Arial", 10))
        
        # Main application layout frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill="both", expand=True)

        header = ttk.Label(main_frame, text="Multigas Calibration & Processing", font=("Arial", 16, "bold"))
        header.pack(pady=10)

        # Control panel for Step 1: Data ingestion
        ctrl_frame = ttk.LabelFrame(main_frame, text=" Step 1: Data Loading ", padding="10")
        ctrl_frame.pack(fill="x", pady=5)

        self.load_btn = ttk.Button(ctrl_frame, text="Load Pre-registered Data (CSV/Excel)", command=self.load_file)
        self.load_btn.pack(side="left", padx=10)

        self.status_var = tk.StringVar(value="Status: Waiting for data...")
        ttk.Label(ctrl_frame, textvariable=self.status_var, foreground="blue").pack(side="left", padx=20)

        # Section for Phase 2 & 3: Plateau selection and linear regression
        calib_frame = ttk.LabelFrame(main_frame, text=" Step 2 & 3: Plateau Selection & Regression ", padding="10")
        calib_frame.pack(fill="x", pady=5)

        ttk.Label(calib_frame, text="Select Target Gas:").pack(side="left", padx=5)
        self.gas_selector = ttk.Combobox(calib_frame, values=["CO2", "SO2", "H2S"], width=10)
        self.gas_selector.pack(side="left", padx=5)
        self.gas_selector.set("SO2")

        self.process_btn = ttk.Button(calib_frame, text="Process & Calibrate", command=self.start_calibration_process, state="disabled")
        self.process_btn.pack(side="left", padx=10)

        self.save_btn = ttk.Button(calib_frame, text="Save to Archive (JSON/CSV)", command=self.save_results, state="disabled")
        self.save_btn.pack(side="left", padx=10)

        # Data preview components (Treeview widget)
        self.tree_frame = ttk.Frame(main_frame)
        self.tree_frame.pack(fill="both", expand=True, pady=10)
        
        self.tree = ttk.Treeview(self.tree_frame, show='headings')
        self.tree.pack(fill="both", expand=True)

    # ---------------- Data Loading & Ingestion Logic ----------------
    def load_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Data Files", "*.csv *.xlsx *.xls")])
        if file_path:
            try:
                if file_path.endswith('.csv'):
                    self.raw_df = pd.read_csv(file_path)
                else:
                    self.raw_df = pd.read_excel(file_path)
                
                # Initial preprocessing: Remove completely blank rows
                self.raw_df.dropna(how='all', inplace=True)
                
                self.status_var.set(f"Loaded: {os.path.basename(file_path)}")
                self.process_btn.config(state="normal")
                self.show_preview()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")

    def show_preview(self):
        """Render the first 15 rows of the loaded dataset inside the UI Treeview for verification."""
        self.tree["columns"] = list(self.raw_df.columns)
        for col in self.raw_df.columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=100)
        
        for i in self.tree.get_children(): 
            self.tree.delete(i)
        for _, row in self.raw_df.head(15).iterrows():
            self.tree.insert("", "end", values=list(row))

    def clean_pressure_transients(self, df):
        """Mathematically eliminate severe pressure transients using a 2-sigma threshold 
        to prevent regression bias in unstable environmental conditions."""
        if df is None: return None
        df_cleaned = df.copy()
        if 'Pressure' in df_cleaned.columns:
            df_cleaned['P_delta'] = df_cleaned['Pressure'].diff().abs()
            threshold = df_cleaned['P_delta'].std() * 2
            return df_cleaned[df_cleaned['P_delta'] < threshold].copy()
        return df_cleaned

    # ---------------------------------------------------------
# Phase 1: Multi-Point Sensor Calibration Management
    # ---------------------------------------------------------
    def start_calibration_process(self):
        """Execute the complete multi-step calibration workflow: inputting reference concentrations, triggering the span selector, and performing ordinary least squares (OLS) regression."""
        target_gas = self.gas_selector.get()
        num_points = simpledialog.askinteger(
            "Calibration Steps", 
            f"Enter number of concentration levels for {target_gas}:",
            initialvalue=2, minvalue=2, maxvalue=10
        )
        
        if not num_points: return

        # Apply the pressure transient filter if the cleaning routine is available
        cleaned_df = self.clean_pressure_transients(self.raw_df) if hasattr(self, 'clean_pressure_transients') else self.raw_df
        self.temp_points = [] 

        for i in range(num_points):
            ref_conc = simpledialog.askfloat(
                f"Level {i+1}", 
                f"Enter Reference Concentration (ppm) for point {i+1}:"
            )
            if ref_conc is None: break
            
            # Invoke the interactive visualization window for regional span selection
            self.open_span_selector_sync(cleaned_df, target_gas, ref_conc, i+1)
            
        if len(self.temp_points) >= 2:
            self.finalize_multi_point_regression(target_gas, self.temp_points)
        else:
            messagebox.showwarning("Incomplete Data", "Calibration requires at least 2 points.")

    # Phase 2: Interactive Regional Span Selector (Plateau Identification)
    # ---------------------------------------------------------
    def open_span_selector_sync(self, df, gas, ref_conc, step_num):
        """Open a blocking graphical window allowing the user to select stable gas signal plateaux 
        by dragging the mouse cursor across the timeseries chart."""
        plot_win = tk.Toplevel(self.root)
        plot_win.title(f"VGR Selector - Level {step_num}")
        plot_win.grab_set()

        fig, ax = plt.subplots(figsize=(8, 4), dpi=100)
        ax.plot(df[gas].values, color='#2c3e50', label=f"Raw {gas}")
        ax.set_title(f"Drag to select stable range for {ref_conc} ppm")
        ax.legend()

        self.selected_range = None # Internal variable to track bounding horizontal coordinates

        def onselect(xmin, xmax):
            # Extract boundary indices from coordinate floats
            start, end = int(xmin), int(xmax)
            self.selected_range = (start, end)
            
            # Graphically highlight the selected target region on the canvas
            for patch in ax.patches: patch.remove()
            ax.axvspan(xmin, xmax, color='red', alpha=0.3)
            canvas.draw()

        from matplotlib.widgets import SpanSelector
        span = SpanSelector(ax, onselect, 'horizontal', useblit=True, 
                            props=dict(alpha=0.5, facecolor='red'))

        canvas = FigureCanvasTkAgg(fig, master=plot_win)
        canvas.get_tk_widget().pack(fill="both", expand=True)

        def confirm():
            if self.selected_range:
                start, end = self.selected_range
                actual_start, actual_end = sorted([start, end])
                
                # Calculate the mean electrochemical response within the defined plateau window
                avg_signal = df[gas].iloc[actual_start:actual_end].mean()
                
                # Intelligent parsing for environmental reference pressure tracking
                if 'Pressure' in df.columns:
                    # Extract mean barometric pressure from the same temporal window
                    avg_p = df['Pressure'].iloc[actual_start:actual_end].mean()
                else:
                    # Prompt user manually if the barometric log column is absent
                    user_p = simpledialog.askfloat("Pressure Missing", 
                                                   "Column 'Pressure' not found in Excel.\nPlease enter the reference pressure (hPa):",
                                                   initialvalue=780)
                    # Adopt user-defined value or fallback to standard pressure (780)
                    avg_p = user_p if user_p is not None else 780
                
                # Append the parsed data point to the active multi-point matrix
                self.temp_points.append({'x': ref_conc, 'y': avg_signal, 'p': avg_p})
                plt.close(fig)
                plot_win.destroy()
            else:
                messagebox.showwarning("Warning", "Please drag on the chart to select a range!")

        tk.Button(plot_win, text="Confirm & Next", command=confirm, bg="#27ae60", fg="white").pack(pady=5)
        self.root.wait_window(plot_win)

    # Phase 3: Ordinary Least Squares (OLS) Regression and Parameter Archiving
    # ---------------------------------------------------------
    def finalize_multi_point_regression(self, gas, data_points):
        """Fit the parsed plateau data to an OLS linear model to calculate the final 
        sensor calibration coefficients: Gain (ai), Intercept (bi), and the R-squared value."""
        
        # 1. Structure arrays for the linear model (X: standard concentration, y: sensor response)
        X = np.array([p['x'] for p in data_points]).reshape(-1, 1) # Independent variable: Standard gas concentration (ppm)
        y = np.array([p['y'] for p in data_points])                # Dependent variable: Raw sensor response
        mean_p = np.mean([p['p'] for p in data_points])            # Averaged reference barometric pressure (hPa)

        # 2. Fit the ordinary least squares regression line
        model = LinearRegression().fit(X, y)
        gain_ai = model.coef_[0]
        offset_bi = model.intercept_
        r2 = model.score(X, y)

        # 3. Retrieve cross-interference calibration factors (strictly enforced for H2S sensors)
        # ---------------------------------------------------------
        aij_value = 0.0
        if gas == "H2S": 
            aij_value = simpledialog.askfloat("Cross Interference", 
                                            "Enter aij for H2S (SO2 interference):", 
                                            initialvalue=0.12)
            if aij_value is None: aij_value = 0.0
        else:
            # Cross-interference is default-zeroed for direct ultraviolet or non-interfered targets
            aij_value = 0.0
        # ---------------------------------------------------------

        # 4. Cache structured metadata and parameters for structural export
        self.results[gas] = {
            "gas_type": gas,
            "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "gain_ai": float(gain_ai),
            "offset_bi": float(offset_bi),
            "r_squared": float(r2),
            "calib_pressure": float(mean_p),
            "cross_interference": float(aij_value) # Archived cross-talk factor for dynamic plume correction
        }
        
        self.save_btn.config(state="normal")
        
        # Compile and render a comprehensive analytic summary report to the user
        report = (f"Calibration Report: {gas}\n"
                  f"{'-'*30}\n"
                  f"Gain (ai): {gain_ai:.6f}\n"
                  f"Offset (bi): {offset_bi:.4f}\n"
                  f"R-Squared: {r2:.4f}\n"
                  f"Ref Pressure: {mean_p:.1f} hPa\n"
                  f"Cross-Int (aij): {aij_value}")
        
        messagebox.showinfo("VGR Analytics", report)
        
    def save_results(self):
        """Archive calculated calibration coefficients to JSON and CSV formats 
        while mitigating NumPy-to-Python type conversion errors (TypeError)."""
        try:
            os.makedirs("models", exist_ok=True)
            os.makedirs("calibration_results", exist_ok=True)
            json_path = "models/calibration_coeffs.json"
            
            # Ingest existing structural database or initialize a clean memory map
            archive = {}
            if os.path.exists(json_path):
                with open(json_path, "r") as f:
                    try: archive = json.load(f)
                    except: archive = {}

            clean_entries = []
            for g, data in self.results.items():
                # Sanitize data types: Explicitly convert NumPy primitives to standard Python floats/ints
                safe_data = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in data.items()}
                if g not in archive: archive[g] = {"history": []}
                archive[g]["history"].append(safe_data)
                clean_entries.append(safe_data)
            
            with open(json_path, "w") as f:
                json.dump(archive, f, indent=4)

            if clean_entries:
                path = f"calibration_results/report_{datetime.datetime.now().strftime('%Y%m%d_%H%M')}.csv"
                pd.DataFrame(clean_entries).to_csv(path, index=False)
                messagebox.showinfo("Success", f"Saved to {path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))
            
selected_range = []
done = False

def onselect(xmin, xmax):
    global selected_range, done
    
    # Determine upper boundary limits dynamically to prevent out-of-bounds indexing exceptions
    try:
        max_limit = len(df) # نام متغیر دیتافریم یا لیست داده‌هایت را اینجا چک کن
    except NameError:
        max_limit = int(xmax)

    # Bounding constraint execution: Restrict index range strictly within [0, max_limit]
    start = max(0, min(int(xmin), max_limit))
    end = max(0, min(int(xmax), max_limit))

    # Enforce left-to-right positional sorting to handle inverted graphical selections
    actual_start, actual_end = sorted([start, end])

    selected_range = [actual_start, actual_end]
    done = True
    
    plt.close()
    
    # Trigger a non-blocking asynchronous user alert displaying the finalized boundary coordinates
    root.after(100, lambda: messagebox.showinfo(
        "Range Selected", 
        f"Selected range: {selected_range}"
    ))


class GasSmoother:
    def best_t90_filter(self, data, compare_gas, t90_range=(1, 60)):
     best_corr = -1
     best_t90 = 1
     best_smoothed = data.copy()

     for t90 in range(t90_range[0], t90_range[1] + 1):
         # Calculate the exponential smoothing factor (alpha) from the t90 response time
         alpha = 1 - np.exp(-2.2 / t90)
        
         # Define the numerator (b) and denominator (a) coefficients for the single-pole digital filter
         b = [alpha]
         a = [1, -(1 - alpha)]
        
         # Execute zero-phase forward-backward digital filtering to completely eliminate sensor response lag
         try:
             # The filtfilt routine applies the filter coefficients in both directions to eliminate temporal shifting
             smoothed = filtfilt(b, a, data)
             
             min_len = min(len(smoothed), len(compare_gas))
             corr = np.corrcoef(smoothed[:min_len], compare_gas[:min_len])[0, 1]
        
             if not np.isnan(corr) and corr > best_corr:
                 best_corr = corr
                 best_t90 = t90
                 best_smoothed = smoothed
         except:
             continue

     return best_smoothed, best_t90, best_corr

    def apply_selected_filter(self, data, compare_gas=None, selected_range=None, gas_name="this gas"):
        # Inject standard Central Moving Average to the interactive filter selection matrix
        filter_choice = simpledialog.askinteger(
            f"Smoothing Filter for {gas_name}",
            "Which filter would you like to apply?\n"
            "1. Savitzky-Golay\n"
            "2. Single pole Low-pass (t90)\n"
            "3. Both SG & t90\n"
            "4. Moving Average\n"
            "5. None (No filtering)",
            minvalue=1, maxvalue=5
        )

        if filter_choice is None or filter_choice == 5:
            return data.copy(), 0, len(data)

        final_data = data.copy()
        total_length = len(data)
        
        # --- Array Slicing Alignment: Coordinate boundary verification without programmatic alteration ---
        if selected_range is not None and len(selected_range) == 2:
            start, end = int(selected_range[0]), int(selected_range[1])
            # Bound slicing indices strictly within the actual timeline array length
            start = max(0, min(start, total_length - 1))
            end = max(1, min(end, total_length))
        else:
            # Fallback routine: Apply a conservative internal safe window if no boundary is specified
            start, end = 5, total_length - 5

        # Verify the segment fulfills the minimum statistical degrees of freedom (N >= 3)
        if end - start < 3:
            messagebox.showwarning("Warning", "Selected range is too short.")
            return final_data, start, end

        # 1. Single pole Low-pass (t90)
        if filter_choice in [2, 3]:
            if compare_gas is not None:
                g1 = data[start:end]
                g2 = compare_gas[start:end]
                min_len = min(len(g1), len(g2))
                if min_len >= 5:
                    smoothed_part, t90, corr = self.best_t90_filter(g1, g2)
                    final_data[start:start+min_len] = smoothed_part[:min_len]
                    messagebox.showinfo("t90 Result", f"Best t90: {t90}s\nCorr: {corr:.3f}")

        # 2. Savitzky-Golay
        if filter_choice in [1, 3]:
            window = simpledialog.askinteger("SG", "Window length (odd number >= 3):", minvalue=3)
            if window:
                if window % 2 == 0: window += 1
                poly = simpledialog.askinteger("SG", "Polynomial order:", minvalue=1)
                if poly and poly < window:
                    try:
                       # Compute Savitzky-Golay parameters strictly inside the regional slice bounds
                        final_data[start:end] = savgol_filter(
                            final_data[start:end],
                            window,
                            poly
                        )
                    except Exception as e:
                        messagebox.showerror("Error", f"SG Error: {e}")

        # 3. Moving Average 
        if filter_choice == 4:
            
            window_ma = simpledialog.askinteger("Moving Average", "Enter window size:", initialvalue=21, minvalue=2)
            
            if window_ma:
                try:
                    # Wrap the target numerical array into a Pandas Series object
                    series = pd.Series(final_data)
                    
                    # Compute a center-aligned moving average, matching spreadsheet engine specifications
                    ma_smoothed = series.rolling(window=window_ma,center=True, min_periods=1).mean()
                    
                    # Backfill and forward-fill remaining NaN edge values to prevent edge propagation
                    ma_smoothed = ma_smoothed.ffill().bfill()
                    
                    # Update the analytical array with the finalized smoothed series
                    final_data = ma_smoothed.to_numpy()
                    
                    messagebox.showinfo("Success", f"Moving Average (window={window_ma}) applied.")
                    
                    # Critical compliance fix: Return the 3-element tuple expected by the calibration routine
                    return final_data, start, end

                except Exception as e:
                    messagebox.showerror("Error", f"MA processing failed: {str(e)}")
                    # Exception handling: Return a copy of the raw segment to prevent core application crashes
                    return data.copy(), start, end

        return final_data, start, end

if hasattr(sys, '_MEIPASS'):
    base_path = sys._MEIPASS
else:
    base_path = os.path.abspath(".")

background_path = os.path.join(base_path, "backgr.png")
logo_path = os.path.join(base_path, "logo.png")

# Functions for gas ratio analysis
def analyze_h2s_so2():
    run_analysis("H2S", "SO2")
    
def analyze_co2_h2s():
    run_analysis("CO2", "H2S")

def analyze_h2o_co2():
    run_analysis("H2O", "CO2")

def analyze_h2o_h2s():
    run_analysis("H2O", "H2S")

def analyze_h2o_so2():
    run_analysis("H2O", "SO2")

def analyze_co2_so2():
    run_analysis("CO2", "SO2")

# Request initial range from user
def get_initial_range():
    start_row = simpledialog.askinteger("Input", "Enter start row for data processing:", minvalue=0, maxvalue=1000000)
    end_row = simpledialog.askinteger("Input", "Enter end row for data processing:", minvalue=start_row, maxvalue=1000000)
    return start_row, end_row

# Request delay correction range from user
def get_delay_correction_range():
    start_time = simpledialog.askinteger("Input", "Enter start time (in seconds) for delay correction:", minvalue=0, maxvalue=1000)
    end_time = simpledialog.askinteger("Input", "Enter end time (in seconds) for delay correction:", minvalue=start_time, maxvalue=1000)
    return start_time, end_time

# Load Excel or CSV file
def load_excel_file():
    global current_file_path
    current_file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx"), ("CSV files", "*.csv")])
    
    if current_file_path:  
        try:
            if current_file_path.endswith('.csv'):
                data = pd.read_csv(current_file_path)
            else:
                data = pd.read_excel(current_file_path, engine='openpyxl')
            return data
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")
            return None

def apply_calibration_smart(gas_name, gas_data, all_data_df=None):
    """
    1. Intelligent archive management (gracefully falls back to manual adjustment if the archive is missing).
    2. Redundant prompt elimination following archive-based calibration.
    3. Strict preservation of geochemical conversion formulae and cross-interference corrections.
    """
    # Step 0: Prompt user to determine whether archive-based calibration is required
    use_archive = messagebox.askyesno("Calibration", f"Does {gas_name} require calibration from archive?")
    
    if not use_archive:
        # If the user selects 'No', the routine bypasses the archive and proceeds directly to manual conversion/formula entry.
        return apply_cconversion(gas_name, gas_data)

    # 1. Load archive
    calib_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", "calibration_coeffs.json")
    
    # If file does not exist, redirect to manual calibration instead of exiting
    if not os.path.exists(calib_file):
        messagebox.showwarning("Archive Not Found", "No calibration file found. Switching to manual adjustment.")
        return apply_cconversion(gas_name, gas_data)

    try:
        with open(calib_file, "r") as f:
            archive = json.load(f)
    except Exception as e:
        messagebox.showerror("File Error", f"Could not read JSON: {e}")
        return apply_cconversion(gas_name, gas_data)

    # If no records exist for this gas, redirect to manual calibration
    if gas_name not in archive or not archive[gas_name].get("history"):
        messagebox.showinfo("No History", f"No calibration records found for {gas_name}. Switching to manual adjustment.")
        return apply_cconversion(gas_name, gas_data)

    # # 2. Extract history and display selection window
    history = archive[gas_name]["history"]
    selection_win = tk.Toplevel()
    selection_win.title(f"Select Calibration for {gas_name}")
    selection_win.geometry("650x450")
    selection_win.grab_set()

    cols = ("Date", "Gain (ai)", "Offset (bi)", "Pressure", "Cross-Gas")
    tree = ttk.Treeview(selection_win, columns=cols, show="headings")
    for col in cols: 
        tree.heading(col, text=col)
        tree.column(col, width=120)
    tree.pack(fill="both", expand=True, padx=10, pady=10)

    for entry in history:
        tree.insert("", "end", values=(
            entry.get("timestamp", "N/A")[:16],
            f"{entry.get('gain_ai', 1.0):.5f}",
            f"{entry.get('offset_bi', 0.0):.4f}",
            f"{entry.get('calib_pressure', 1013):.1f}",
            entry.get("cross_interference", "None")
        ))

    selected_params = {"ai": 1.0, "bi": 0.0, "aij": 0.0}
    status = {"confirmed": False} # To check for final confirmation

    def on_confirm():
        sel = tree.selection()
        if sel:
            idx = tree.index(sel[0])
            record = history[idx]
            selected_params["ai"] = record.get("gain_ai", 1.0)
            selected_params["bi"] = record.get("offset_bi", 0.0)
            selected_params["aij"] = record.get("cross_interference", 0.0)
            status["confirmed"] = True
            selection_win.destroy()

    ttk.Button(selection_win, text="Apply Selected Record", command=on_confirm).pack(pady=10)
    selection_win.wait_window()

    # If the user closed the window without clicking the confirm button
    if not status["confirmed"]:
        return apply_cconversion(gas_name, gas_data)

    # 3. Mathematical calculations (strictly according to the formula in the paper)
    try:
        raw_signal = np.array(gas_data, dtype=float)
        ai = float(selected_params["ai"])
        bi = float(selected_params["bi"])
        aij = float(selected_params["aij"])

        # Core calibration formula
        calibrated_conc = (raw_signal - bi) / ai 
        
        # Cross-interference correction for H2S on SO2 (if total dataset is available)
        if gas_name == "H2S" and aij != 0 and all_data_df is not None:
             if "SO2" in all_data_df.columns:
                 so2_data = all_data_df["SO2"].values
                 # Inverse matrix formula: Ci = (Ri - bi - aij * Cj) / ai
                 calibrated_conc = (raw_signal - bi - (aij * so2_data)) / ai

        calibrated_conc[calibrated_conc < 0] = 0
        
        # Since the data is successfully calibrated using archive constants, return directly.
        return calibrated_conc

    except Exception as e:
        messagebox.showerror("Processing Error", f"Error in calculation: {e}")
        return gas_data
def apply_cconversion(gas_name, gas_data):
    calibrate = messagebox.askyesno("Calibration", f"Do you need to manually calibrate or adjust {gas_name} levels?" )
    if calibrate:
        formula = simpledialog.askstring("Calibration Formula", f"Please enter the formula to apply to {gas_name} (e.g., x * 2 for multiplication, x / 1.5 for division):")
        try:
            gas_data = eval(formula, {"x": gas_data, "np": np})
        except Exception as e:
            messagebox.showerror("Error", f"Invalid formula: {e}")
            return gas_data
    return gas_data

def correct_saturation(gas_name, gas_data, x_data):
    """Corrects gas sensor saturation using real or synthetic slopes based on healthy point intervals."""
    
    if not messagebox.askyesno("Saturation Check", f"Is {gas_name} saturated?"):
        return gas_data

    y_saturation = np.max(gas_data)
    saturated_indices = np.where(gas_data >= y_saturation)[0]

    if not len(saturated_indices):
        return gas_data

    threshold_healthy = 2  # Minimum required healthy points for boundary slope definition

    # Initial boundary markers
    n1 = saturated_indices[0] - 1 if saturated_indices[0] > 0 else 0
    n2 = saturated_indices[-1] + 1 if saturated_indices[-1] < len(gas_data) - 1 else len(gas_data) - 1

    # --- Boundary adjustment n1 (pre-saturation evaluation)
    healthy_count_n1 = 0
    while n1 >= 0:
        if gas_data[n1] < y_saturation:
            healthy_count_n1 += 1
            if healthy_count_n1 >= threshold_healthy:
                break
        else:
            healthy_count_n1 = 0  # Reset counter if saturation threshold is re-encountered
        n1 -= 1
    if n1 < 0:
        n1 = 0

    # --- Boundary adjustment n2 (post-saturation evaluation)
    healthy_count_n2 = 0
    while n2 < len(gas_data):
        if gas_data[n2] < y_saturation:
            healthy_count_n2 += 1
            if healthy_count_n2 >= threshold_healthy:
                break
        else:
            healthy_count_n2 = 0  # Reset counter if saturation threshold is re-encountered
        n2 += 1
    if n2 >= len(gas_data):
        n2 = len(gas_data) - 1

    # Calculation of real or synthetic sensor response slopes
    if n1 > 0:
        df_dx_n1_real = (gas_data[n1] - gas_data[n1-1]) / (x_data[n1] - x_data[n1-1])
    else:
        df_dx_n1_real = None

    if n2 + 1 < len(gas_data):
        df_dx_n2_real = (gas_data[n2+1] - gas_data[n2]) / (x_data[n2+1] - x_data[n2])
    else:
        df_dx_n2_real = None

    df_dx_n1 = df_dx_n1_real if df_dx_n1_real is not None else (df_dx_n2_real * 0.8 if df_dx_n2_real is not None else 0.001)
    df_dx_n2 = df_dx_n2_real if df_dx_n2_real is not None else (df_dx_n1_real * 0.8 if df_dx_n1_real is not None else 0.001)

    # Reconstruct peak coordinates using line intersections
    if (df_dx_n1 - df_dx_n2) != 0:
        x_peak = (df_dx_n1 * x_data[n1] - df_dx_n2 * x_data[n2] + gas_data[n2] - gas_data[n1]) / (df_dx_n1 - df_dx_n2)
        y_peak = df_dx_n1 * (x_peak - x_data[n1]) + gas_data[n1]
    else:
        x_peak = (x_data[n1] + x_data[n2]) / 2
        y_peak = (gas_data[n1] + gas_data[n2]) / 2

    corrected_gas = gas_data.copy()

    # Synthetic correction of saturated data segment
    half = len(saturated_indices) // 2
    first_half = saturated_indices[:half+1]
    second_half = saturated_indices[half+1:]

    for i in first_half:
        corrected_gas[i] = df_dx_n1 * (x_data[i] - x_data[n1]) + gas_data[n1]

    for i in second_half:
        corrected_gas[i] = df_dx_n2 * (x_data[i] - x_data[n2]) + gas_data[n2]

    return corrected_gas
    

# Run analysis
def run_analysis(gas1_name, gas2_name):
    global processed_data, map_button, intercept, selected_range, raw_multigas_data
    selected_range = []
    try:
        data = load_excel_file()
        if data is None:
            return
        raw_multigas_data = data 
        if gas1_name not in data.columns or gas2_name not in data.columns:
            messagebox.showerror("Error", f"Columns {gas1_name} and/or {gas2_name} not found in data.")
            return

        gas1 = np.array(data[gas1_name])
        gas2 = np.array(data[gas2_name])

        x_data = np.array(data.index)

      
        gas1 = apply_calibration_smart(gas1_name, gas1, data)
        gas2 = apply_calibration_smart(gas2_name, gas2, data)

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(x_data, gas1, label=f'{gas1_name} (Raw)', color="red")
        plt.xlabel('Index')
        plt.ylabel(f'{gas1_name}')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.plot(x_data, gas2, label=f'{gas2_name} (Raw)', color="orange")
        plt.xlabel('Index')
        plt.ylabel(f'{gas2_name}')
        plt.legend()
        plt.show()

        # Saturation evaluation and reconstruction for Gas 1
        gas1 = correct_saturation(gas1_name, gas1, x_data)

        # Saturation evaluation and reconstruction for Gas 2
        gas2 = correct_saturation(gas2_name, gas2, x_data)

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(data.index, gas1, label=f'{gas1_name}')
        plt.xlabel('Index')
        plt.ylabel(f'{gas1_name}')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.plot(data.index, gas2, label=f'{gas2_name}')
        plt.xlabel('Index')
        plt.ylabel(f'{gas2_name}')
        plt.legend()
        plt.show()

        # Render evaluation plots
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        axs[0].plot(data.index, gas1, label=f'{gas1_name}')
        axs[0].set_xlabel('Index')
        axs[0].set_ylabel(f'{gas1_name}')
        axs[0].legend()

        axs[1].plot(data.index, gas2, label=f'{gas2_name}')
        axs[1].set_xlabel('Index')
        axs[1].set_ylabel(f'{gas2_name}')
        axs[1].legend()

        # Interactive span selector assigned to the primary gas timeline
        span = SpanSelector(axs[0], onselect, 'horizontal', useblit=True,
                    props=dict(alpha=0.5, facecolor='red'))


        plt.suptitle("Select the range on the first plot (gas1), or close the window to enter manually")       
        plt.show(block=True)  

        if not done or not selected_range:
           messagebox.showinfo("No Selection", "No range was selected. Please enter it manually.")
           start_row, end_row = get_initial_range()
        else:
            # If sub-range is interactively chosen, prioritize selected_range array boundaries
            start_row, end_row = selected_range[0], selected_range[1]           

        # Slice time-series arrays to targeted analytical interval
        gas1 = gas1[start_row:end_row]
        gas2 = gas2[start_row:end_row]
      

        start_time, end_time = 0, 30

        b = []
        for i in range(start_time, end_time + 1):
            if i >= len(gas2):
                break
            gas2_aligned = gas2[i:]
            gas1_aligned = gas1[:len(gas2_aligned)]
            correlation = np.corrcoef(gas1_aligned, gas2_aligned)[0, 1]
            b.append((i, correlation))

        best_alignment = max(b, key=lambda x: x[1])
        gas2_aligned = gas2[best_alignment[0]:]
        gas1_aligned = gas1[:len(gas2_aligned)]

        messagebox.showinfo("Best Correlation", f"Best correlation found at {best_alignment[0]} seconds with correlation {best_alignment[1]:.2f}")

        # Detrending prompt for the primary gas
        detrend_x = messagebox.askyesno("Detrending", f"Do you want to perform detrending for {gas1_name}?")
        if detrend_x:
            x = signal.detrend(gas1_aligned) + np.mean(gas1_aligned)
        else:
            x = gas1_aligned

        # Detrending prompt for the secondary gas
        detrend_y = messagebox.askyesno("Detrending", f"Do you want to perform detrending for {gas2_name}?")
        if detrend_y:
            y = signal.detrend(gas2_aligned) + np.mean(gas2_aligned)
        else:
            y = gas2_aligned

        # Display comparative plots (only if at least one gas sequence was detrended)
        if detrend_x or detrend_y:
            plt.figure(figsize=(12, 6))
            
            # Primary gas time-series subplot
            plt.subplot(1, 2, 1)
            plt.plot(gas1_aligned, label=f'Original {gas1_name}', alpha=0.5)
            plt.plot(x, label=f'Processed {gas1_name}', color='blue')
            plt.title(f"Detrending: {gas1_name}")
            plt.legend()

            # Secondary gas time-series subplot
            plt.subplot(1, 2, 2)
            plt.plot(gas2_aligned, label=f'Original {gas2_name}', alpha=0.5)
            plt.plot(y, label=f'Processed {gas2_name}', color='green')
            plt.title(f"Detrending: {gas2_name}")
            plt.legend()
            
            plt.tight_layout()
            plt.show()

        # 2. Apply exponential smoothing using newly processed x and y datasets
        smoother = GasSmoother()
        x_smoothed, start, end = smoother.apply_selected_filter(x, compare_gas=y, selected_range=selected_range, gas_name=gas1_name)
        y_smoothed, start, end = smoother.apply_selected_filter(y, compare_gas=x, selected_range=selected_range, gas_name=gas2_name)
        
        # --- Step 2: Display Smoothing Results ---
        plt.figure(figsize=(14, 10))

        # Gas 1 comparison plot (Original vs Smoothed)
        plt.subplot(2, 2, 1)
        plt.plot(x, label=f'Original {gas1_name}', alpha=0.5)
        plt.plot(x_smoothed, label=f'Smoothed {gas1_name}', linewidth=1.5)
        plt.title(f"{gas1_name} Signal Smoothing")
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Gas 2 comparison plot (Original vs Smoothed)
        plt.subplot(2, 2, 2)
        plt.plot(y, label=f'Original {gas2_name}', alpha=0.5)
        plt.plot(y_smoothed, label=f'Smoothed {gas2_name}', linewidth=1.5)
        plt.title(f"{gas2_name} Signal Smoothing")
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Scatter plot of finalized smoothed observations
        plt.subplot(2, 1, 2)
        plt.scatter(y_smoothed, x_smoothed, color='green', alpha=0.5, s=10)
        plt.xlabel(f"{gas2_name} (Smoothed)")
        plt.ylabel(f"{gas1_name} (Smoothed)")
        plt.title(f"Correlation: {gas1_name} vs {gas2_name}")
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()
       # --- Step 3: Cross-Interference Correction via Scientific Formula ---
        
       # 1. Define correction function (accepts interference factor as input parameter)
        def correct_h2s(h2s_meas, so2_meas, factor):
            """Corrects SO₂ cross-interference on H₂S data using standard matrix formula."""
            return np.maximum(h2s_meas - (factor * so2_meas), 0)

        # 2. Build local dataframe environment
        df = pd.DataFrame({'x': x_smoothed, 'y': y_smoothed})
        
        # 3. Validate dual presence condition for SO2 and H2S variables
        if {gas1_name, gas2_name} == {"SO2", "H2S"}:
            
            # Prompt user to initiate cross-interference sequence
            if messagebox.askyesno("Cross-Interference", "Do you want to correct SO2 interference on H2S?"):
                
                # Dynamic column mapping to identify H2S and SO2 variables (x or y)
                h2s_col = "x" if gas1_name == "H2S" else "y"
                so2_col = "y" if gas1_name == "H2S" else "x"
                
                # Fetch F-factor input parameter from user
                from tkinter import simpledialog
                c_factor = simpledialog.askfloat("Correction Factor", 
                                                 "Enter F (Interference factor, e.g. 0.1-0.2):", 
                                                 initialvalue=0.15, minvalue=0.0, maxvalue=1.0)
                
                if c_factor is not None:
                    # Execute cross-interference evaluation function
                    corrected_values = correct_h2s(df[h2s_col], df[so2_col], c_factor)
                    
                    # Store in dedicated column for evaluation plotting
                    df["h2s_corrected"] = corrected_values

                    # **Plot Evaluation Layout**
                    plt.figure(figsize=(12, 8))
                    plt.plot(df.index, df[h2s_col], label="H₂S (Measured)", linestyle="--", color="red", alpha=0.7)
                    plt.plot(df.index, df["h2s_corrected"], label=f"H₂S (Corrected, F={c_factor})", linestyle="-", color="green", linewidth=2)
                    plt.plot(df.index, df[so2_col], label="SO₂ (Interfering Gas)", linestyle="-.", color="blue", alpha=0.5)
                    
                    plt.xlabel("Index")
                    plt.ylabel("Gas Concentration (PPM)")
                    plt.title(f"Effect of SO₂ on H₂S and Correction (F={c_factor})")
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    plt.show()

                    # Map corrected datasets back to original workspaces for subsequent regression
                    if gas1_name == "H2S":
                        df["x"] = df["h2s_corrected"]
                    else:
                        df["y"] = df["h2s_corrected"]
                        
     # --- Step 4: Regression Direction Configuration ---
     # Prompt User for Ratio Direction Setup ---
        ratio_choice = messagebox.askyesno("Regression Direction", 
            f"Select the Ratio Direction:\n\n"
            f"YES: {gas1_name} / {gas2_name} (Axis Y: {gas1_name})\n"
            f"NO: {gas2_name} / {gas1_name} (Axis Y: {gas2_name})")

        # --- Map DataFrame Vectors to Regression Model Variables ---
        # Note: df["x"] has been seamlessly overwritten with corrected values if steps were executed
        if ratio_choice: #Selection route matching Gas1/Gas2 ratio direction criteria
            dep_var = df['x'].values      # Extract finalized target dependent (Y) values from local storage arrays
            indep_var = df['y'].values    # Extract smoothed independent (X) source metrics from local storage arrays
            y_label, x_label = gas1_name, gas2_name
        else: # Direct Gas Ratio Configuration (gas2 / gas1)
            # Define the dependent variable (Y) for the regression model
            dep_var = df['y'].values
            # Define the independent variable (X) for the regression model 
            indep_var = df['x'].values    
            # Assign dynamic axis labels representing the target gas pairs
            y_label, x_label = gas2_name, gas1_name

        # Build localized analytical dataframe filtering null fields automatically
        reg_df = pd.DataFrame({'X': indep_var, 'Y': dep_var}).dropna()
        msk = np.random.rand(len(reg_df)) < 0.8
        train = reg_df[msk]
        test = reg_df[~msk]

        # Segment remaining 20% of records into unseen validation testing partition
        regr = linear_model.LinearRegression()
        train_x = np.asanyarray(train[['X']])
        train_y = np.asanyarray(train[['Y']])
        regr.fit(train_x, train_y)

        # Calculate dual-partition model predictions for comparative variance analysis 
        train_y_pred = regr.predict(train_x)
        test_x = np.asanyarray(test[['X']])
        test_y = np.asanyarray(test[['Y']])
        test_y_pred = regr.predict(test_x)

        # Extract geostatistical descriptors and model parameters
        coef = regr.coef_[0][0]
        intercept = regr.intercept_[0]
        # Evaluate Partition Performance and Model Generalization Gap
        r2_train = r2_score(train_y, train_y_pred) # Evaluate training subset goodness-of-fit
        r2_test = r2_score(test_y, test_y_pred)   # Evaluate validation testing subset performance
        r2_diff = abs(r2_train - r2_test)         # Quantify generalization variance gap to detect overfitting
        #Compute Residual Error Metrics and Linear Correlation
        mae = np.mean(np.absolute(test_y_pred - test_y))  # Calculate Mean Absolute Error on validation data
        mse = np.mean((test_y_pred - test_y) ** 2)        # Calculate Mean Squared Error for variance penalization
        corr = np.corrcoef(reg_df['X'], reg_df['Y'])[0, 1]# Extract Pearson correlation coefficient between processed gas pairs 

        # Render Cross-Validation Diagnostic Scatter Plot
        plt.figure(figsize=(9, 7))
        
        # Plot training data subset in low opacity blue
        plt.scatter(train.X, train.Y, color='blue', alpha=0.3, label=f'Train Data ($R^2={r2_train:.3f}$)')
        
        # Plot testing data subset in high contrast green for unseen sample verification
        plt.scatter(test.X, test.Y, color='green', alpha=0.7, edgecolors='k', label=f'Test Data ($R^2={r2_test:.3f}$)')
        
        # Plot the primary linear regression trendline in distinct solid red
        plt.plot(train_x, coef * train_x + intercept, '-r', linewidth=2.5, label='Regression Model')
        
        # Construct dynamic summary text box showing regression parameters and generalization gap
        stats_label = (f"Coefficient (Slope): {coef:.4f}\n"
                       f"Train $R^2$: {r2_train:.3f}\n"
                       f"Test $R^2$: {r2_test:.3f}\n"
                       f"Diff: {r2_diff:.3f}")
        
        plt.gca().text(0.05, 0.95, stats_label, transform=plt.gca().transAxes, 
                       fontsize=11, verticalalignment='top', 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

        plt.xlabel(f"{x_label} (Processed)")
        plt.ylabel(f"{y_label} (Processed)")
        plt.title(f"VGR Analysis: {y_label} / {x_label}")
        plt.legend(loc='lower right')
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.show()

       # Compile Statistical Log Report with Stability Warning Indicators ---
        status = "STABLE" if r2_diff < 0.1 else "POTENTIAL HETEROGENEITY (Check for secondary sources)"
        
        results_message = (
            f"--- VGR Processor Report ---\n"
            f"Target Ratio: {y_label} / {x_label}\n"
            f"---------------------------\n"
            f"Coefficient (Slope): {coef:.4f}\n"
            f"Train R2: {r2_train:.4f}\n"
            f"Test R2: {r2_test:.4f}\n"
            f"R2 Difference: {r2_diff:.4f}\n"
            f"---------------------------\n"
            f"Intercept: {intercept:.4f}\n"
            f"Correlation Coefficient (r): {corr:.4f}\n"
            f"Mean Absolute Error: {mae:.4f}\n"
            f"Mean Squared Error: {mse:.4f}\n"
        )

        messagebox.showinfo("Analysis Results", results_message)
        # Smart Storage Protocol Using Dynamically Resolved Source Prefixes
        if messagebox.askyesno("Save Results", "Would you like to save these results to a text file?"):
            try:
                import os
                # Extract file base name from root global path reference built during initialization
                base_name = os.path.basename(file_path)
                file_name_no_ext = os.path.splitext(base_name)[0]
                # Format automatic suggested export file identifier using base name and gas pairs
                suggested_name = f"Results_{file_name_no_ext}_{y_label}_{x_label}.txt"
            except NameError:
                # Fallback template name if base file path pointer is unresolved
                suggested_name = f"VGR_Processor_{y_label}_{x_label}.txt"

            save_file_path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                initialfile=suggested_name,
                filetypes=[("Text files", "*.txt")],
                title="Save Your Analysis"
            )
            
            if save_file_path:
                with open(save_file_path, "w", encoding="utf-8") as f:
                    f.write(results_message)
                messagebox.showinfo("Saved", f"Results saved successfully as:\n{os.path.basename(save_file_path)}")

        # Cache Active Analytical Slice to Global Variables for Cartography
        global processed_data, processed_data_start_index, processed_data_end_index
        processed_data = pd.DataFrame({'x': indep_var, 'y': dep_var}) 
        processed_data_start_index = start
        processed_data_end_index = end
        
        # Evaluate total record count within currently processed timestamp range
        num_processed = end - start + 1
            
        messagebox.showinfo(
    "Generate Map Instructions",
    f"{num_processed} rows were processed.\nTo generate the map, load all coordinates.\nThe software will automatically select the range from row {start} to {end}."
)
        
        # Enable map generation trigger widget upon validation
        map_button.config(state=tk.NORMAL)
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred during analysis: {e}")   
        
def save_shapefile_pyshp(data, gas_ratio, output_path):
    """Save processed data and gas ratio as a shapefile."""
    with shapefile.Writer(output_path) as shp:
        # Define attribute schema field names and floating-point data types
        shp.field('Gas_Ratio', 'F')  
        
        for index, row in data.iterrows():
            # Populate vector geometry by injecting longitude and latitude coordinate pairs
            shp.point(row['Longitude'], row['Latitude'])
            # Commit the corresponding gas ratio attribute to the Shapefile database record
            shp.record(gas_ratio[index])  

    print(f"Shapefile saved to: {output_path}")

def generate_map():
    try:
        # --- Geostatistical Validation and Data Source Checks ---
        if 'processed_data' not in globals():
            messagebox.showerror("Error", "No processed data available.")
            return
        if 'processed_data_start_index' not in globals() or 'processed_data_end_index' not in globals():
            messagebox.showerror("Error", "Processed data range not defined.")
            return
        if 'raw_multigas_data' not in globals():
            messagebox.showerror("Error", "Raw multigas data not found.")
            return

        start_index = processed_data_start_index
        end_index = processed_data_end_index

        from xml.etree import ElementTree as ET
        from scipy.interpolate import interp1d, PchipInterpolator

        # --- Launch Location Tracking Log File Dialog ---
        coordinates_file = filedialog.askopenfilename(
            filetypes=[("KML files", "*.kml"), ("Excel files", "*.xlsx")]
        )
        if not coordinates_file:
            return

        # --- Core Parser Routine to Translate KML Tracks into DataFrames ---
        def read_kml_to_df(path):
            ns = {'kml': 'http://www.opengis.net/kml/2.2', 'gx': 'http://www.google.com/kml/ext/2.2'}
            tree = ET.parse(path)
            root = tree.getroot()

            lats, lons, whens = [], [], []

            for track in root.findall('.//gx:Track', ns):
                coord_list = track.findall('gx:coord', ns)
                time_list  = track.findall('kml:when', ns)
                n = min(len(coord_list), len(time_list))
                for i in range(n):
                    parts = coord_list[i].text.strip().split()
                    if len(parts) >= 2:
                        lon, lat = float(parts[0]), float(parts[1])
                        lons.append(lon)
                        lats.append(lat)
                        whens.append(time_list[i].text.strip())

            for pm in root.findall('.//kml:Placemark', ns):
                ts = pm.find('.//kml:TimeStamp/kml:when', ns)
                pt = pm.find('.//kml:Point/kml:coordinates', ns)
                if ts is not None and pt is not None and pt.text:
                    coords_txt = pt.text.strip()
                    try:
                        lonlatalt = coords_txt.split(',')
                        lon = float(lonlatalt[0])
                        lat = float(lonlatalt[1])
                        lons.append(lon)
                        lats.append(lat)
                        whens.append(ts.text.strip())
                    except Exception:
                        pass

            if not whens or not lats or not lons:
                return None

            df = pd.DataFrame({
                'Latitude': lats,
                'Longitude': lons,
                'Time': pd.to_datetime(whens, errors='coerce', utc=True)
            })
            df.dropna(subset=['Time'], inplace=True)
            if df.empty:
                return None

            df.sort_values('Time', inplace=True)
            df['Time'] = df['Time'].dt.tz_convert(None)
            return df

        # --- Spatial Data Ingestion and Format-Specific Validation ---
        if coordinates_file.lower().endswith('.kml'):
            coordinates_data = read_kml_to_df(coordinates_file)
            if coordinates_data is None or coordinates_data.empty:
                messagebox.showerror("Error", "No valid GPS data found in KML file.")
                return
        else:  # Excel
            coordinates_data = pd.read_excel(coordinates_file)
            for col in ['Latitude', 'Longitude', 'Time']:
                if col not in coordinates_data.columns:
                    messagebox.showerror("Error", f"Excel file must contain '{col}' column.")
                    return
            coordinates_data['Time'] = pd.to_datetime(coordinates_data['Time'], errors='coerce', utc=True)
            coordinates_data.dropna(subset=['Time'], inplace=True)
            if coordinates_data.empty:
                messagebox.showerror("Error", "Excel file contains no valid times.")
                return
            coordinates_data.sort_values('Time', inplace=True)
            coordinates_data['Time'] = coordinates_data['Time'].dt.tz_convert(None)

        # --- Process Selected Vectors via Extension Validation ---
        gps_seconds = (coordinates_data['Time'] - coordinates_data['Time'].iloc[0]).dt.total_seconds().values
        gas_seconds = (pd.to_datetime(raw_multigas_data['Time']) - pd.to_datetime(raw_multigas_data['Time']).iloc[0]).dt.total_seconds().values

        # --- Apply Multi-Rate 1D Interpolation for Sparse GPS Sampling Density ---
        total_multigas_length = len(raw_multigas_data)
        if len(coordinates_data) < total_multigas_length:
            choice = simpledialog.askstring(
                "Interpolation Required",
                f"Coordinates file is too short ({len(coordinates_data)} < {total_multigas_length}).\n\n"
                "Select interpolation method:\n"
                "1. Linear  → Suitable for scattered data / data with time gaps\n"
                "2. PCHIP   → Suitable for continuous data (without sharp fluctuations)"
            )

            if choice not in ["1", "2"]:
                messagebox.showerror("Error", "Invalid choice.")
                return

            lat_old = coordinates_data['Latitude'].values
            lon_old = coordinates_data['Longitude'].values

            if choice == "1":
                f_lat = interp1d(gps_seconds, lat_old, kind='linear', fill_value="extrapolate")
                f_lon = interp1d(gps_seconds, lon_old, kind='linear', fill_value="extrapolate")
            else:
                f_lat = PchipInterpolator(gps_seconds, lat_old, extrapolate=True)
                f_lon = PchipInterpolator(gps_seconds, lon_old, extrapolate=True)

            lat_new = f_lat(gas_seconds)
            lon_new = f_lon(gas_seconds)
            coordinates_data = pd.DataFrame({
                'Latitude': lat_new,
                'Longitude': lon_new,
                'Time': pd.to_datetime(raw_multigas_data['Time'])
            })

        # --- Align Geospatial Tracking Coordinate Arrays with Target Metrics Bound ---
        expected_length = len(processed_data)
        if len(coordinates_data) < expected_length:
            messagebox.showerror("Error", f"Coordinates file too short. Expected {expected_length}, got {len(coordinates_data)}.")
            return

        processed_data['Latitude'] = coordinates_data['Latitude'].iloc[:expected_length].values
        processed_data['Longitude'] = coordinates_data['Longitude'].iloc[:expected_length].values
        data = processed_data.dropna(subset=['x', 'y', 'Latitude', 'Longitude']).copy()

    except Exception as e:
        messagebox.showerror("Error", f"Error processing coordinates: {e}")
        return

    # --- Spatial Kriging Interpolation and Cartography Production Module ---
    try:
        ratio_choice = simpledialog.askstring("Gas Ratio", "Select gas ratio to plot:\n1. Original (x/y)\n2. Inverted (y/x)")
        ratio_label = "Original (x/y)"
        use_inverse = False
        if ratio_choice=="2":
            ratio_label="Inverted (y/x)"
            use_inverse=True

        total_points = len(data)
        group_size = max(30, min(300, total_points//50))
        num_groups = total_points // group_size

        slope_points = []
        for i in range(num_groups):
            group = data.iloc[i*group_size:(i+1)*group_size]
            if len(group) >= 5:
                model = LinearRegression()
                X = group[['y']].values
                y_local = group['x'].values
                model.fit(X, y_local)
                slope = model.coef_[0]
                inv_slope = 1/slope if slope!=0 else np.nan
                slope_points.append({
                    'lat': group['Latitude'].mean(),
                    'lon': group['Longitude'].mean(),
                    'slope': slope,
                    'inv_slope': inv_slope
                })

        if not slope_points:
            messagebox.showinfo("Info","No valid groups for regression.")
            return

        lons = np.array([pt['lon'] for pt in slope_points])
        lats = np.array([pt['lat'] for pt in slope_points])
        slopes = np.array([pt['inv_slope'] if use_inverse else pt['slope'] for pt in slope_points])

        OK = OrdinaryKriging(lons, lats, slopes, variogram_model='linear', verbose=False, enable_plotting=False)
        grid_lon = np.linspace(min(lons), max(lons), 100)
        grid_lat = np.linspace(min(lats), max(lats), 100)
        gridx, gridy = np.meshgrid(grid_lon, grid_lat)
        z, _ = OK.execute("grid", grid_lon, grid_lat)

        plt.figure(figsize=(10,8))
        levels = np.linspace(np.nanmin(z), np.nanmax(z), 10)
        contourf = plt.contourf(gridx, gridy, z, levels=levels, cmap='Greys')
        plt.contour(gridx, gridy, z, levels=levels, colors='black', linewidths=0.3)
        cbar = plt.colorbar(contourf)
        cbar.set_label(f'Gas Ratio (Slope: {ratio_label})')
        plt.xlabel("Longitude"); plt.ylabel("Latitude")

        map_title = simpledialog.askstring("Map Title", "Enter the title of the map:")
        if map_title: plt.title(map_title, fontsize=14)

        map_output_path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF files","*.pdf")])
        if map_output_path: plt.savefig(map_output_path,dpi=300)
        plt.close()

        slope_shp_path = filedialog.asksaveasfilename(defaultextension=".shp", title="Save slope points shapefile", filetypes=[("Shapefiles","*.shp")])
        if slope_shp_path:
            with shapefile.Writer(slope_shp_path, autoBalance=1) as shp:
                shp.field('Slope','F',decimal=6)
                for pt in slope_points:
                    value = pt['inv_slope'] if use_inverse else pt['slope']
                    shp.point(pt['lon'], pt['lat'])
                    shp.record(Slope=round(float(value),6))

        print("Map and shapefile saved successfully.")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while generating the map: {e}")

# --- Workspace Windows and Asset Management Helpers ---

def exit_program():
    root.quit()


def resource_path(relative_path):
    """ Resolve local asset storage pathways for PyInstaller bundles and IDE environments. """
    try:
        # PyInstaller generates dynamic runtime folders and stores the path in the _MEIPASS property
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def open_manual():
    # Call runtime path locator to retrieve valid local user manual PDF location
    file_path = resource_path("User Manual.pdf") 
    
    if os.path.exists(file_path):
        # Open documentation using default host application environment hooks
        webbrowser.open(f'file:///{file_path}')
    else:
        messagebox.showwarning("File Not Found", f"The manual was not found at:\n{file_path}")

def open_calibration_window():
    """ Secondary Window Handler for System Calibration Tools """
    calib_win = tk.Toplevel(root)
    # Instantiate calibration toolbox context matching core layout declarations
    app = VGR_CalibrationToolbox(calib_win)

# --- Main Tkinter Application Window Initialization ---

root = tk.Tk()
root.title("VGR - Volcanic Gas Ratio Processor")
root.geometry("900x700")

# ---  Application Visual Branding and Graphic Assets ---
# Load and configure the main background imagery and institutional logo for the window environment
script_dir = os.path.dirname(os.path.abspath(__file__))
bg_image_path = os.path.join(script_dir, "backgr2.png")
logo_path = os.path.join(script_dir, "logo.png")

try:
    original_bg = Image.open(bg_image_path)
    bg_image_tk = ImageTk.PhotoImage(original_bg)
    background_label = tk.Label(root, image=bg_image_tk)
    background_label.place(x=0, y=0, relwidth=1, relheight=1)

    logo_image = Image.open(logo_path).resize((180, 180), Image.LANCZOS)
    logo_image_tk = ImageTk.PhotoImage(logo_image)
    logo_label = tk.Label(root, image=logo_image_tk, bg="#F0F0F0")
    logo_label.pack(pady=10)
except Exception as e:
    print(f"Image Error: {e}")

# ---  Application Greeting Interface Component ---
# Instantiate and render the primary header banner to display the platform's title at the top of the interface
welcome_label = tk.Label(root, text="Volcanic Gas Ratio Processor (VGR)", 
                         font=("Helvetica", 18, "bold"), bg="#F0F0F0", fg="#333333")
welcome_label.pack(pady=10)

# --- Structural Container Layout for Control Grid ---
# Create a dedicated layout frame to systematically organize operational widgets using grid coordinates
buttons_frame = tk.Frame(root, bg="#F0F0F0")
buttons_frame.pack(pady=10)

# ---  Admin Controls Left Flank Layout ---
# Instantiate and position calibration tool trigger on the left flank of the initial grid row
calibration_button = tk.Button(
    buttons_frame, text="Calibration Tool", command=open_calibration_window, 
    bg='sky blue', font=("Helvetica", 11, "bold"), width=18
)
calibration_button.grid(row=0, column=0, padx=10, pady=10, sticky="w")

# --- Top-Level Admin Controls Layout ---
# Instantiate and position user guide reference launcher on the right flank of the initial grid row
manual_button = tk.Button(
    buttons_frame, text="User Manual", command=open_manual, 
    bg='light green', font=("Helvetica", 11, "bold"), width=18
)
manual_button.grid(row=0, column=1, padx=10, pady=10, sticky="e")

# ---  Mid-Level UI Grid Layout for Gas Ratio Computations ---
# Define structural collection mapping dynamic gas analysis titles to their respective logic routines
analysis_list = [
    ("Analyze H2S / SO2", analyze_h2s_so2),
    ("Analyze CO2 / H2S", analyze_co2_h2s),
    ("Analyze H2O / CO2", analyze_h2o_co2),
    ("Analyze H2O / H2S", analyze_h2o_h2s),
    ("Analyze H2O / SO2", analyze_h2o_so2),
    ("Analyze CO2 / SO2", analyze_co2_so2)
]

for i, (txt, cmd) in enumerate(analysis_list):
    r = (i // 2) + 1
    c = i % 2
    tk.Button(buttons_frame, text=txt, command=cmd, bg='silver', 
              font=("Helvetica", 10), width=22).grid(row=r, column=c, padx=5, pady=5, sticky="ew")

# ---  Final UI Grid Layout Configuration ---
# Position spatial cartography action trigger in the bottom row of the button container
map_button = tk.Button(buttons_frame, text="Generate Map", command=generate_map, 
                       bg='gold', font=("Helvetica", 12, "bold"), state="disabled")
map_button.grid(row=4, column=0, columnspan=2, padx=5, pady=15, sticky="ew")

# ---  Exit Interface and System Shutdown ---
# Position exit command trigger at the base of the user interface
exit_button = tk.Button(root, text="Exit Program", font=("Helvetica", 12, "bold"), 
                        bg="#FF6347", fg="white", padx=20, command=exit_program)
exit_button.pack(side="bottom", pady=20)


root.mainloop()
