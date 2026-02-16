def colors(columns_names):
    cell_colors_map = {
        "Cortical_neurons": "#1f77b4",      # blue
        "Cortical_neurons_microglia": "#ff7f0e",  # orange
        "organoids": "#d62728",              # red
        "Microglia": "#2ca02c"             # green
    }    
    treatment_colors_map = {
        "mock":     "#FFC0CB",     # pink
        "SARS-CoV": "#000000"               # black (fill)
    }
       # Step 6: Define color and marker mappings
    treatment_markers = {
        "mock": "o",        # circle
        "SARS-CoV": "s"     # square
    }
    # Extract metadata from column names
    # Parse cell type and treatment from sample names
    cell_types = []
    treatments = []
    for col in columns_names:
        if "Cortical_neurons_microglia" in col:
            cell_type = "Cortical_neurons_microglia"
        elif "Cortical_neurons" in col:
            cell_type = "Cortical_neurons"
        elif "Microglia" in col:
            cell_type = "Microglia"
        elif "organoids" in col:
            cell_type = "organoids"
        else:
            cell_type = "unknown"     
        if "mock" in col:
            treatment = "mock"
        elif "SARS-CoV" in col or "SARS_CoV" in col:
            treatment = "SARS-CoV"
        else:
            treatment = "unknown"       
        cell_types.append(cell_type)
        treatments.append(treatment)
    return cell_colors_map, treatment_colors_map, treatment_markers, cell_types, treatments

def sort_table(df):
    # Define custom order for cell type and treatment combinations
    cell_treatment_order = [
        "Cortical_neurons_mock",
        "Cortical_neurons_SARS-CoV",
        "organoids_mock",
        "organoids_SARS-CoV",
        "Cortical_neurons_microglia_mock",
        "Cortical_neurons_microglia_SARS-CoV",
        "Microglia_mock",
        "Microglia_SARS-CoV"]
    sorted_cols = []
    for sample in cell_treatment_order:
        matching_cols = [col for col in df.columns if sample in col]
        sorted_cols.extend(matching_cols)
    df = df[sorted_cols]
    return df       

if __name__ == "__main__":
    main() 