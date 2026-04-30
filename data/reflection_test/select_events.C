#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

/**
 * Event selection based on charge threshold.
 * @param input_file Path to the input ROOT file.
 * @param threshold_ch0 Threshold for channel 0 in pC.
 * @param threshold_ch1 Threshold for channel 1 in pC.
 * @param output_file Path to the output ROOT file for selected events.
 */
void select_events(const char* input_file = "output.root", 
                   double threshold_ch0 = 10.0, 
                   double threshold_ch1 = 10.0,
                   const char* output_file = "selected_events.root") {
    
    // Use RDataFrame for efficient data processing
    ROOT::RDataFrame df("T_Charge", input_file);

    // Filter events where both channels exceed the threshold
    // Note: charge_pC is a double[2], so we access indices 0 and 1.
    auto df_filtered = df.Filter([threshold_ch0, threshold_ch1](const ROOT::VecOps::RVec<double>& charge) {
        return charge[0] > threshold_ch0 && charge[1] > threshold_ch1;
    }, {"charge_pC"});

    // Count the selected events
    auto count = df_filtered.Count();
    std::cout << "Selected events: " << *count << std::endl;

    // Save the filtered tree to a new file
    df_filtered.Snapshot("T_Charge_Selected", output_file);
    
    std::cout << "Saved selected events to " << output_file << std::endl;
}
