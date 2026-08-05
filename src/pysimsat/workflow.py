from .genRSP import gen_rsp
from .genBKG import gen_background
from .genOBS import gen_observation

def run_all(instrument, source, instrument_json, source_file, output_dir, spec_dir, resp_dir):
    """
    Run the complete PySimSat simulation pipeline.

    Steps:
    1. Generate response files
    2. Generate background spectrum
    3. Generate observation spectrum
    """
    # Generate response files
    print("\nGenerating response files...\n")
    chars, sourcechars, background, mission = gen_rsp(instrument = instrument, source = source, instrument_json = instrument_json, source_file = source_file, output_dir = output_dir, resp_dir = resp_dir)

    # Generate background spectrum
    print("\nGenerating background spectrum...\n")
    gen_background(background, mission, sourcechars["exposure"], chars["num_det_pixels"], output_dir, spec_dir, resp_dir)

    # Generate observation spectrum
    print("\nGenerating observation spectrum...\n")
    gen_observation(mission, sourcechars, output_dir, spec_dir, resp_dir)