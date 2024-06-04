import xml.etree.ElementTree as ET
import sys
import os
from tqdm import tqdm

def calculate_outgoing_protons(photon_data, initial_proton_momentum):
    total_initial_energy, total_initial_px, total_initial_py, total_initial_pz = initial_proton_momentum

    photon1_data, photon2_data = photon_data  # Unpack photon data for photon 1 and photon 2

    # Initialize variables for photon 1
    photon1_px, photon1_py, photon1_pz, photon1_energy = 0, 0, 0, 0

    # Extract photon 1 data
    _, _, _, _, _, _, photon1_px, photon1_py, photon1_pz, photon1_energy, _, _, _ = photon1_data.split()
    photon1_px = float(photon1_px)
    photon1_py = float(photon1_py)
    photon1_pz = float(photon1_pz)
    photon1_energy = float(photon1_energy)

    # Calculate momentum of outgoing proton 1 based on photon 1 data
    outgoing_proton1_energy = total_initial_energy - photon1_energy
    outgoing_proton1_px = -photon1_px
    outgoing_proton1_py = -photon1_py
    if photon1_pz > 0:
        outgoing_proton1_pz = total_initial_pz - photon1_pz
    else:
        outgoing_proton1_pz = - total_initial_pz + abs(photon1_pz)

    # Initialize variables for photon 2
    photon2_px, photon2_py, photon2_pz, photon2_energy = 0, 0, 0, 0

    # Extract photon 2 data
    _, _, _, _, _, _, photon2_px, photon2_py, photon2_pz, photon2_energy, _, _, _ = photon2_data.split()
    photon2_px = float(photon2_px)
    photon2_py = float(photon2_py)
    photon2_pz = float(photon2_pz)
    photon2_energy = float(photon2_energy)

    # Calculate momentum of outgoing proton 2 based on photon 2 data
    outgoing_proton2_energy = total_initial_energy - photon2_energy
    outgoing_proton2_px = -photon2_px
    outgoing_proton2_py = -photon2_py
    if photon2_pz > 0:
        outgoing_proton2_pz = total_initial_pz - photon2_pz
    else:
        outgoing_proton2_pz = - total_initial_pz + abs(photon2_pz)

    outgoing_proton1_momentum = (outgoing_proton1_energy, outgoing_proton1_px, outgoing_proton1_py, outgoing_proton1_pz)
    outgoing_proton2_momentum = (outgoing_proton2_energy, outgoing_proton2_px, outgoing_proton2_py, outgoing_proton2_pz)

#    print("Incoming Proton 4-Momentum:")
#    print("Energy:", total_initial_energy)
#    print("px:", total_initial_px)
#    print("py:", total_initial_py)
#    print("pz:", total_initial_pz)
#
#    print("\nOutgoing Proton 1 4-Momentum:")
#    print("Energy:", outgoing_proton1_energy)
#    print("px:", outgoing_proton1_px)
#    print("py:", outgoing_proton1_py)
#    print("pz:", outgoing_proton1_pz)
#
#    print("\nOutgoing Proton 2 4-Momentum:")
#    print("Energy:", outgoing_proton2_energy)
#    print("px:", outgoing_proton2_px)
#    print("py:", outgoing_proton2_py)
#    print("pz:", outgoing_proton2_pz)
#
#    print("\nPhoton 1 4-Momentum:")
#    print("px:", photon1_px)
#    print("py:", photon1_py)
#    print("pz:", photon1_pz)
#    print("Energy:", photon1_energy)
#
#    print("\nPhoton 2 4-Momentum:")
#    print("px:", photon2_px)
#    print("py:", photon2_py)
#    print("pz:", photon2_pz)
#    print("Energy:", photon2_energy)

    return outgoing_proton1_momentum, outgoing_proton2_momentum


def add_outgoing_protons_to_event(event, proton1_momentum, proton2_momentum):
    proton1_data = (
        f"       2212 1 1 2 0 0 "
        f"{proton1_momentum[1]:+14.10e} {proton1_momentum[2]:+14.10e} {proton1_momentum[3]:+14.10e} "
        f"{proton1_momentum[0]:+14.10e} 938.2720813 0.0000e+00 1.0000e+00"
    )
    proton2_data = (
        f"       2212 1 1 2 0 0 "
        f"{proton2_momentum[1]:+14.10e} {proton2_momentum[2]:+14.10e} {proton2_momentum[3]:+14.10e} "
        f"{proton2_momentum[0]:+14.10e} 938.2720813 0.0000e+00 1.0000e+00\n"  # Add line break
    )

    event_lines = event.text.strip().split('\n')
    
    # Find the index of the last photon line
    photon_indices = [i for i, line in enumerate(event_lines) if line.strip().split()[0] == last_outgoing]
    insert_index = max(photon_indices) + 1 if photon_indices else 1

    # Insert proton data into the event at the appropriate index
    updated_event_lines = event_lines[:insert_index] + [proton1_data, proton2_data] + event_lines[insert_index:]
    event.text = '\n'.join(updated_event_lines)

def parse_and_modify_lhe_file(input_file_path, output_file_path):
    tree = ET.parse(input_file_path)
    root = tree.getroot()

    initial_proton_momentum = (7000.0, 0.0, 0.0, 7000.0)

    events = root.findall('event')
    total_events = len(events)

    with tqdm(total=total_events, desc="Processing events", unit="event", ncols=100) as pbar:
        for event in events:
            lines = event.text.strip().split('\n')
            
            photon1_data = None
            photon2_data = None
            
            for line in lines:
                # Check if the line meets the criteria for photon data
                if line.strip().split()[0] == pdgId and line.strip().split()[1] == '-1':
                    photon1_data = line
                elif line.strip().split()[0] == '-'+pdgId and line.strip().split()[1] == '-1':
                    photon2_data = line
                    break  # Exit the loop after the second detection

            if photon1_data and photon2_data:
                # Process photon data
                proton1_momentum, proton2_momentum = calculate_outgoing_protons([photon1_data, photon2_data], initial_proton_momentum)
                add_outgoing_protons_to_event(event, proton1_momentum, proton2_momentum)

            pbar.update(1)

    with open(output_file_path, 'wb') as f:
        # Add XML declaration and root tag
        f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(b'<LesHouchesEvents version="1.0">\n')
        
        for child in root:
            event_str = ET.tostring(child, encoding='utf-8')
            f.write(event_str)
        
        # Close the root tag
        f.write(b'</LesHouchesEvents>\n')

def add_line_break_after_event(output_file_path):
    with open(output_file_path, 'r') as file:
        content = file.readlines()

    # Modify the content
    for i, line in enumerate(content):
        if line.strip().startswith("<event>"):
            # Split the line into two lines with a line break after "<event>"
            content[i] = line.replace("<event>", "<event>\n")

    # Write the modified content back to the file
    with open(output_file_path, 'w') as file:
        file.writelines(content)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path_to_lhe_file> <pdgId-no-sign> <pdgId-last-outgoing>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    pdgId = sys.argv[2]
    last_outgoing = sys.argv[3]
    output_file_path = os.path.join(os.path.dirname(input_file_path), f"protonsADD_{os.path.basename(input_file_path)}")
    parse_and_modify_lhe_file(input_file_path, output_file_path)
    add_line_break_after_event(output_file_path)

