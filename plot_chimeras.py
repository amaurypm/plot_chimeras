#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl

def read_fasta(file_path):
    """Parses a FASTA file and returns a dictionary of sequences."""
    seqs = {}
    with open(file_path, 'r') as f:
        name, seq = "", []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name: seqs[name] = "".join(seq)
                name = line[1:]
                seq = []
            elif line:
                seq.append(line)
        if name: seqs[name] = "".join(seq)
    return seqs

def read_features(file_path):
    """
    Parses a CSV text file containing features.
    Supports both:
      - 3 columns: FeatureName, Start, End (Applies to ALL targets)
      - 4 columns: TargetName, FeatureName, Start, End (Applies to specific target)
    """
    features = {'ALL': []}
    if not file_path or not os.path.exists(file_path):
        return features
        
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = [p.strip() for p in line.strip().split(',')]
                
                # 3 columns: Applies to all targets
                if len(parts) == 3:
                    name, start, end = parts
                    features['ALL'].append((int(start), int(end), name))
                    
                # 4 columns: Applies to a specific target
                elif len(parts) >= 4:
                    target_name, name, start, end = parts[:4]
                    if target_name.upper() == 'ALL':
                        features['ALL'].append((int(start), int(end), name))
                    else:
                        if target_name not in features:
                            features[target_name] = []
                        features[target_name].append((int(start), int(end), name))
    return features

def infer_segments(target_seq, parentals_dict):
    """
    Uses Dynamic Programming (Viterbi algorithm) to find the sequence of parental 
    origins that minimizes the total number of segment switches (crossovers).
    Returns segments in the coordinates of the gapless target sequence.
    """
    valid_states = [] # valid_states[i] will hold all parentals matching the target at pos i
    
    # 1. Map out all matching parentals at each position
    gapless_pos = 1
    for i, res in enumerate(target_seq):
        if res == '-':
            continue # Skip alignment gaps in the target sequence
        
        matches = []
        for p_name, p_seq in parentals_dict.items():
            if i < len(p_seq) and p_seq[i] == res:
                matches.append(p_name)
                
        if not matches:
            matches = ["Unknown/Mut"]
            
        valid_states.append(matches)
        gapless_pos += 1
        
    total_length = gapless_pos - 1
    if total_length == 0:
        return [], 0
        
    # 2. Dynamic Programming Initialization
    # costs[p] stores the minimum number of switches to reach the current position ending in parental p
    # paths[p] stores the sequence of parental choices (the path) that achieved that minimum cost
    costs = {p: 0 for p in valid_states[0]}
    paths = {p: [p] for p in valid_states[0]}
    
    # 3. Dynamic Programming Iteration
    for i in range(1, total_length):
        current_states = valid_states[i]
        new_costs = {}
        new_paths = {}
        
        for p in current_states:
            min_cost = float('inf')
            best_prev_state = None
            
            # Calculate the cost to transition to state 'p' from any valid previous state 'q'
            for q, prev_cost in costs.items():
                # Cost is 0 if we stay on the same parental, 1 if we switch
                transition_cost = prev_cost + (0 if p == q else 1)
                
                if transition_cost < min_cost:
                    min_cost = transition_cost
                    best_prev_state = q
            
            new_costs[p] = min_cost
            new_paths[p] = paths[best_prev_state] + [p]
            
        costs = new_costs
        paths = new_paths
        
    # 4. Find the global optimal path at the end of the sequence
    # Find the final state that resulted in the lowest overall number of switches
    best_final_state = min(costs, key=costs.get)
    optimal_path = paths[best_final_state]
    
    # 5. Convert the optimal AA-by-AA path into a list of segmented blocks
    segments = []
    current_source = optimal_path[0]
    start_pos = 1
    
    for i in range(1, total_length):
        if optimal_path[i] != current_source:
            # We hit a switch, record the segment
            segments.append((start_pos, i, current_source))
            start_pos = i + 1
            current_source = optimal_path[i]
            
    # Append the final segment
    segments.append((start_pos, total_length, current_source))
    
    return segments, total_length

def main():
    parser = argparse.ArgumentParser(description="Generate combined chimeric protein sequence maps from an MSA.")
    parser.add_argument("-i", "--fasta", required=True, help="Input MSA FASTA file")
    parser.add_argument("-t", "--targets", required=True, nargs='+', help="List of target sequence IDs (space-separated)")
    parser.add_argument("-p", "--parentals", required=True, nargs='+', help="List of parental sequence IDs (space-separated)")
    parser.add_argument("-f", "--features", required=False, help="Optional CSV file with features")
    parser.add_argument("-o", "--output", required=False, default="combined_chimeras_map.png", help="Output PNG file path")
    
    args = parser.parse_args()

    # Define color palette
    cmap_colors = mpl.colormaps['tab20'].colors
    parental_colors = {pid: cmap_colors[i % len(cmap_colors)] for i, pid in enumerate(args.parentals)}
    parental_colors["Unknown/Mut"] = "#CCCCCC"

    # Parse files
    alignment = read_fasta(args.fasta)
    all_features = read_features(args.features)
    
    # Extract parentals
    parentals_dict = {}
    for pid in args.parentals:
        if pid in alignment:
            parentals_dict[pid] = alignment[pid]
        else:
            print(f"Error: Parental '{pid}' not found in the alignment file.")
            return

    # Process all targets first to gather data
    processed_targets = []
    max_length = 0
    for target in args.targets:
        if target not in alignment:
            print(f"Warning: Target '{target}' not found in the alignment. Skipping.")
            continue
            
        target_seq = alignment[target]
        segments, total_length = infer_segments(target_seq, parentals_dict)
        processed_targets.append({'name': target, 'segments': segments, 'length': total_length})
        max_length = max(max_length, total_length)

    if not processed_targets:
        print("No valid targets were processed. Exiting.")
        return

    # ==========================================
    # PLOTTING COMBINED FIGURE
    # ==========================================
    num_targets = len(processed_targets)
    fig, axes = plt.subplots(nrows=num_targets, ncols=1, figsize=(14, 2.5 * num_targets), sharex=True)
    
    if num_targets == 1:
        axes = [axes]

    global_legend_patches = {}

    for i, (target_data, ax) in enumerate(zip(processed_targets, axes)):
        target_name = target_data['name']
        
        # 1. Draw segments
        for start, end, source in target_data['segments']:
            width = end - start + 1
            color = parental_colors.get(source, '#CCCCCC') 
            
            rect = mpatches.Rectangle((start, 0), width, 1, facecolor=color, 
                                      edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)
            
            if source not in global_legend_patches:
                global_legend_patches[source] = mpatches.Patch(color=color, label=source)

        # 2. Draw Annotations (Combine 'ALL' features with target-specific features)
        target_features = all_features.get('ALL', []).copy()
        target_features.extend(all_features.get(target_name, []))
        
        y_levels = [1.3, 1.9, 2.5] 
        for j, (start, end, feat_name) in enumerate(target_features):
            level = y_levels[j % len(y_levels)]
            mid = (start + end) / 2
            
            ax.plot([start, end], [level, level], color='black', lw=1.2)
            ax.plot([start, start], [level - 0.1, level], color='black', lw=1.2)
            ax.plot([end, end], [level - 0.1, level], color='black', lw=1.2)
            ax.text(mid, level + 0.05, feat_name, ha='center', va='bottom', fontsize=9, fontweight='bold')

        # 3. Axis Formatting
        ax.set_xlim(0, max_length + 10)
        ax.set_ylim(-0.5, 3.5)
        ax.set_yticks([]) 
        
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
            
        ax.set_title(target_name, fontsize=12, fontweight='bold', loc='left', pad=10)

        if i == num_targets - 1:
            ax.set_xlabel("Amino Acid Position", fontsize=12, fontweight='bold')
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# 4. Global Legend and Layout
    # Organize the sequence tracks and add vertical breathing room between them
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.6) 

    # Place the legend at the very bottom edge of the figure
    fig.legend(handles=list(global_legend_patches.values()), loc='upper center', 
               bbox_to_anchor=(0.5, 0.0), ncol=min(len(global_legend_patches), 6), 
               frameon=False, fontsize=11)
    
    # bbox_inches='tight' will automatically expand the saved PNG to encompass the legend
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nSaved combined sequence map to: {args.output}")

if __name__ == "__main__":
    main()
