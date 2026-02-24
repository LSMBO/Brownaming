import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
import re

def generate_combined_figure(stats, output_file='combined_results.png'):
    total_nb_query = stats['Step 1']['nb_query']
    stats_dict = {
        'step': ['0'],
        'rank': ['NA'],
        'taxon_name': ['NA'],
        'taxon_id': ['NA'],
        'dbsize': ['NA'],
        'prots_with_hit': [0],
        '%_prots_with_hit': [0],
        'elapsed_time_min': [0.0],
        'elapsed_time_str': ['0m']
    }

    
    for step in stats:
        if step.startswith("Step"):
            elapsed_time_min = float(stats[step]['elapsed_time'])
            stats_dict['step'].append(step.replace("Step ", ""))
            stats_dict['rank'].append(stats[step]['rank'])
            stats_dict['taxon_name'].append(stats[step]['taxon_name'])
            stats_dict['taxon_id'].append(stats[step]['taxon_id'])
            stats_dict['dbsize'].append(stats[step]['dbsize'])
            stats_dict['prots_with_hit'].append(stats[step]['prots_with_hit'])
            stats_dict['%_prots_with_hit'].append(round(100 * stats[step]['prots_with_hit'] / total_nb_query))
            stats_dict['elapsed_time_min'].append(elapsed_time_min)
            stats_dict['elapsed_time_str'].append(format_elapsed_time(elapsed_time_min))
    
    # si entre 0 et 20 -> height=10 et ratio=40-60
    # si entre 21 et 30 -> ??
    fig = plt.figure(figsize=(12, 10))
    
    gs = fig.add_gridspec(2, 1, height_ratios=[0.6, 0.4])
    
    ax_table = fig.add_subplot(gs[0, 0])
    create_table(stats_dict, ax_table)
    
    ax_plot = fig.add_subplot(gs[1, 0])
    create_plot(stats_dict, ax_plot)
       
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Combined figure saved as {output_file}")
    
    return fig

def format_elapsed_time(minutes, pos=None):
    if minutes >= 60:
        hours = int(minutes // 60)
        remaining_minutes = minutes % 60
        if remaining_minutes.is_integer():
            return f"{hours}h {remaining_minutes:.0f}m"
        return f"{hours}h {remaining_minutes:.1f}m"
    else:
        if minutes.is_integer():
            return f"{minutes:.0f}m"
        return f"{minutes:.1f}m"
    
def create_table(stats_dict, ax):
    table_data = []
    headers = [
        "Step", 
        "Taxonomic Rank", 
        "Taxon Name", 
        "Taxon ID", 
        "Database size\n(#seqs)", 
        "# Named proteins",
        "Cumulative\nelapsed time"
    ]

    for i in range(len(stats_dict['step'])):
        table_data.append([
            stats_dict['step'][i],
            stats_dict['rank'][i],
            stats_dict['taxon_name'][i],
            stats_dict['taxon_id'][i],
            stats_dict['dbsize'][i],
            f"{stats_dict['prots_with_hit'][i]} ({stats_dict['%_prots_with_hit'][i]}%)",
            stats_dict['elapsed_time_str'][i]            
        ])

    ax.axis('tight')
    ax.axis('off')
    table = ax.table(
        cellText=table_data,
        colLabels=headers,
        cellLoc='left',
        loc='center',
        colWidths=[0.1, 0.2, 0.3, 0.1, 0.2, 0.2, 0.2]
    )
    table.scale(1, 1.5)
    
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('darkgreen')
            cell.set_height(cell.get_height() * 1.5)
        else:
            if row % 2:
                cell.set_facecolor('lightgrey')    
    return table


def create_plot(stats_dict, ax):
    color = 'tab:green'
    ax.plot(stats_dict['step'], stats_dict['%_prots_with_hit'], 'o-', color=color, label='Named proteins (%)')
    ax.set_ylabel('Named proteins (%)', color=color)
    ax.tick_params(axis='y', labelcolor=color)

    ax2 = ax.twinx()
    color = 'tab:blue'
    ax2.plot(stats_dict['step'], stats_dict['elapsed_time_min'], 'o-', color=color, label=f'Cumulative Elapsed Time')
    ax2.yaxis.set_major_formatter(FuncFormatter(format_elapsed_time))
    
    ax2.set_ylabel('Cumulative Elapsed Time', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    
    ax.set_xlabel('Step (refer to table above)')

    return ax
