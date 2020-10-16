# ======== merge data =============================================================

import pandas as pd
import numpy as np

chain_A = pd.read_csv('MEF2C_DNA_hydrate_A_scores.csv').drop_duplicates('mutation', keep='last')
chain_B = pd.read_csv('MEF2C_DNA_hydrate_B_scores.csv').drop_duplicates('mutation', keep='last')
chain_A['resi'] = chain_A.mutation.apply(lambda v: int(v[1:-1]) if v[1:-1].isdigit() else float('inf'))
chain_A = chain_A.sort_values('resi')
chain_B['resi'] = chain_B.mutation.apply(lambda v: int(v[1:-1]) if v[1:-1].isdigit() else float('inf'))
chain_B = chain_B.sort_values('resi')
# .set_index('mutation')

combo = pd.concat([chain_A, chain_B]).groupby('mutation').max().sort_values('resi')
combo = combo.assign(homodimer_ddG=np.max(combo[['homodimer_B_interface_ddG', 'homodimer_B_interface_ddG']].values, axis=1) )
combo = combo[['complex_ddG', 'DNA_interface_ddG', 'homodimer_ddG', 'FA_RMSD', 'CA_RMSD']]

# ======== add proximity data =============================================================


import pymol2

closeness = []

with pymol2.PyMOL() as pymol:
    pymol.cmd.load('MEF2C_long.rx.pdb', 'MEF2C')
    combo = combo.assign(
                    is_dna_within_4Å=combo.index.map(lambda mutation: pymol.cmd.select(selection=f'(MEF2C and chain A+B and resi {mutation[1:-1]} around 4) and chain L+K',name='none') != 0),
                    is_dimer_within_4Å=combo.index.map(lambda mutation: pymol.cmd.select(selection=f'(MEF2C and chain A and resi {mutation[1:-1]} around 4) and chain B', name='none') +\
                                                                        pymol.cmd.select(selection=f'(MEF2C and chain B and resi {mutation[1:-1]} around 4) and chain A', name='none') != 0)
                   )

# ========== add grouping data ===========================================================

import json

grouping = json.load(open('grouping.json'))
combo = combo.assign(group=combo.index.map(lambda x: ','.join(grouping[x])))

# ========== scatter plot =============================================================================

def colorise(value):
    if 'DDD' in value:
        return 'green'
    elif 'ClinVar' in value:
        return 'blue'
    elif 'gnomAD_MEF2C' in value:
        return 'dimgrey'
    else:
        return 'gainsboro'


def symbolify(row):
    if row.is_dna_within_4Å:
        return 'x'
    #     elif row.is_dimer_within_4Å:
    #         return 'diamond'
    else:
        return 'circle'


import plotly.graph_objects as go


series = []
for group in [
             'gnomAD_MEF2A',
             'gnomAD_MEF2A,gnomAD_MEF2B',
             'gnomAD_MEF2A,gnomAD_MEF2D',
             'gnomAD_MEF2B',
             'gnomAD_MEF2B,gnomAD_MEF2D',
             'gnomAD_MEF2D',
             'gnomAD_MEF2C',
             'gnomAD_MEF2C,gnomAD_MEF2B',
                'ClinVar',
                 'DDD']:
    sub = combo.loc[combo.group == group]
    series.append(go.Scatter(
                            x=sub.complex_ddG.values,
                            y=sub.DNA_interface_ddG.values,
                            mode='markers+text' if group in ['DDD', 'ClinVar'] else 'markers',
                            marker=dict(color=sub.group.apply(colorise).values,
                                        symbol=sub.apply(symbolify, axis=1).values),
                            name=group,
                            text=sub.index.values,
                            textposition='top right'
                            )
                 )

fig = go.Figure(data=series, layout=go.Layout(title='∆∆G of MEF2C mutations<br>(saltire: within 4 Å of DNA, circle: farther than 4 Å of DNA)',
                                             xaxis=dict(title='∆∆G complex (kcal/mol)'),
                                             yaxis=dict(title='∆∆G DNA interface (kcal/mol)')))
fig.show(renderer="notebook")

# ====== heatmap ====================================================================

import plotly.graph_objects as go

def get_variant_order(df:pd.DataFrame) -> list:
    return sorted(df.columns.values, key=lambda v: int(v[1:-1]))

def get_annotation_for_heatmap(df:pd.DataFrame) -> list:
    d = df.to_numpy()
    annotations = []
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            if abs(d[i,j]) > 2:
                text = f"<b>{d[i,j]:.1f}</b>"
            else:
                text = f"{d[i,j]:.1f}"
            annotations.append(dict(
                                    x=j,
                                    y=i,
                                    xref="x",
                                    yref="y",
                                    text=text,
                                    showarrow=False
                                    ))
    return annotations

#ddG = scores.pivot_table(index=['model'], columns=['mutation'], values=f'{n}_ddG').round(1)
ddG = data[['mutation','complex_ddG','DNA_interface_ddG','homodimer_interface_ddG']]
ddG = ddG.set_index('mutation').astype(float) #.transpose()
ddG = ddG.transpose()[reversed(get_variant_order(ddG.transpose()))].transpose()

fig = go.Figure(data=go.Heatmap(
    z=ddG,
    x=ddG.columns,
    y=list(map(lambda v: f'{v} ({",".join(grouping[v])})', ddG.index.values)),
    hoverongaps=False,
    colorscale='temps'))
fig.data[0].update(zmin=-10, zmax=10)
fig.update_layout(title=f'ΔΔG (fold and interface)',
                  annotations=get_annotation_for_heatmap(ddG),
                  xaxis = {'title': 'Variant'},
                  yaxis = {'title': 'Model'})
#fig.write_image(f'heatmap.png', scale=3)
fig.show()

