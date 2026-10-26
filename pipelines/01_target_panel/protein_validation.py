"""Independent GSE267552 RNA/ADT availability analysis, with donors as units.

Primary cells are selected from original author metadata, independently of the
existing CD127 table; the old table is used solely for a membership/value audit.
No ADT positivity threshold is defined. This is a paired tissue phenotype study,
not a therapy-response study. Tissue and technical lane allocation differ.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import platform

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

TARGETS = ['PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','BTLA','CD160','VSIR',
           'ENTPD1','NT5E','ADORA2A','TNFRSF9','TNFRSF4','TNFRSF18','ICOS',
           'CD27','CD28','CXCR6','ITGAE','IL7R','IL2RA','CCR8','CD38','SLAMF6']
PROTEINS = dict(PDCD1='Hu.CD279', LAG3='Hu.CD223', TIGIT='Hu.TIGIT',
                CTLA4='Hu.CD152', BTLA='Hu.CD272', ENTPD1='Hu.CD39',
                NT5E='Hu.CD73', TNFRSF9='Hu.CD137', TNFRSF4='Hu.CD134',
                ICOS='HuMsRt.CD278', CD27='Hu.CD27', CD28='Hu.CD28',
                ITGAE='Hu.CD103', IL7R='Hu.CD127', IL2RA='Hu.CD25',
                CD38='Hu.CD38_HIT2', SLAMF6='Hu.CD352')
DONORS = ['BC07','BC08','BC09','BC10','BC12']
TISSUES = ['Tumour','Peri_Healthy']
STATES = ['TEM','TEMRA','TEX','TRM']
SEED = 20261025


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def decode(values):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in values])


def exact_signflip(difference):
    x = np.asarray(difference, float)
    assert np.isfinite(x).all()
    signs = np.array(list(itertools.product([-1., 1.], repeat=len(x))))
    null = np.abs((signs * x).mean(axis=1))
    return float(np.mean(null >= abs(x.mean()) - 1e-12))


def rho(x, y):
    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan
    return float(stats.spearmanr(x, y).statistic)


def write_csv(frame, directory, name):
    frame.to_csv(directory / name, index=False, float_format='%.12g')


def analyze(source, out):
    data_dir = out / 'data'
    results_dir = out / 'results'
    for p in (data_dir, results_dir):
        p.mkdir(parents=True, exist_ok=True)
    metadata_path = source / 'raw/GSE267552_GEO_metadata.csv.gz'
    prior_path = source / 'analysis/primary_CD8_RNA_protein_TCR.csv.gz'
    metadata = pd.read_csv(metadata_path).rename(columns={'Unnamed: 0': 'cell_id',
                                                         'CellType_major': 'author_state'})
    assert metadata.cell_id.is_unique
    metadata['lane'] = metadata['orig.ident'].str.extract(r'^Lane ([ABC])$', expand=False)
    assert metadata.lane.notna().all()
    metadata['barcode'] = metadata.cell_id.str.replace(r'_\d+$', '', regex=True)
    assert not metadata.duplicated(['lane', 'barcode']).any()
    selected = metadata[metadata.donor.isin(DONORS) & metadata.tissue.isin(TISSUES)
                        & metadata.author_state.isin(STATES)].copy()
    prior = pd.read_csv(prior_path)
    assert len(selected) == len(prior) == 5061
    assert set(selected.cell_id) == set(prior.cell_id)
    check = selected.merge(prior, on='cell_id', suffixes=('', '_prior'), validate='one_to_one')
    for a, b in [('donor', 'patient'), ('tissue', 'tissue_prior'),
                 ('lane', 'lane_prior'), ('barcode', 'barcode_prior'),
                 ('author_state', 'author_state_prior')]:
        assert check[a].eq(check[b]).all(), (a, b)
    blocks, panel_rows, lane_audit = [], [], []
    sources = [metadata_path, prior_path, source / 'analysis/antibody_panels.csv']
    for lane in 'ABC':
        h5path = source / f'raw/GSE267552_lane{lane}-CITE-GEX.h5'
        sources.append(h5path)
        q = selected[selected.lane.eq(lane)].copy()
        with h5py.File(h5path, 'r') as h:
            m = h['matrix']
            names = decode(m['features/name'][:])
            types = decode(m['features/feature_type'][:])
            ids = decode(m['features/id'][:])
            barcodes = decode(m['barcodes'][:])
            assert len(barcodes) == len(set(barcodes))
            indices = pd.Index(barcodes).get_indexer(q.barcode)
            assert (indices >= 0).all()
            x = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]),
                                  shape=tuple(m['shape'][:]))[:, indices]
        rna_mask = types == 'Gene Expression'
        adt_mask = types == 'Antibody Capture'
        rna = x[rna_mask, :].tocsc()
        gene_names = names[rna_mask]
        adt = x[adt_mask, :].toarray().T.astype(float)
        adt_names, adt_ids = names[adt_mask], ids[adt_mask]
        assert len(set(adt_names)) == len(adt_names)
        iso_mask = np.array([name.startswith('Isotype_') for name in adt_names])
        assert iso_mask.sum() == 7
        isotype_background = np.log1p(adt[:, iso_mask]).mean(axis=1)
        log_adt_mean = np.log1p(adt).mean(axis=1)
        total_adt = adt.sum(axis=1)
        total_rna = np.asarray(rna.sum(axis=0)).ravel()
        q['rna_total_umi'] = total_rna
        q['rna_n_genes'] = np.diff(rna.indptr)
        mt = np.array([name.startswith('MT-') for name in gene_names])
        q['rna_mito_fraction'] = np.asarray(rna[mt, :].sum(axis=0)).ravel() / np.maximum(total_rna, 1)
        q['numeric_qc'] = (q.rna_total_umi >= 800) & (q.rna_n_genes >= 250) & (q.rna_mito_fraction <= .1)
        q['adt_total'] = total_adt
        q['isotype_log_background'] = isotype_background
        for name, vector in zip(adt_names[iso_mask], adt[:, iso_mask].T):
            q[f'{name}_adt_raw'] = vector.astype(int)
        for target_index, target in enumerate(TARGETS):
            if target_index % 8 == 0:
                q = q.copy()  # consolidate blocks before adding assay columns
            gene_rows = gene_names == target
            assert gene_rows.sum() >= 1, (lane, target)
            raw = np.asarray(rna[gene_rows, :].sum(axis=0)).ravel()
            q[f'{target}_rna_raw'] = raw.astype(int)
            q[f'{target}_rna_logcp10k'] = np.log1p(10000 * raw / np.maximum(total_rna, 1))
            entry = dict(lane=lane, target=target, rna_feature_rows=int(gene_rows.sum()),
                         rna_feature_ids=';'.join(ids[rna_mask][gene_rows]),
                         protein_assayed=target in PROTEINS,
                         adt_feature_name=PROTEINS.get(target), adt_feature_id=None)
            if target in PROTEINS:
                ai = np.flatnonzero(adt_names == PROTEINS[target])
                assert len(ai) == 1, (lane, target, ai)
                raw_adt = adt[:, ai[0]]
                q[f'{target}_adt_raw'] = raw_adt.astype(int)
                q[f'{target}_adt_corrected'] = np.log1p(raw_adt) - isotype_background
                q[f'{target}_adt_clr'] = np.log1p(raw_adt) - log_adt_mean
                q[f'{target}_adt_logcp10k'] = np.log1p(10000 * raw_adt / np.maximum(total_adt, 1))
                entry['adt_feature_id'] = adt_ids[ai[0]]
            panel_rows.append(entry)
        blocks.append(q)
        lane_audit.append(dict(lane=lane, raw_matrix_cells=len(barcodes),
                               author_cells=int(metadata.lane.eq(lane).sum()),
                               primary_cells=len(q), adt_features=int(adt_mask.sum()),
                               isotypes=int(iso_mask.sum()), all_primary_barcodes_matched=True))
        print(json.dumps(lane_audit[-1]), flush=True)
    cells = pd.concat(blocks, ignore_index=True).sort_values('cell_id').reset_index(drop=True)
    assert cells.cell_id.is_unique and len(cells) == 5061
    audit_prior = cells.merge(prior, on='cell_id', suffixes=('', '_prior'), validate='one_to_one')
    audit_errors = {}
    for a, b in [('IL7R_adt_raw', 'CD127_raw'), ('IL7R_adt_corrected', 'CD127_protein'),
                 ('IL7R_rna_raw', 'IL7R_counts'), ('IL7R_rna_logcp10k', 'IL7R_RNA'),
                 ('ENTPD1_adt_raw', 'CD39_raw'), ('NT5E_adt_raw', 'CD73_raw'),
                 ('isotype_log_background', 'isotype_background'), ('rna_total_umi', 'GEX_UMI')]:
        audit_errors[f'{a}_vs_{b}'] = float(np.max(np.abs(audit_prior[a] - audit_prior[b])))
    assert max(audit_errors.values()) < 1e-10, audit_errors
    write_csv(cells, data_dir, 'protein_primary_cells.csv.gz')
    write_csv(pd.DataFrame(panel_rows), results_dir, 'protein_assay_coverage.csv')
    write_csv(pd.DataFrame(lane_audit), results_dir, 'protein_lane_audit.csv')
    denominators = cells.groupby(['donor', 'tissue', 'lane', 'author_state']).size().rename('n_cells').reset_index()
    write_csv(denominators, results_dir, 'protein_cell_denominators.csv')
    # Per-lane and per-donor/tissue results retain the technical allocation.
    lane_summaries, donor_summaries, correlations = [], [], []
    for keys, frame, dest in [(['donor','tissue','lane'], cells, lane_summaries),
                              (['donor','tissue'], cells, donor_summaries)]:
        for key, z in frame.groupby(keys, sort=True):
            for target in TARGETS:
                row = dict(zip(keys, key))
                row.update(target=target, n_cells=len(z), n_numeric_qc=int(z.numeric_qc.sum()),
                           rna_detection_fraction=float(z[f'{target}_rna_raw'].gt(0).mean()),
                           rna_mean_logcp10k=float(z[f'{target}_rna_logcp10k'].mean()),
                           protein_assayed=target in PROTEINS,
                           n_lanes=int(z.lane.nunique()), lanes=';'.join(sorted(z.lane.unique())))
                if target in PROTEINS:
                    row.update(adt_mean_corrected=float(z[f'{target}_adt_corrected'].mean()),
                               adt_median_corrected=float(z[f'{target}_adt_corrected'].median()),
                               adt_mean_clr=float(z[f'{target}_adt_clr'].mean()),
                               adt_mean_logcp10k=float(z[f'{target}_adt_logcp10k'].mean()),
                               adt_raw_nonzero_fraction=float(z[f'{target}_adt_raw'].gt(0).mean()))
                dest.append(row)
    lane_summary = pd.DataFrame(lane_summaries)
    donor_summary = pd.DataFrame(donor_summaries)
    write_csv(lane_summary, results_dir, 'protein_lane_tissue_summary.csv')
    write_csv(donor_summary, results_dir, 'protein_donor_tissue_summary.csv')
    for keys in [['donor'], ['donor','tissue'], ['donor','tissue','lane']]:
        for key, z in cells.groupby(keys, sort=True):
            key = key if isinstance(key, tuple) else (key,)
            for target in PROTEINS:
                row = dict(zip(keys, key))
                row.update(target=target, grouping='_'.join(keys), n_cells=len(z),
                           rna_adt_rho=rho(z[f'{target}_rna_logcp10k'], z[f'{target}_adt_corrected']),
                           n_rna_detected=int(z[f'{target}_rna_raw'].gt(0).sum()))
                correlations.append(row)
    cor = pd.DataFrame(correlations)
    write_csv(cor, results_dir, 'protein_rna_correlations.csv')
    # Identical bootstrap donor indices across targets support reproducibility.
    bootstrap_indices = np.random.default_rng(SEED).integers(0, len(DONORS), (5000, len(DONORS)))
    paired_rows, effect_rows = [], []
    for mode in ['primary_cell_weighted', 'equal_technical_lane', 'numeric_qc', 'clr_sensitivity', 'adt_depth_sensitivity']:
        for target in PROTEINS:
            value_col = f'{target}_adt_corrected'
            if mode == 'clr_sensitivity':
                value_col = f'{target}_adt_clr'
            if mode == 'adt_depth_sensitivity':
                value_col = f'{target}_adt_logcp10k'
            z = cells[cells.numeric_qc] if mode == 'numeric_qc' else cells
            if mode == 'equal_technical_lane':
                values = z.groupby(['donor','tissue','lane'])[value_col].mean().groupby(['donor','tissue']).mean()
            else:
                values = z.groupby(['donor','tissue'])[value_col].mean()
            pivot = values.unstack('tissue').loc[DONORS]
            assert pivot[TISSUES].notna().all().all()
            differences = pivot.Tumour.to_numpy() - pivot.Peri_Healthy.to_numpy()
            ci = np.quantile(differences[bootstrap_indices].mean(axis=1), [.025, .975])
            for donor, difference in zip(DONORS, differences):
                paired_rows.append(dict(mode=mode, target=target, donor=donor,
                                        tumour_mean=float(pivot.loc[donor, 'Tumour']),
                                        peritumour_mean=float(pivot.loc[donor, 'Peri_Healthy']),
                                        tumour_minus_peritumour=float(difference)))
            effect_rows.append(dict(mode=mode, target=target, n_donors=len(DONORS),
                                    mean_tumour_minus_peritumour=float(differences.mean()),
                                    median_tumour_minus_peritumour=float(np.median(differences)),
                                    bootstrap_ci_low=float(ci[0]), bootstrap_ci_high=float(ci[1]),
                                    exact_signflip_p=exact_signflip(differences),
                                    positive_donors=int((differences > 0).sum()),
                                    negative_donors=int((differences < 0).sum())))
    effects = pd.DataFrame(effect_rows)
    effects['bh_q'] = effects.groupby('mode').exact_signflip_p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
    write_csv(pd.DataFrame(paired_rows), results_dir, 'protein_paired_donor_differences.csv')
    write_csv(effects, results_dir, 'protein_paired_effects.csv')
    figure = donor_summary.copy()
    primary_effect = effects[effects['mode'].eq('primary_cell_weighted')].drop(columns='mode')
    figure = figure.merge(primary_effect, on='target', how='left', validate='many_to_one')
    donor_cor = cor[cor.grouping.eq('donor')][['donor','target','rna_adt_rho']].rename(columns={'rna_adt_rho': 'donor_pooled_rna_adt_rho'})
    figure = figure.merge(donor_cor, on=['donor','target'], how='left', validate='many_to_one')
    write_csv(figure, results_dir, 'protein_figure5_source.csv')
    # A target-level view shows transcript detection and continuous antibody signal.
    concordance = []
    for target in PROTEINS:
        z = donor_summary[donor_summary.target.eq(target)]
        t = z[z.tissue.eq('Tumour')]
        p = z.pivot(index='donor', columns='tissue', values='rna_detection_fraction').loc[DONORS]
        rna_differences = p.Tumour - p.Peri_Healthy
        protein_differences = np.array([r['tumour_minus_peritumour'] for r in paired_rows
                                       if r['mode'] == 'primary_cell_weighted' and r['target'] == target])
        pooled_rhos = donor_cor[donor_cor.target.eq(target)].donor_pooled_rna_adt_rho
        tissue_rhos = cor[cor.grouping.eq('donor_tissue') & cor.target.eq(target)].rna_adt_rho
        concordance.append(dict(target=target,
                               tumour_median_donor_rna_detection=float(t.rna_detection_fraction.median()),
                               tumour_median_donor_adt_median=float(t.adt_median_corrected.median()),
                               mean_rna_detection_tumour_minus_peritumour=float(rna_differences.mean()),
                               mean_protein_tumour_minus_peritumour=float(protein_differences.mean()),
                               tissue_effect_same_sign=bool(np.sign(rna_differences.mean()) == np.sign(protein_differences.mean())),
                               within_donor_median_rho=float(pooled_rhos.median()),
                               within_donor_positive_rhos=int(pooled_rhos.gt(0).sum()),
                               within_donor_finite_rhos=int(pooled_rhos.notna().sum()),
                               within_donor_tissue_median_rho=float(tissue_rhos.median()),
                               within_donor_tissue_positive_rhos=int(tissue_rhos.gt(0).sum()),
                               within_donor_tissue_finite_rhos=int(tissue_rhos.notna().sum())))
    concordance = pd.DataFrame(concordance)
    write_csv(concordance, results_dir, 'protein_rna_concordance.csv')
    # Known identities test the actual statistical functions used above.
    assert exact_signflip([1, 1, 1, 1, 1]) == .0625
    assert exact_signflip([0, 0, 0, 0, 0]) == 1
    assert exact_signflip([1, -1, 0, 0, 0]) == 1
    assert rho([1, 2, 3], [3, 2, 1]) == -1
    assert np.isnan(rho([0, 0, 0], [1, 2, 3]))
    assert len(donor_summary) == 24 * 5 * 2
    assert len(primary_effect) == 17
    assert effects.exact_signflip_p.ge(.0625).all()
    assert not donor_summary.loc[~donor_summary.protein_assayed, 'adt_mean_corrected'].notna().any()
    summary = dict(study='GSE267552 breast CD8 RNA/ADT availability', source_root=str(source),
                   primary_cells=len(cells), independent_donors=DONORS, n_donors=5,
                   tumour_cells=int(cells.tissue.eq('Tumour').sum()),
                   peritumour_cells=int(cells.tissue.eq('Peri_Healthy').sum()),
                   author_state_counts=cells.author_state.value_counts().to_dict(),
                   selection='Author-retained TEM/TEMRA/TEX/TRM, five paired donors, Tumour/Peri_Healthy; no RNA or ADT target gate',
                   targets=TARGETS, protein_assayed=list(PROTEINS),
                   protein_not_assayed=[target for target in TARGETS if target not in PROTEINS],
                   normalization=dict(rna='log1p(10000*gene count / full GEX UMI)',
                                      adt='log1p(target ADT) - mean(log1p(seven isotype ADTs))',
                                      clr_sensitivity='log1p(target ADT) - mean(log1p(all 137 ADTs))',
                                      protein_positivity_threshold=None),
                   inference=dict(unit='independent donor', contrast='Tumour minus matched Peri_Healthy',
                                  primary='cell-weighted donor/tissue mean corrected ADT, equally weighted donors',
                                  exact_test='All 32 donor sign flips; two-sided absolute mean difference',
                                  minimum_two_sided_p=.0625, multiplicity='BH separately for 17 targets in each analysis mode',
                                  bootstrap='5000 donor bootstrap draws, percentile 95% CI; five donors give coarse uncertainty',
                                  seed=SEED, correlation='Descriptive per donor and per donor/tissue/lane Spearman; no cell-level significance test'),
                   interpretation='Continuous antibody signal supports measured protein availability; raw ADT nonzero counts and positive corrected signals do not define protein-positive cells. Tumour/peritumour contrasts are descriptive because tissue allocation differs across technical lanes. No treatment response or functional efficacy is measured.',
                   audit=dict(independent_metadata_membership_matches_prior=True, one_to_one_lane_barcode=True,
                              value_audit_max_abs_errors=audit_errors, known_statistical_identities_passed=True,
                              missing_proteins_not_encoded_zero=True),
                   positive_median_within_donor_rna_protein_correlations=concordance.loc[concordance.within_donor_median_rho.gt(0), 'target'].tolist(),
                   tissue_effect_same_sign_targets=concordance.loc[concordance.tissue_effect_same_sign, 'target'].tolist(),
                   all_five_protein_tissue_differences_positive=primary_effect.loc[primary_effect.positive_donors.eq(5),'target'].tolist(),
                   all_five_protein_tissue_differences_negative=primary_effect.loc[primary_effect.negative_donors.eq(5),'target'].tolist(),
                   number_bh_q_below_05=int(primary_effect.bh_q.lt(.05).sum()),
                   versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                                 scipy=scipy.__version__, h5py=h5py.__version__))
    (results_dir / 'protein_analysis_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    sources.append(Path(__file__))
    manifest = dict(source_files=[dict(path=str(path), size_bytes=path.stat().st_size, sha256=digest(path)) for path in sources],
                    output_files=[dict(path=str(path.relative_to(out)), size_bytes=path.stat().st_size, sha256=digest(path))
                                  for path in sorted(data_dir.glob('protein_*')) + sorted(results_dir.glob('protein_*'))
                                  if path.name != 'protein_source_manifest.json'])
    (results_dir / 'protein_source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('D:/CD8_Immunotherapy_20261023/analysis_extensions/GSE267552_protein_validation'))
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    analyze(args.source, args.out)
