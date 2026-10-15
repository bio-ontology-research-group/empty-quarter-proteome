"""Extract bounded family claims from pinned compact evidence, not a raw-data rerun."""
import csv
import json
import statistics
from pathlib import Path

BASE = Path(__file__).resolve().parent
REPRO = BASE.parents[1]
PREFIX = 'darkmatter/family_followup_20261015/'
OLD = '../family_discovery_20261015/'
SOURCES = set()


def table(path, **filters):
    with open(BASE / path) as handle:
        return [r for r in csv.DictReader(handle, delimiter='\t')
                if all(r[k] == str(v) for k, v in filters.items())]


def one(path, field, **filters):
    data = table(path, **filters)
    assert len(data) == 1, (path, filters, len(data))
    return float(data[0][field])


def claims():
    result = []

    def add(key, desc, value, source, fn, tol=0):
        source_path = str((BASE / source).resolve().relative_to(REPRO))
        SOURCES.add(source_path)
        result.append(('family.' + key, 'Family follow-up', desc, value,
                       'FROZEN', 'reproducibility/' + source_path, fn, tol))

    meta = 'results/f3_comparison/taxonomy.tsv'
    full = OLD + 'root_audit/f345_comparison/Helio/metadata.tsv'
    add('helio.complete', 'Complete EQ heliorhodopsin profile candidates', 65, full,
        lambda: len(table(full, source='EQ')))
    for source, n in [('EQ', 45), ('SMAG', 35), ('GEM', 73), ('curated_reference', 25)]:
        add('helio.site.' + source, 'Aligned-lysine-conserving HeR records: ' + source, n, meta,
            lambda source=source: len(table(meta, source=source)))
    finer = 'results/f3_finer/taxonomy_unique.tsv'
    for source, n in [('EQ', 22), ('GEM', 10), ('SMAG', 2)]:
        add('helio.sacch.' + source, 'Unique Saccharimonadia HeRs: ' + source, n, finer,
            lambda source=source: len(table(finer, source=source, **{'class': 'Saccharimonadia'})))
    pubsum = 'results/published_helio/summary.json'
    add('helio.published_n', 'Profile/site-qualified published HeRs', 1566, pubsum,
        lambda: json.load(open(BASE / pubsum))['site_conserving_pass'])
    hits = 'results/published_helio/closest_published_hits.tsv'
    add('helio.published_hits', 'EQ HeRs with coverage-qualified published hits', 45, hits,
        lambda: sum(bool(r['target']) for r in table(hits)))
    add('helio.identity', 'Median identity of best-bit-score published hit', 57.3, hits,
        lambda: statistics.median(float(r['identity']) for r in table(hits)), .05)
    add('helio.patesci_n', 'EQ Patescibacteriota HeR loci', 24, hits,
        lambda: len(table(hits, phylum='Patescibacteriota')))
    add('helio.patesci_identity', 'Patescibacteriota median published hit identity', 39.4, hits,
        lambda: statistics.median(float(r['identity']) for r in table(hits, phylum='Patescibacteriota')), .05)
    rare = 'results/published_clusters/rarefaction.tsv'
    for kind, value, clusters in [('full', 7.36, 13), ('core', 7.33, 12)]:
        filters = dict(rank='class', taxon='Saccharimonadia', comparison='EQ+GEM', identity='50', kind=kind)
        for source, v in [('EQ', value), ('GEM', 0)]:
            add('helio.rarefaction.' + kind + '.' + source, 'Mean added 50%-identity clusters, 10 sequences', v, rare,
                lambda filters=filters, source=source: one(rare, 'mean_external_clusters', source=source, **filters), .005)
        add('helio.extra_clusters.' + kind, 'All EQ Saccharimonadia clusters absent expanded reference', clusters, rare,
            lambda filters=filters: one(rare, 'clusters_absent_expanded_reference', source='EQ', **filters))
    sensitivity = 'results/expanded_analysis/read_supported_rarefaction.tsv'
    add('helio.read_rarefaction', 'Mean core clusters in read-supported EQ sensitivity', 7.31, sensitivity,
        lambda: one(sensitivity, 'mean_external_clusters', kind='core', identity='50', source='EQ'), .005)
    dist = 'results/expanded_analysis/native_distances.tsv'
    add('helio.GEM_overlap', 'GEM records exactly present in published panel', 8, dist,
        lambda: sum(int(r['exact_sequence_in_published_panel']) for r in table(dist, source='GEM')))
    projection = 'results/expanded_analysis/summary.json'
    for key, n in [('projection_unique_sequences', 1726), ('input_records', 1744)]:
        add('helio.' + key, 'Expanded ProstT5 ' + key, n, projection,
            lambda key=key: json.load(open(BASE / projection))[key])
    read = 'results/read_support/loci.tsv'
    for track, n in [('F3_Helio', 43), ('F4_EctCD', 5)]:
        add('reads.' + track, 'Loci passing depth/breadth rule: ' + track, n, read,
            lambda track=track: sum(int(r['coverage_support_gate']) for r in table(read, track=track)))
    for key, fun, n in [('min', min, 4), ('max', max, 19)]:
        add('fusion.templates.' + key, key + ' junction-spanning templates', n, read,
            lambda fun=fun: fun(int(r['junction_spanning_templates']) for r in table(read, track='F4_EctCD')))
    structures = 'results/structure_audit/structures.tsv'
    for field, threshold, value in [('mean_CA_pLDDT', 70, 93.13), ('query_TM', .5, .865)]:
        add('helio.structure_pass.' + field, 'EQ models above ' + field + ' threshold', 45, structures,
            lambda field=field, threshold=threshold: sum(float(r[field]) >= threshold for r in table(structures, source='EQ')))
        add('helio.structure_median.' + field, 'Median EQ ' + field, value, structures,
            lambda field=field: statistics.median(float(r[field]) for r in table(structures, source='EQ')), .005 if field == 'mean_CA_pLDDT' else .0005)
    cpf = OLD + 'root_audit/f1_full/candidates.tsv'
    add('cpf.genes', 'Complete CPF profile candidates', 914, cpf, lambda: len(table(cpf, partial='00')))
    add('cpf.mags', 'MAGs with complete CPF candidates', 775, cpf,
        lambda: len({r['mag'] for r in table(cpf, partial='00')}))
    cpfsum = OLD + 'root_audit/f1_full/summary.json'
    add('prodigal.proteins', 'Caller-specific family-screen Prodigal translations', 20451986, cpfsum,
        lambda: json.load(open(BASE / cpfsum))['proteins'])
    eco = 'results/ecology/median_ecology_tests.tsv'
    for trait, key, n, tol in [('CPF', 'sites', 37, 0), ('CPF', 'visits', 69, 0),
                              ('CPF', 'median_difference_pp', 8.47, .005),
                              ('CPF', 'p_Holm_six_contrasts', .00131, .000005),
                              ('EctABC_strict', 'median_difference_pp', .56, .005),
                              ('EctABC_strict', 'p_raw', .0974, .00005),
                              ('EctABC_strict', 'p_Holm_six_contrasts', .389, .0005)]:
        add('ecology.' + trait + '.' + key, trait + ' ' + key, n, eco,
            lambda trait=trait, key=key: one(eco, key, trait=trait), tol)
    dec = 'results/ecology/decomposition_summary.tsv'
    for metric, n in [('total_pp', 8.71), ('between_pp', 5.31), ('within_pp', 1.03), ('exclusive_pp', 2.38)]:
        add('cpf.decomposition.' + metric, 'Genus grouping CPF ' + metric, n, dec,
            lambda metric=metric: one(dec, 'mean', trait='CPF', rank='genus', metric=metric), .005)
    add('cpf.unclassified_mass', 'Mean mass in individually grouped genus-unclassified MAGs, percent', 25.7, dec,
        lambda: 100 * one(dec, 'mean', trait='CPF', rank='genus', metric='unclassified_mass'), .05)
    for metric, low, high in [('total_pp', 4.19, 13.20), ('between_pp', 1.21, 9.02),
                              ('exclusive_pp', .90, 3.89), ('within_pp', .68, 1.40)]:
        for field, value in [('mean_ci_low', low), ('mean_ci_high', high)]:
            add('cpf.decomposition.' + metric + '.' + field, 'CPF descriptive site-bootstrap interval', value, dec,
                lambda metric=metric, field=field: one(dec, field, trait='CPF', rank='genus', metric=metric), .005)
    f4 = 'results/f4_integrity/loci.tsv'
    for key, n, filters in [('complete', 940, {}), ('WEYH', 822, {'WEYH': '1'}),
                           ('ABC', 616, {'strict_ABC': '1'}), ('ABC_cov70', 612, {'strict_ABC_cov70': '1'}),
                           ('edge', 164, {'WEYH': '1', 'context_status': 'edge_censored'}),
                           ('internal', 42, {'WEYH': '1', 'context_status': 'internal_without_complete_ABC'})]:
        add('ect.' + key, 'EQ EctC ' + key, n, f4,
            lambda filters=filters: len(table(f4, source='EQ_MAG', **filters)))
    for field, n in [('strict_ABC', 608), ('strict_ABC_cov70', 604)]:
        add('ect.mags.' + field, 'MAGs carrying ' + field, n, f4,
            lambda field=field: len({r['mag'] for r in table(f4, source='EQ_MAG', **{field: '1'})}))
    add('ect.distinct_AB', 'Strict EQ loci with separate A/B ORFs', 609, f4,
        lambda: sum(r['complete_A'] != r['complete_B'] for r in table(f4, source='EQ_MAG', strict_ABC='1')))
    add('ect.dual_AB', 'Strict EQ loci with a single A/B-profile ORF', 7, f4,
        lambda: sum(r['complete_A'] == r['complete_B'] for r in table(f4, source='EQ_MAG', strict_ABC='1')))
    add('ect.fusion_nonWEYH', 'Complete EQ EctCD fusions outside WEYH criterion', 1, f4,
        lambda: len(table(f4, source='EQ_MAG', fusion_CD='1', WEYH='0')))
    add('cpf.without_K01669', 'CPF-carrier MAGs lacking K01669', 28, eco,
        lambda: one(eco, 'carriers', trait='CPF_without_K01669'))
    overlap = 'results/review_checks/published_overlap_unique.tsv'
    add('helio.sacch.GEM_overlap', 'Unique GEM Saccharimonadia exactly present in published panel', 2, overlap,
        lambda: sum(int(r['exact_published']) for r in table(overlap, source='GEM', **{'class': 'Saccharimonadia'})))
    free = 'results/review_checks/reference_free_rarefaction.tsv'
    for identity, source, value in [('50', 'EQ', 8.04), ('50', 'GEM', 4), ('80', 'EQ', 9.61), ('80', 'GEM', 10)]:
        add('helio.reference_free.' + identity + '.' + source, 'Mean total core clusters without references', value, free,
            lambda identity=identity, source=source: one(free, 'mean_total_clusters', kind='core', identity=identity, source=source), .0051)
    assembly = 'results/structure_audit/EQ_assembly_verification.tsv'
    add('helio.translation_matches', 'EQ HeR translations independently matching original contigs', 45, assembly,
        lambda: sum(int(r['translation_match']) for r in table(assembly)))
    return result
