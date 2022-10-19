import pandas as pd
from glob import glob


def get_best_solution_keys():
    return [
        'replicon', 'hit_id', 'gene_name',
        'hit_pos', 'model_fqn', 'sys_id', 'sys_loci', 'locus_num',
        'sys_wholeness', 'sys_score', 'sys_occ', 'hit_gene_ref',
        'hit_status', 'hit_seq_len', 'hit_i_eval', 'hit_score',
        'hit_profile_cov', 'hit_seq_cov', 'hit_begin_match',
        'hit_end_match', 'counterpart', 'used_in'
    ]


def get_hmmer_extract_keys():

    return ['hit_id', 'replicon_name', 'position_hit', 'hit_sequence_length', 'gene_name', 'i_eval',
            'score', 'profile_coverage', 'sequence_coverage', 'begin', 'end']


def is_file_empty(path):
    prev_line = ''
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                prev_line = line
            else:
                break
    if prev_line.startswith('# No Systems found'):
        return True
    return False


def read_best_solution(path):
    col_names = get_best_solution_keys()

    if is_file_empty(path):
        return pd.DataFrame(columns=col_names)
    else:
        return pd.read_table(path, comment='#')


def read_hmmer_extract(path):
    col_names = get_hmmer_extract_keys()

    if is_file_empty(path):
        return pd.DataFrame(columns=col_names)
    else:
        return pd.read_table(path, comment='#', names=col_names)


def gather_best_solutions(tmp_dir):
    min_df = pd.DataFrame(columns=get_best_solution_keys())
    return pd.concat([min_df] + [read_best_solution(i) for i in glob(f'{tmp_dir}/*/best_solution.tsv')]).reset_index(
        drop=True)


def gather_hmmer_results(tmp_dir):
    min_df = pd.DataFrame(columns=get_hmmer_extract_keys())
    return pd.concat([min_df] + [read_hmmer_extract(i) for i in glob(f'{tmp_dir}/*/hmmer_results/*.res_hmm_extract')]).reset_index(
        drop=True)


def export_tables(tmp_dir, outdir):
    genes = gather_best_solutions(tmp_dir)
    genes[['type', 'subtype']] = genes['model_fqn'].str.split('/', expand=True)[[2, 3]]
    genes.to_csv(f'{outdir}/defense_finder_genes.tsv', sep='\t', index=False)

    systems = genes.groupby('sys_id').agg(type=('type', lambda x: x.iloc[0]),
                                          subtype=('subtype', lambda x: x.iloc[0]),
                                          protein_in_syst=('hit_id', lambda x: ','.join(x)),
                                          genes_count=('hit_id', 'count'),
                                          name_of_profiles_in_sys=('gene_name', lambda x: ','.join(x))
                                          )

    systems.to_csv(f'{outdir}/defense_finder_systems.tsv', sep='\t')

    hmms = gather_hmmer_results(tmp_dir)
    hmms.to_csv(f'{outdir}/defense_finder_hmmer.tsv', sep='\t', index=False)