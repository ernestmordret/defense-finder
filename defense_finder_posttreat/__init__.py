from defense_finder_posttreat import gather_results


def run(tmp_dir, outdir):
    # Run analysis
    gather_results.export_tables(tmp_dir, outdir)

    print("Analysis results were written to {}".format(outdir))
