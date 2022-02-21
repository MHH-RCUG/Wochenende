import pandas as pd
from pathlib import Path
#from tqdm import tqdm
from model import compute_ptr
import click

ROOT_SAVE_PATH = 'fit_results/'


def evaluate_all_data(sample_folder_path, experiment_name, save_plots=True):
    """
    iterate through all csv files in sample_folder_path and compute
    growth rate for each
    :param sample_folder_path: folder containing samples
    :param experiment_name: folder name to save results to
    :return:
    """
    save_path = Path(ROOT_SAVE_PATH) / experiment_name / sample_folder_path.name
    # create result folder
    Path(save_path).mkdir(parents=True, exist_ok=True)
    growth_rates = []
    sample_names = []
    reads_all = []
    used_bins = []
    initial_bins = []
    error_codes = []
    bins_fit_error = []
    growth_rate_class = []
    bins_initial = 50

    # iterate through all csv files in directory
    files = list(sorted(sample_folder_path.glob('*.csv')))
    png_name = ""
    # Remove tqdm progress bar
    #for file in tqdm(files, total=len(files)):
    for file in files:
        png_name = file.stem
        # shorten filename
        png_name = png_name.replace(".fix.s", "").replace(".dup", "").replace(".ns", "").replace(".mm", "")
        png_name = png_name.replace(".mq20", "").replace(".mq30", "").replace(".calmd", "").replace(".filt", "")
        png_name = png_name.replace("chromosome", "").replace("complete", "").replace("genome", "")
        save_name = save_path / (png_name + '.png')
        sample_names.append(file.stem)
        gr, reads, bins, errs, fit_err = compute_ptr(file, n_windows=bins_initial,
                                                     save_name=save_name, save_plot=save_plots,
                                                     high_res_plot=False)
        used_bins.append(bins)
        initial_bins.append(bins_initial)
        growth_rates.append(gr)
        reads_all.append(reads)
        error_codes.append(errs)
        bins_fit_error.append(fit_err)

        if gr < 1 or len(errs) > 0:
            growth_rate_class.append('failed')
        elif 1.0 <= gr <= 1.1:
            growth_rate_class.append('no growth')
        elif 1.1 < gr <= 1.3:
            growth_rate_class.append('slow')
        elif 1.3 < gr < 2.0:
            growth_rate_class.append('moderate')
        else:
            growth_rate_class.append('fast')

    df = pd.DataFrame(data={"Name": sample_names, 'Growth_class': growth_rate_class,
                            'Growth_Rate': growth_rates, "No_Reads": reads_all,
                            'Initial_Bins': initial_bins, 'Used_Bins': used_bins,
                            'Fit_Err': bins_fit_error, 'Error_Codes': error_codes})
    # change filename from generic results_summary.csv to include sample name, also one level up
    filename = save_path / ".." / (sample_folder_path.name + '_results.csv')
    #filename = save_path / 'results_summary.csv'
    # save results to a single csv file
    df.to_csv(filename, sep=',', index=False, header=True, float_format='%.2f')


@click.command()
@click.option('--dir_path', '-p',
              help='Path to the sample directory. Calculate reproduction rate for all samples in the'
                   ' given dir', required=True)
def main(dir_path):
    sample_folder_path = Path(dir_path)
    experiment_name = 'output'
    evaluate_all_data(sample_folder_path=sample_folder_path, experiment_name=experiment_name,
                      save_plots=True)


if __name__ == '__main__':
    main()

    # sample_folder_path = Path('Dundee029_S94_R1.ndp.trm.s.mm.dup.mq30.calmd_subsamples')
    # sample_folder_path = Path(dir_path)
    # experiment_name = 'output'
    # evaluate_all_data(sample_folder_path=sample_folder_path, experiment_name=experiment_name,
    #                   save_plots=True)
