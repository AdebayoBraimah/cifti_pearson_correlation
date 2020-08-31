# cifti_pearson_correlation

Performs timeseries extraction of CIFTI timeseries and computes the Pearson correlation between the timeseries of two ROIs.

NOTE: `Python` environmental issues may arise. If so, try this: `export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/fslpython/envs/fslpython/lib`.

```
usage: corr_comp.py [-h] -i CIFTI.dtseries.nii -s CIFTI.dscalar.nii -a
                    CIFTI.dscalar.nii -o CIFTI.dscalar.nii [-t FLOAT] [-l LOG]
                    [--debug] [--dry-run] [-v] [--keep-tmp]

Computes the Pearson correlation coefficient between two masks (one being a
seed mask and the other being a statistics mask). The Pearson correlation
coefficient is written to an output file ending with '.pear_corr.txt'.

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i CIFTI.dtseries.nii, -in CIFTI.dtseries.nii, --input CIFTI.dtseries.nii
                        CIFTI-2 dense timeseries file (e.g. subject's fMRI
                        timeseries mapped to some surface).
  -s CIFTI.dscalar.nii, -seed CIFTI.dscalar.nii, --seed-mask CIFTI.dscalar.nii
                        CIFTI-2 dense scalar file (used as a seed/mask in a
                        previous analysis).
  -a CIFTI.dscalar.nii, -stat CIFTI.dscalar.nii, --stat-mask CIFTI.dscalar.nii
                        CIFTI-2 dense scalar file (statistics file from a
                        previous statistical analysis, to be thresholded).
  -o CIFTI.dscalar.nii, -out CIFTI.dscalar.nii, --output-prefix CIFTI.dscalar.nii
                        CIFTI-2 dense scalar file (statistics file from a
                        previous statistical analysis, to be thresholded).

Optional arguments:
  -t FLOAT, -thresh FLOAT, --thresh FLOAT
                        Cluster threshold. [default: 1.77]
  -l LOG, -log LOG, --log-file LOG
                        Log file name. [default: 'log_file.log']
  --debug               Enables printing of diagnostic messages. [default:
                        'disabled']
  --dry-run             Preforms dry-run (e.g. no files are created).
                        [default: 'disabled']
  -v, --verbose         Enables printing of verbose messages. [default:
                        'disabled']
  --keep-tmp            Keep temporary working directory. [default:
                        'disabled']
```
