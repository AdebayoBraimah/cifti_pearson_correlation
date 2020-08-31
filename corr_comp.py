#!/usr/bin/env python

'''
Python script with a collection of wrapper functions for FSL and Connectome Workbench binaries to perform timeseries extraction of
CIFTI files and to compute the Pearson correlation between two ROIs.
'''

# Import packages/modules
import subprocess
import logging
import os
import numpy as np
import random
import shutil
import platform
import sys

# Import modules/packages argument parser
import argparse

# Define class(es)
class Command(object):
    '''
    Creates a command and an empty command list for UNIX command line programs/applications. Primary use and
    use-cases are intended for the subprocess module and its associated classes (i.e. Popen/call/run).
    
    Attributes (class and instance attributes):
        command (instance): Command to be performed on the command line.
        cmd_list (instance): Mutable list that can be appended to.
    
    Modules/Packages required:
        - os
        - logging
        - subprocess
    '''

    def __init__(self,command):
        '''
        Init doc-string for Command class. Initializes a command to be used on UNIX command line.
        The input argument is a command (string), and a mutable list is returned (, that can later
        be appended to).
        
        Usage:
            echo = Command("echo")
            echo.cmd_list.append("Hi!")
            echo.cmd_list.append("I have arrived!")
        
        Arguments:
            command (string): Command to be used. Note: command used must be in system path
        Returns:
            cmd_list (list): Mutable list that can be appended to.
        '''
        self.command = command
        self.cmd_list = [f"{self.command}"]
        
    def log(self,log_file="log_file.log",log_cmd=""):
        '''
        Log function for logging commands and messages to some log file.
        
        Usage:
            # Initialize the `log` function command
            log_msg = Command("log")
            
            # Specify output file and message
            log_msg.log("sub.log","test message 1")
            
            # Record message, however - no need to re-initialize `log` funcion command or log output file
            log_msg.log("test message 2")
        
        NOTE: The input `log_file` only needs to be specified once. Once specified,
            this log is written to each time this or the `run` function is invoked.
        
        Arguments:
            log_file(file): Log file to be written to. 
            log_cmd(str): Message to be written to log file
        '''
        
        # Set-up logging to file
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%d-%m-%y %H:%M:%S',
                            filename=log_file,
                            filemode='a')
        
        # Define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        
        # Add the handler to the root logger
        logging.getLogger().addHandler(console)
        
        # Define logging
        logger = logging.getLogger(__name__)
        
        # Log command/message
        logger.info(f"{log_cmd}")
        
    def run(self,log_file="",debug=False,dryrun=False,env=None,stdout="",shell=False):
        '''
        Uses python's built-in subprocess class to execute (run) a command from an input command list.
        The standard output and error can optionally be written to file.
        
        Usage:
            echo.run() # This will return tuple (returncode,log,None,None), but will echo "Hi!" to screen.
            
        NOTE: 
            - The contents of the 'stdout' output file will be empty if 'shell' is set to True.
            - Once the log file name 'log_file' has been set, that value is stored and cannot be changed.
                - This log file will continue to be appended to for each invocation of this class.

        Arguments:
            log_file(file): Output log file name.
            debug(bool): Sets logging function verbosity to DEBUG level
            dryrun(bool): Dry run -- does not run task. Command is recorded to log file.
            env(dict): Dictionary of environment variables to add to subshell.
            stdout(file): Output file to write standard output to.
            shell(bool): Use shell to execute command.
        Returns:
            p.returncode(int): Return code for command execution should the 'log_file' option be used.
            log_file(file): Output log file with appended information should the 'log_file' option be used.
            stdout(file): Standard output writtent to file should the 'stdout' option be used.
            stderr(file): Standard error writtent to file should the 'stdout' option be used.
        '''
        
        # Define logging
        logger = logging.getLogger(__name__)
        cmd = ' '.join(self.cmd_list) # Join list for logging purposes
        
        if debug:
            logger.debug(f"Running: {cmd}")
        else:
            logger.info(f"Running: {cmd}")
        
        if dryrun:
            logger.info("Performing command as dryrun")
            return 0
        
        # Define environment variables
        merged_env = os.environ
        if env:
            merged_env.update(env)
        
        # Execute/run command
        p = subprocess.Popen(self.cmd_list,shell=shell,env=merged_env,
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        # Write log files
        out,err = p.communicate()
        out = out.decode('utf-8')
        err = err.decode('utf-8')

        # Write std output/error files
        if stdout:
            stderr = os.path.splitext(stdout)[0] + ".err"
            with open(stdout,"w") as f_out:
                with open(stderr,"w") as f_err:
                    f_out.write(out)
                    f_err.write(err)
                    f_out.close(); f_err.close()
        else:
            stdout = None
            stderr = None

        if p.returncode:
            logger.error(f"command: {cmd} \n Failed with returncode {p.returncode}")

        if len(out) > 0:
            if debug:
                logger.debug(out)
            else:
                logger.info(out)

        if len(err) > 0:
            if debug:
                logger.info(err)
            else:
                logger.warning(err)
        return p.returncode,log_file,stdout,stderr

# Define functions
def threshold_cifti(nii,out,thresh,log_file="",debug=False,dryrun=False,env=None,stdout="",shell=False,verbose=False):
    '''
    Performs thresholding of NIFTI-1 files (that were converted from CIFTI-2 files),
    using FSL's cluster binary.
    
    Arguments:
        nii(file): Input NIFTI-1 image file from wb_command -cifit-convert
        out(file): Output file name for thresholded NIFTI-1 file
        thresh(float): All values below this are set to 0
        log_file(log): Log file to be written to. 
            - NOTE: if the log function has been used previously, then this argument need not be assigned.
        debug(bool): Turn on logging's diagnostic messaging
        dryrun(bool): Perform dryrun (i.e. does not generate any files)
        env(dict): Dictionary of environmental variables
        stdout(file): Standard output file to be written to.
            - NOTE: This file can only be written to if `shell` is set to False.
        shell(bool): Run the command using a shell.
        verbose(bool): Turn on verbose/diagnostic messages for UNIX command
    Returns:
        out(file): Output file name for thresholded NIFTI-1 file
    '''
    
    # Init UNIX command
    cluster = Command("cluster")
    cluster.cmd_list.append(f"--in={nii}")
    cluster.cmd_list.append(f"--thresh={thresh}")
    cluster.cmd_list.append(f"--oindex={out}")
    cluster.cmd_list.append("--no_table")
    
    if verbose:
        cluster.cmd_list.append("--verbose")
    
    # Execute command
    [exit_status,log_file,stdout,stderr] = cluster.run(log_file=log_file,
                                                       debug=debug,
                                                       dryrun=dryrun,
                                                       env=env,
                                                       stdout=stdout,
                                                       shell=shell)
    return out

def cifti_to_nifti(cii,out,thresh=0,log_file="",debug=False,dryrun=False,env=None,stdout="",shell=False,verbose=False):
    '''
    Performs conversion of input CIFTI-2 file to NIFTI-1 file via
    wb_command -cifti-convert.
    
    Arguments:
        cii(file): Input CIFTI-2 file
        out(file): Output file name for NIFTI-1 file
        thresh(float): All values below this are set to 0
        log_file(log): Log file to be written to. 
            - NOTE: if the log function has been used previously, then this argument need not be assigned.
        debug(bool): Turn on logging's diagnostic messaging
        dryrun(bool): Perform dryrun (i.e. does not generate any files)
        env(dict): Dictionary of environmental variables
        stdout(file): Standard output file to be written to.
            - NOTE: This file can only be written to if `shell` is set to False.
        shell(bool): Run the command using a shell.
        verbose(bool): Turn on verbose/diagnostic messages for UNIX command
    Returns:
        out(file): Output file name for NIFTI-1 file
    '''
    
    # Format variable
    if '.nii.gz' in out:
        out_tmp = out[:-7] + ".tmp.nii.gz"
        out_txt = out[:-7] + ".cluster.txt"
    elif '.nii' in out:
        out_tmp = out[:-4] + ".tmp.nii"
        out_txt = out[:-4] + ".cluster.txt"
    else:
        out_tmp = out + ".tmp"
        out_txt = out + ".cluster.txt"
    
    # Init UNIX command
    cii_to_nii = Command("wb_command")
    cii_to_nii.cmd_list.append("-cifti-convert")
    cii_to_nii.cmd_list.append("-to-nifti")
    cii_to_nii.cmd_list.append(f"{cii}")
    cii_to_nii.cmd_list.append(f"{out_tmp}")
    
    # Execute command
    [exit_status,log_file,stdout,stderr] = cii_to_nii.run(log_file=log_file,
                                                          debug=debug,
                                                          dryrun=dryrun,
                                                          env=env,
                                                          stdout=stdout,
                                                          shell=shell)
    
    # Threshold convert CIFTI-2 file
    if thresh:
        out = threshold_cifti(nii=out_tmp,
                              out=out,
                              thresh=thresh,
                              log_file=log_file,
                              debug=debug,
                              dryrun=dryrun,
                              env=env,
                              stdout=out_txt,
                              shell=shell,
                              verbose=verbose)
    else:
        os.rename(out_tmp,out)
    return out

def meants(nii,out,mask,log_file="",debug=False,dryrun=False,env=None,stdout="",shell=False,verbose=False):
    '''
    Performs conversion of input CIFTI-2 file to NIFTI-1 file via
    wb_command -cifti-convert.
    
    Arguments:
        nii(file): Input NIFTI-1 file
        out(file): Output file name for NIFTI-1 mean timeseries
        mask(file): Input NIFTI-1 mask file (dimensions must match input NIFTI-1 file)
        log_file(log): Log file to be written to. 
            - NOTE: if the log function has been used previously, then this argument need not be assigned.
        debug(bool): Turn on logging's diagnostic messaging
        dryrun(bool): Perform dryrun (i.e. does not generate any files)
        env(dict): Dictionary of environmental variables
        stdout(file): Standard output file to be written to.
            - NOTE: This file can only be written to if `shell` is set to False.
        shell(bool): Run the command using a shell.
        verbose(bool): Turn on verbose/diagnostic messages for UNIX command
    Returns:
        out(file): Output file name for NIFTI-1 file
    '''
    
    # Init UNIX command
    mean_ts = Command("fslmeants")
    mean_ts.cmd_list.append("-i"); mean_ts.cmd_list.append(f"{nii}")
    mean_ts.cmd_list.append("-o"); mean_ts.cmd_list.append(f"{out}")
    mean_ts.cmd_list.append("-m"); mean_ts.cmd_list.append(f"{mask}")
    
    if verbose:
        mean_ts.cmd_list.append("--verbose")
    
    # Execute command 
    [exit_status,log_file,stdout,stderr] = mean_ts.run(log_file=log_file,
                                                      debug=debug,
                                                      dryrun=dryrun,
                                                      env=env,
                                                      stdout=stdout,
                                                      shell=shell)
    return out

def cii_meants(cii,out_prefix,mask,thresh=0,log_file="",debug=False,dryrun=False,env=None,stdout="",shell=False,verbose=False):
    '''
    Computes mean timeseries for some input CIFTI-2 file with some CIFTI-2 mask file.
    
    Arguments:
        cii(file): Input CIFTI-2 file
        out_prefix(file): Output file name prefixes for intermediate files
        mask(file): Input NIFTI-1 mask file (dimensions must match input NIFTI-1 file)
        thresh(float): All values below this are set to 0
        log_file(log): Log file to be written to. 
            - NOTE: if the log function has been used previously, then this argument need not be assigned.
        debug(bool): Turn on logging's diagnostic messaging
        dryrun(bool): Perform dryrun (i.e. does not generate any files)
        env(dict): Dictionary of environmental variables
        stdout(file): Standard output file to be written to.
            - NOTE: This file can only be written to if `shell` is set to False.
        shell(bool): Run the command using a shell.
        verbose(bool): Turn on verbose/diagnostic messages for UNIX command
    Returns:
        out(file): Output file name for NIFTI-1 file
    '''
    
    # Convert CIFTI-2 to NIFTI-1
    nii_data = cifti_to_nifti(cii=cii,
                              out=out_prefix + ".nii.gz",
                              log_file=log_file,
                              thresh=0,
                              debug=debug,
                              dryrun=dryrun,
                              env=env,
                              stdout=stdout,
                              shell=shell,
                              verbose=verbose)
    mask_data = cifti_to_nifti(cii=mask,
                               out=out_prefix + ".mask.nii.gz",
                               log_file=log_file,
                               thresh=thresh,
                               debug=debug,
                               dryrun=dryrun,
                               env=env,
                               stdout=stdout,
                               shell=shell,
                               verbose=verbose)
    
    # Compute mean timeseries
    mean_ts = meants(nii=out_prefix + ".nii.gz",
                     out=out_prefix + ".mat.txt",
                     mask=out_prefix + ".mask.nii.gz",
                     log_file=log_file,
                     debug=debug,
                     dryrun=dryrun,
                     env=env,
                     stdout=stdout,
                     shell=shell,
                     verbose=verbose)
    return mean_ts

def remove_diagonal(A):
    '''
    Removes values from the main diagonal.
    
    Information/function from: 
    https://stackoverflow.com/questions/46736258/deleting-diagonal-elements-of-a-numpy-array
    
    Arguments:
        A(matrix, numpy array): Input matrix
    Returns:
        A(matrix, numpy array): Output matrix with main diagonal removed
    '''
    
    m = A.shape[0]
    strided = np.lib.stride_tricks.as_strided
    s0,s1 = A.strides
    return strided(A.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)

def pearson_corr(file1,file2,log_file=""):
    '''
    Computes Pearson correlation between two N x 1 matrices/arrays
    that are stored as text files.
    
    Arguments:
        file1(file): Input file containing N x 1 matrix
        file2(file): Input file containing N x 1 matrix
        log_file(log): Log file to be written to. 
    Returns:
        Pearson correlation coefficient(float): Pearson correlation coefficient
    '''
    
    # Log message
    log_msg = Command("log")
    log_msg.log(log_file=log_file,
                log_cmd="Computing Pearson correlation")
    
    # Load files
    A = np.loadtxt(file1)
    B = np.loadtxt(file2)
    
    # Compute Pearson correlation (assumes A & B are N x 1 matrices/arrays)
    return remove_diagonal(np.tril(np.corrcoef(A,B),k=0)).flatten(order='C')[1]

def write_to_file(out_file,text=""):
    '''
    Writes text to file.
    
    Arguments:
        out_file(file): Output file name
        text: Text to be written to file (i.e. strings, floats, ints etc.)
    Returns:
        out_file(file): Output file name
    '''
    
    # Convert text to string if not string
    if not isinstance(text,str):
        text = str(text)
        
    # Write text to file
    with open(out_file,"w") as file:
        file.write(text + "\n")
        file.close()
    return out_file

def corr_comp(cii,seed_mask,stat_mask,out_prefix,thresh=0,log_file="file.log",debug=False,dryrun=False,env=None,stdout="",shell=False,verbose=False,keep_tmp_dir=False):
    '''
    Computes mean timeseries for some input CIFTI-2 file with some CIFTI-2 mask file.
    
    Arguments:
        cii(file): Input CIFTI-2 file
        seed_mask(file): Input CIFTI-2 seed mask file (dimensions must match input CIFTI-2 file)
        stat_mask(file): Input CIFTI-2 stat mask file (dimensions must match input CIFTI-2 file)
        out_prefix(file): Output file name prefixes for intermediate files
        thresh(float): All values below this are set to 0
        log_file(log): Log file to be written to. 
        debug(bool): Turn on logging's diagnostic messaging
        dryrun(bool): Perform dryrun (i.e. does not generate any files)
        env(dict): Dictionary of environmental variables
        stdout(file): Standard output file to be written to.
            - NOTE: This file can only be written to if `shell` is set to False.
        shell(bool): Run the command using a shell.
        verbose(bool): Turn on verbose/diagnostic messages for UNIX command
    Returns:
        out(file): Output file name for NIFTI-1 file
    '''
    
    # Ascertain absolute file paths
    cii = os.path.abspath(cii)
    seed_mask = os.path.abspath(seed_mask)
    stat_mask = os.path.abspath(stat_mask)
    
    if '' in (os.path.dirname(out_prefix)):
        out_dir = os.getcwd()
    else:
        out_dir = os.path.abspath(os.path.dirname(out_prefix))
        out_name = os.path.basename(out_prefix)
    
    # Log message
    log_msg = Command("log")
    log_msg.log(log_file=log_file,
                log_cmd=f"Processing: {os.path.basename(cii)} \n \
                Seed mask: {os.path.basename(seed_mask)} \n \
                Stat mask: {os.path.basename(stat_mask)}")
    
    # Create temporary directory and filenames
    log_msg.log(log_file=log_file,
                log_cmd="Creating temporary directory")
    cwd = os.getcwd()
    n = 10000 # maximum N for random number generator
    tmp_dir = os.path.join(out_dir, 'tmp_dir_' + str(random.randint(0, n)))

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    os.chdir(tmp_dir)
    
    # Compute mean timeseries
    mean_ts_1 = cii_meants(cii=cii,
                           out_prefix="mask.seed",
                           mask=seed_mask,
                           thresh=0,
                           debug=debug,
                           dryrun=dryrun,
                           env=env,
                           stdout=stdout,
                           shell=shell,
                           verbose=verbose)
    mean_ts_2 = cii_meants(cii=cii,
                           out_prefix="mask.stat",
                           mask=stat_mask,
                           thresh=thresh,
                           debug=debug,
                           dryrun=dryrun,
                           env=env,
                           stdout=stdout,
                           shell=shell,
                           verbose=verbose)
    
    # Compute Pearson correlation coefficient
    corr_coeff = pearson_corr(mean_ts_1,mean_ts_2)
    
    # Clean-up
    if not keep_tmp_dir:
        log_msg.log(log_file=log_file,
                    log_cmd="Temporory directory and file clean-up")
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)
        # os.removedirs(tmp_dir)
    else:
        os.chdir(cwd)
        
    # Write result to file
    text_file = out_prefix + ".pear_corr.txt" 
    text_file = write_to_file(out_file=text_file,text=corr_coeff)
    
    return corr_coeff,text_file

# Write main function
def main():
    '''
    main function
    - Checks system platform
    - Checks for required external software dependencies outside of python
    - Parses arguments
    '''

    # Check system
    if platform.system().lower() == 'windows':
        print("")
        print("\tThe required software (FSL) is not installable on Windows platforms. Exiting.")
        print("")
        sys.exit(1)

    # Check for external dependencies
    if not shutil.which("fslmeants"):
        print("")
        print("\tThe required software (FSL) is not installed or on the system path. Exiting.")
        print("")
        sys.exit(1)

    if not shutil.which("wb_command"):
        print("")
        print("\tThe required software (Connectome Workbench) is not installed or on the system path. Exiting.")
        print("")
        sys.exit(1)

    # Argument parser
    parser = argparse.ArgumentParser(
        description="Computes the Pearson correlation coefficient between two masks (one being a seed mask and the other being a statistics mask). \
                    The Pearson correlation coefficient is written to an output file ending with '.pear_corr.txt'.")

    # Parse Arguments
    # Required Arguments
    reqoptions = parser.add_argument_group('Required arguments')
    reqoptions.add_argument('-i', '-in', '--input',
                            type=str,
                            dest="cii_file",
                            metavar="CIFTI.dtseries.nii",
                            required=True,
                            help="CIFTI-2 dense timeseries file (e.g. subject's fMRI timeseries mapped to some surface).")
    reqoptions.add_argument('-s', '-seed', '--seed-mask',
                            type=str,
                            dest="seed_mask",
                            metavar="CIFTI.dscalar.nii",
                            required=True,
                            help="CIFTI-2 dense scalar file (used as a seed/mask in a previous analysis).")
    reqoptions.add_argument('-a', '-stat', '--stat-mask',
                            type=str,
                            dest="stat_mask",
                            metavar="CIFTI.dscalar.nii",
                            required=True,
                            help="CIFTI-2 dense scalar file (statistics file from a previous statistical analysis, to be thresholded).")
    reqoptions.add_argument('-o', '-out', '--output-prefix',
                            type=str,
                            dest="out_prefix",
                            metavar="CIFTI.dscalar.nii",
                            required=True,
                            help="CIFTI-2 dense scalar file (statistics file from a previous statistical analysis, to be thresholded).")

    # Optional Arguments
    optoptions = parser.add_argument_group('Optional arguments')
    optoptions.add_argument('-t', '-thresh', '--thresh',
                            type=float,
                            dest="thresh",
                            metavar="FLOAT",
                            default=1.77,
                            required=False,
                            help="Cluster threshold. [default: 1.77]")
    optoptions.add_argument('-l', '-log', '--log-file',
                            type=str,
                            dest="log_file",
                            metavar="LOG",
                            default="log_file.log",
                            required=False,
                            help="Log file name. [default: 'log_file.log']")
    optoptions.add_argument('--debug',
                            dest="debug",
                            action="store_true",
                            required=False,
                            default=False,
                            help="Enables printing of diagnostic messages. [default: 'disabled']")
    optoptions.add_argument('--dry-run',
                            dest="dryrun",
                            action="store_true",
                            required=False,
                            default=False,
                            help="Preforms dry-run (e.g. no files are created). [default: 'disabled']")
    optoptions.add_argument('-v', '--verbose',
                            dest="verbose",
                            action="store_true",
                            required=False,
                            default=False,
                            help="Enables printing of verbose messages. [default: 'disabled']")
    optoptions.add_argument('--keep-tmp',
                            dest="keep_tmp",
                            action="store_true",
                            required=False,
                            default=False,
                            help="Keep temporary working directory. [default: 'disabled']")

    args = parser.parse_args()

    # Print help message in the case
    # of no arguments
    try:
        args = parser.parse_args()
    except SystemExit as err:
        if err.code == 2:
            parser.print_help()

    [corr_coeff, text_file] = corr_comp(cii=args.cii_file,
                                        seed_mask=args.seed_mask,
                                        stat_mask=args.stat_mask,
                                        out_prefix=args.out_prefix,
                                        thresh=args.thresh,
                                        log_file=args.log_file,
                                        debug=args.debug,
                                        dryrun=args.dryrun,
                                        env=None,
                                        stdout="",
                                        shell=False,
                                        verbose=args.verbose,
                                        keep_tmp_dir=args.keep_tmp)

if __name__ == "__main__":
    main()
