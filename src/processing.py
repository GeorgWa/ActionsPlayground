#!/usr/bin/env python3

import collections
import queue
import sys, getopt
import argparse
import os
import pathlib
import shutil
import subprocess
import pandas as pd
from datetime import datetime
import xml.etree.ElementTree as ET
import pymzml
import numpy as np
from tqdm import tqdm
from scipy.sparse import coo_matrix

import matplotlib.pyplot as plt
from scipy import stats

import matplotlib.pyplot as plt

import matplotlib.colors as colors
import matplotlib.cm as cmx
import multiprocessing

from diannbase.processing import MatchReport

def generate_parser():
        
    # instantiate the parser
    parser = argparse.ArgumentParser(description=('Command line tool for feature detection in shotgun MS experiments. Can be used together '
    'with DIA-NN to provide additional information on the peptide like features identified in the MS1 spectra.'))
    
    # Required command line arguments
    parser.add_argument('report', type=pathlib.Path, default=os.getcwd(), help="Path pointing to report.tsv output from DIA-NN.")
    
    parser.add_argument("--raw-parser-location", type=pathlib.Path, required=True, help="Path pointing to the ThermoRawFileParser executeable.")
    
    # optional command line arguments
    parser.add_argument("--dinosaur-location", help="Path pointing to the dinosaur jar executeable.")

    parser.add_argument("-m","--mono", default=False, action='store_true',  help="Use mono for ThermoRawFileParser under Linux and OSX.")

    parser.add_argument("-d","--delete", default=False, action='store_true',  help="Delete generated mzML and copied raw files after successfull feature generation.")
    
    parser.add_argument("-v","--verbose", default=False, action='store_true',  help="Show verbose output.")
    
    parser.add_argument("-t","--temporary-folder", type=pathlib.Path, help="Input Raw files will be temporarilly copied to this folder. Required for use with Google drive.")

    parser.add_argument("--no-feature-detection", default=False, action='store_true', help="All steps are performed as usual but Dinosaur feature detection is skipped. No features.tsv file will be generated.")

    parser.add_argument("--no-fill-times", default=False, action='store_true', help="All steps are performed as usual but fill times are not extracted. No fill_times.tsv file will be generated.")

    parser.add_argument("--no-tic", default=False, action='store_true', help="All steps are performed as usual but binned TIC is not extracted. No tic.tsv file will be generated.")

    parser.add_argument("--no-sn", default=False, action='store_true', help=('Signal to Noise ratio is not estimated for precursors'))

    parser.add_argument("--no-mzml-generation", default=False, action='store_true', help=('Raw files are not converted to .mzML. '
    'Nevertheless, mzML files are expected in their theoretical output location and loaded. Should be only be carefully used for repeated calulcations or debugging'))       

    parser.add_argument("--mz-bin-size", default=10, help="Bin size over the mz dimension for TIC binning.")

    parser.add_argument("--tic-dense", default=False, action='store_true', help="Return dense TIC matrices. If this option is selected dense TIC matrices are exported for every dataset seperately. The default behavior is a single sparse output file.")

    parser.add_argument("--resolution", default=70000, help="Set the resolution used for estimating counts from S/N data")

    parser.add_argument("-p", "--processes", default=10, help="Number of Processes")

    parser.add_argument("--isotopes-sn", default=False, action='store_true', help="Use all isototopes from the same scan as the highest intensity datapoint for estimating the SN and copy number.")

    return parser

class FeatureDetection():

    def __init__(self):
        pass

    def get_timestamp(self):
        # datetime object containing current date and time
        now = datetime.now()

        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  
        return "[" + dt_string + "] "

    def log(self, msg):
        print(self.get_timestamp() + msg)

    

    def __call__(self):

        parser = generate_parser()
        self.args = parser.parse_args()       

        self.output_folder = pathlib.Path(self.args.report).parent.resolve()
        self.report_tsv = pd.read_csv(self.args.report, sep='\t')
        self.experiment_files = list(set(self.report_tsv['File.Name']))
        self.experiment_files = [pathlib.Path(file) for file in self.experiment_files]
        self.experiment_names = [path.stem for path in self.experiment_files]

        self.log('The following experiments were found:')
        for exp in self.experiment_names:
            self.log(exp)

        self.mzml_files = [os.path.join(self.output_folder,'.'.join([name,'mzML'])) for name in self.experiment_names]
        self.mzml_profile_files = [os.path.join(self.output_folder,'.'.join([name,'profile','mzML'])) for name in self.experiment_names]
        
        if self.args.temporary_folder is not None:
            #temporary folder is defined, copy files to temporary folder and change input path
            
            for i, experiment in enumerate(self.experiment_files):
        
                input_raw = experiment
                output_raw = os.path.join(self.args.temporary_folder, experiment.name)
                shutil.copy(input_raw, output_raw)
                self.log(f"{experiment.name} copied to {self.args.temporary_folder}")
                
                self.experiment_files[i] = output_raw
        
        # Check if --no-mzml-generation has been set
        if not self.args.no_mzml_generation:
            self.mzml_generation()
            pass
                
        # Check if --no-feature-detection has been set
        if not self.args.no_feature_detection:
            self.feature_detection()
            pass

        # Check if --no-fill-times has been set
        if not self.args.no_fill_times:
            self.fill_times()
            pass

        # Check if --no-tic has been set
        if not self.args.no_tic:
            self.tic()
            pass

        # Check if --no-sn has been set
        if not self.args.no_sn:
            self.sn()           
        
        # delete temporary files if specified
        if self.args.delete:
            for input_mzml in self.mzml_files:
                self.log('delete' + str(input_mzml))
                if os.path.isfile(input_mzml):
                    os.remove(input_mzml)

            for input_mzml in self.mzml_profile_files:
                self.log('delete' + str(input_mzml))
                if os.path.isfile(input_mzml):
                    os.remove(input_mzml)
                    
            if self.args.temporary_folder is not None:
                for experiment in self.experiment_files:
                    self.log(str(experiment))
                    if os.path.isfile(experiment):
                        os.remove(experiment)
    def sn(self):

        
        matcher = MatchReport(self.mzml_files, self.args.report, add_isotopes=self.args.isotopes_sn)
        df = matcher()
        df.to_csv(os.path.join(self.output_folder,'sn.tsv') ,sep='\t', index=False)


    def mzml_generation(self):
        """
        The function converts the provided Thermo raw files to the open mzML format.
        Output files with centroided and profile mode spectra are created ad specified in self.mzml_files and self.mzml_profile_files
        """

        labels = ["-N"]*len(self.experiment_files)+["-p"]*len(self.experiment_files)
        mzmls = self.mzml_files + self.mzml_profile_files
        raw = self.experiment_files + self.experiment_files
        
        queue = list(zip(raw, mzmls, labels))

        with multiprocessing.Pool(processes = self.args.processes) as pool:
            pool.map(self.mzml_job, queue)

    def mzml_job(self, args):

        input_raw, mzml, label = args

        if self.args.mono:
            process = subprocess.Popen(['mono',str(self.args.raw_parser_location),label,'-i',str(input_raw),'-b',str(mzml)],shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            process = subprocess.Popen(['ThermoRawFileParser',label,'-i',str(input_raw),'-b',str(mzml)], executable=str(self.args.raw_parser_location), stdout=subprocess.PIPE, stderr=subprocess.PIPE)         

        process.wait()
        out, err = process.communicate()

        self.log(out.decode())
        self.log(err.decode())

    def feature_detection(self):

        # Dinosaur output columns are defined to make reordering of the Run column easier
        # Needs to be updated if any new columns are added 
        columns = ['mz','mostAbundantMz','charge','rtStart','rtApex','rtEnd','fwhm','nIsotopes','nScans','averagineCorr','mass','massCalib','intensityApex','intensitySum']

        # check if --dinosaur-location has been provided and is valid file path.
        if (self.args.dinosaur_location == None) or (not os.path.isfile(self.args.dinosaur_location)):
            self.log(f"Dinosaur location not valid: {self.args.dinosaur_location}")
            return
        process_list = []

        #perform feature detection based on mzml file
        for input_mzml in self.mzml_profile_files:
            process_list.append(subprocess.Popen(['java',
                                        '-jar',
                                        str(self.args.dinosaur_location),
                                        '--concurrency=10',
                                        
                                        str(input_mzml)], stdout=subprocess.PIPE, stderr=subprocess.PIPE))

        for process in process_list:
            process.wait()
            out, err = process.communicate()

            self.log(out.decode())
            self.log(err.decode())
        
        
        # convert outputs to dataframes
        dfs = []
        for i, experiment_name in enumerate(self.experiment_names):
            feature_file = '.'.join([experiment_name,'profile','features','tsv'])
            input_feature = os.path.join(self.output_folder, feature_file)
            
            df = pd.read_csv(input_feature, sep='\t', header=0)
            
            # create new Run column for DO-MS compatability
            df['Run'] = experiment_name
            dfs.append(df)

        # concatenate outputs to a single dataset    
        df = pd.concat(dfs, ignore_index=True)
        # reorder columns
        df = df.reindex(columns=['Run']+columns)
        
        output_features = os.path.join(self.output_folder, "features.tsv")
        df.to_csv(output_features,index=False, header=True, sep='\t')
        
        
        # delete dinosaur qc folder 
        qc_path = os.path.join(os.getcwd(),'QC')
        if os.path.isdir(qc_path):
            shutil.rmtree(qc_path)

    def fill_times(self):
        # List which will be used to collect datapoints
        df_to_be = []

        
        for input_mzml, name in zip(self.mzml_files, self.experiment_names):
            self.log(f'Collecting fill times for {name}')
            
            try:
                run = pymzml.run.Reader(input_mzml)
                for spectrum in run:

                    ms_level = spectrum.ms_level

                    cv_params = spectrum.get_element_by_path(['scanList', 'scan', 'cvParam'])
                    for el in cv_params:
                        key_val = el.attrib
                        
                        # parse 
                        if key_val['name'] == 'scan start time':
                            scan_start_time = key_val['value']

                        if key_val['name'] == 'ion injection time':
                            injection_time = key_val['value']


                    df_to_be.append([name, ms_level, scan_start_time, injection_time])
            except:
                print(f"{name} failed")

        fill_times_df = pd.DataFrame(df_to_be, columns =['Run', 'Ms.Level', 'RT.Start', 'Fill.Time'])
        output_fill_times = os.path.join(self.output_folder, "fill_times.tsv")
        fill_times_df.to_csv(output_fill_times,index=False, header=True, sep='\t')

    def tic(self):
        collect_dfs = []

        for input_mzml, name in zip(self.mzml_files, self.experiment_names):

            # get tic for single run
            result = self.tic_single(input_mzml, name)

            # only sparse mode returns df
            if isinstance(result, pd.DataFrame):
                collect_dfs.append(result)
        
        # concat TIC dfs from multiple runs and export them as tsv
        if len(collect_dfs) > 0:
            out_df = pd.concat(collect_dfs)
            output_tic = os.path.join(self.output_folder, "tic.tsv")
            out_df.to_csv(output_tic,index=False, header=True, sep='\t')

    def tic_single(self, input_mzml, name):

        run = pymzml.run.Reader(input_mzml)

        ms1_id = -1

        # tic values are collected spectrum wise in a sparse format
        # contains lists of [[ms1_id, current_bin, current_tic], ... ]
        spectra_batch = []

        # contains retention times
        rt_label = []


        for spectrum in run:
            ms_level = spectrum.ms_level
            if ms_level == 1:        
                ms1_id += 1

                # get total ion current meta data
                cv_params = spectrum.get_element_by_path(['cvParam'])
                for el in cv_params:
                    key_val = el.attrib

                    if key_val['name'] == "total ion current":
                        total_ion_current = float(key_val['value'])

                # get total ion current meta data
                cv_params = spectrum.get_element_by_path(['scanList', 'scan', 'cvParam'])
                for el in cv_params:
                    key_val = el.attrib

                    if key_val['name'] == 'scan start time':
                        scan_start_time = float(key_val['value'])    

                rt_label.append(scan_start_time)

                data = np.array(spectrum.peaks("raw"))           
                intensity = data[:,1]
                mz = data[:,0]

                bins = self.calc_sparse_binned_tic(mz, intensity, total_ion_current, ms1_id)
                spectra_batch.append(bins) 

        # list is converted to array and sublists are concatenated
        sparse_arr = np.concatenate(spectra_batch)

        # sparse tic matrix is converted to dense tic matrix
        rt_col = sparse_arr[:,0].astype(int)
        mz_col = sparse_arr[:,1].astype(int)
        tic = sparse_arr[:,2]

        # check if sparse output has been selected
        # Easier to handle in ggplot2
        if not self.args.tic_dense:
            raw_file = [name] * len(tic)
            data = {'Raw.File': raw_file, 
            'RT': rt_col,
            'MZ': mz_col,
            'TIC': tic}
            return pd.DataFrame.from_dict(data)


        sparse_tic = coo_matrix((tic, (rt_col, mz_col)))
        dense_tic = sparse_tic.todense()

        # The following block is needed to generate labels for the mz and RT bins
        # I am surprised how complicated I made it. I'm sure it was late and there is an obvious solution
        min_mz_untransformed = np.min(mz_col)
        min_mz = min_mz_untransformed*self.args.mz_bin_size

        rt_label = np.floor(rt_label)
        unique_label = np.unique(rt_label)

        # perform rt binning. Easy because of dense matrix format
        time_binned_mat = []
        for minute in unique_label:
            idx = np.where(rt_label == minute)
            time_binned = np.sum(dense_tic[idx,:],axis=1)
            time_binned_mat.append(time_binned)
            
        time_binned_mat = np.concatenate(time_binned_mat)

        rt_len, mz_len = time_binned_mat.shape
        rt_index = np.arange(rt_len)
        mz_index = np.arange(mz_len)*self.args.mz_bin_size
        first_index = np.where(min_mz==mz_index)[0][0]
        mz_index = mz_index[first_index:]
        
        # start mz matrix from first non empty mz bin observed
        time_binned_mat = time_binned_mat[:,first_index:]

        # Assemble pandas dataframe for tsv output
        rt_df = pd.DataFrame(data=rt_index, columns=["RT"])
        tic_df = pd.DataFrame(data=time_binned_mat, columns=mz_index)
        df = pd.concat([rt_df, tic_df], axis=1)

        output_tic = os.path.join(self.output_folder, f"tic_{name}.tsv")
        df.to_csv(output_tic,index=False, header=True, sep='\t')

        return None
        

    def calc_sparse_binned_tic(self, mz, intensity, total_ion_current, ms1_id):
        # Resolution, integer, number of mz / bin
        # the integral of the provided spectra is approximated by using the rectangle rule
        # This allows efficient binning by rounding without iterative calculation of integrals
        mz_diff = np.diff(mz)
        mz_diff = np.pad(mz_diff, (0, 1), 'constant')
        
        # The true integral is calculated as the dot product
        integral = np.dot(mz_diff,intensity)

        # rescaled to resample TIC provided in meta information
        scaling_factor = total_ion_current / integral
        scaled_intensity = intensity * scaling_factor
        binned_intensity = scaled_intensity * mz_diff
        
        binned_mz = np.round(mz/self.args.mz_bin_size)
        binned_mz = binned_mz.astype(int)
        
        sparse_binned_tic = []

        current_bin = binned_mz[0]
        current_tic = 0
        
        for mz_bin, int_bin in zip(binned_mz, binned_intensity):
            if mz_bin == current_bin:
                current_tic += int_bin

            else:
                # close last bin, check if bin has tic
                if current_tic > 0:
                    sparse_binned_tic.append((ms1_id, current_bin, current_tic))

                # open new bin
                current_bin = mz_bin
                current_tic = int_bin
                
        return np.array(sparse_binned_tic)

if __name__ == "__main__":
    feature_detection = FeatureDetection()
    feature_detection()
