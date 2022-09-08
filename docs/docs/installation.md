---
layout: default
title: Installation
nav_order: 10
permalink: docs/installation
---

# Installation

With version 2.0 the DO-MS functionality was extended to visulaize both DDA and DIA search engine results. To enrich the information displayed for DIA experiments, a Python preprocessing script is included. The installation instructions are therefore split into two sections, the DO-MS shiny R app and the preprocessing Python script.

To get started, download the latest version of DO-MS from the [release page](https://github.com/SlavovLab/DO-MS/releases/latest).

## DO-MS Shiny App

### Installation

This application has been tested on R >= 3.5.0, OSX >= 10.14 / Windows 7/8/10/11. Make sure you have the mos recent version of R and R Studio installed. R can be downloaded from the main [R Project page](https://www.r-project.org/) or downloaded with the [RStudio Application](https://www.rstudio.com/products/rstudio/download/). All modules are maintained for MaxQuant >= 1.6.0.16.

The application suffers from visual glitches when displayed on unsupported older browsers (such as IE9 commonly packaged with RStudio on Windows). Please use IE >= 11, Firefox, or Chrome for the best user experience.

### Running 

The easiest way to run the app is directly through RStudio, by opening the `DO-MS.Rproj` Rproject file

![]({{site.baseurl}}/assets/images/do-ms-proj.png){: width="70%" .center-image}

and clicking the "Run App" button at the top of the application, after opening the `server.R` file. We recommend checking the "Run External" option to open the application in your default browser instead of the RStudio Viewer.

![]({{site.baseurl}}/assets/images/do-ms-run.png){: width="70%" .center-image}

You can also start the application by running the `start_server.R` script.

## DIA Python Script

### Installation

1. Please make sure that [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) are installed.
Use the provided conda configuration to create the environment with all required dependencies.
```
conda env create -f module/cmd/env.yml
```

2. Activate the environment and check that the command can be run.
```
conda activate doms
python module/cmd/feature_detection.py -h
```

3. For automatic conversion of Thermo Raw files to the open mzML format ThermoRawFileParser (Hulstaert et al. 2020) is required. Download the latest release of the [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) (version v1.4.0 or newer) and note down the location of the ```ThermoRawFileParser.exe``` file. Under OSX and Linux, [Mono](https://www.mono-project.com/download/stable/). COnsult the Please make sure to use the option ```-m``` with the feature detection which will tell the script to use Mono. 

4. For feature detection Dinosaur (Teleman et al. 2016) is used. Download the latest release of the Dinosaur from [Mono](https://github.com/fickludd/dinosaur) and install Java as recommended on your platform. Please note down the location of the ```Dinosaur-xx.jar``` file.

5. Optional, create a custom script for your system.
Using the feature detection requires the correct conda environment, ThermoRawFileParser and Dinosaur. If the tool is used frequently its more convenient to contain the configuration in a script which is added to the system ```PATH```. This will register a local command which can be used everywhere on the system. 


#### Create Custom Script on Windows

1. Create a local folder for example under ```C:\Users\xxx\Documents\bin```.

2. Create a file named ```feature_detection.bat``` with the following content. Make sure that all three file paths are changed to the corresponding locations on your system.
```
@echo off
conda activate doms & ^
python C:\Users\xxx\module\cmd\feature_detection.py %* ^
--dinosaur-location "C:\Users\xxx\dinosaur-1.2\Dinosaur-1.2.0.free.jar" ^
--raw-parser-location "C:\Users\xxx\thermo_raw_file_parser_1.3.4\ThermoRawFileParser.exe" 
```
 
3. Search ```environment variables``` in the windows search and click ```Edit the system environment variables```. 

4. Click ```Environment Variables``` in the bottom right.

5. Select the variable ```Path``` in the upper panel saying ```User variables ...``` and click ```Edit```.

6. Click ```New``` and enter the location of the directory containing the ```feature_detection.bat``` script.
