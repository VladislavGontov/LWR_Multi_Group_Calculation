# LWR-Multi-Group Calculation

This program uses the multigroup method to determine the spectrum of the LWR.

At the input, the program receives homogenized concentrations of elements in the reactor core; element temperatures required for taking into account modifications; core geometry. 
The program is a console application without a user interface. To start the calculation, a command of the following type is entered on the command line:
 
	>reactor.exe input.exe output.exe
 
where reactor.exe – executable file name; input.txt – path to a file with input parameters for calculation, written according to a given template; output.txt ¬– path file to which the results of the program calculation will be written.
Below is the contents of the input.txt file for one of the studied fuel compositions.
The file contains the following necessary core configuration settings:
	POWER – reactor power, MW;
	HEIGHT – reactor core height, cm.;
	RADIUS – reactor core radius, cm.;
	TEMPERATURE_OF_MODERATOR – medium or coolant temperature, K;
	NUMBER_OF_ITERATIONS – number of calculated iterations to determine the neutron flux density spectrum;
	EFFECTIVE_DAYS_STEP – step in effective days when calculating the burnup of nuclear fuel;
	RODS_ON_OR_OF – reactivity compensation enable parameter to maintain critical state of the reactor during a the operating time;
	B-10, B-11 – the isotopic composition of the isotopes 10B and 11B of the compensating and regulating rods as a percentage;
	ELEMENTS – name of the homogenized element contained in the core;
	CONSENTRAION –  homogenized element concentration multiplied by 1024;
	TEMPERATURE – element temperature in the reactor core in the «hot» state, K;

The result of the program is the dependence of changes in homogenized concentrations and various neutron-physical parameters during the operating time of nuclear fuel.
To create the confige template use MakeTemple.exe.

This program was used for the following research: https://earchive.tpu.ru/bitstream/11683/56722/1/conference_tpu-2019-C49_p33-34.pdf
