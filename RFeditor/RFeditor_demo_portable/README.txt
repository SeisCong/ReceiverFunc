1. Source/permission
This demo comes from the receiver functions recorded by EarthScope OIINK network, subsetted for 8 stations for demonstration purpose. Please DO NOT use this database for any other purposes without permission from the project PI (Gary Pavlis, Indiana University, pavlis@indiana.edu).

Compilation of this demo package was done by Xiaotao Yang (formerly at Indiana University from 2011-2016) at the University of Massachusetts Amherst. Contact Xiaotao Yang (xiaotaoyang@geo.umass.edu) for any questions regarding usage/issues of this demo package.

2. Structure
oiinkrfsdemo: database including metadata tables for the OIINK stations, receiver functions in wfprocess tables and related tables (e.g., evlink, sclink). the *.decon table is from deconvolution using tracedecon by Yinzhi Wang at Indiana University.
OIINKRFSDEMOwf: receiver function data (waveforms) saved in three-component format.
pf: parameter files needed by running RFeditor. RFeditor.pf is the default parameter file works for the demo. The editing parameters can be modified after the demo works. seispp_attribute_maps.pf: required to map non-standard attributes (e.g., decon, tredit) when loading the database tables.

3. USAGE:
> Pick trace using middle botton
> The letter with underline or the combination of keys shown on the right of the menu is the keyboard shortcut for that function.
> Sometimes, you need to interactively type in something in the Terminal window that opens up with this program.
> 
(1) RFeditor under review mode: edits will not be saved under this mode.
$ RFeditor oiinkrfsdemo - -rm
This command lets the program run under remove mode (rm) without specifying the output db (-). It uses the default parameter file under pf/RFeditor.pf.

Once the program launches normally, try the following items to get familar with the program while testing it.
* Manual edit modes: Edit -> Manual-Editing Modes ->
single trace (select a trace to kill using the middle botton): keyboard shortcut is 'I'
cutoff mode (traces below the selected trace on the screen will all be droped): shortcut is 'O'
* Clear killed traces in the view: View -> Exclude Killed Traces (keyboard shortcut: 'X')
* Sort traces by different attributes: Sort -> By Correlation With Ref-Trace -> Pick Reference Trace
Explore other sorting functions.
* Get editing statistics: Tools -> Statistics
* Output a trace/waveform to plain text: File -> <Save Picked-Trace To File >. You need to pick a trace for saving using the middle botton of your mouse. There is a matlab function in the Utilities folder that comes with the RFeditor codes, named readrf.m, that can read the save text file.

* Go to the next station: File -> Save & Go Next, keyboard shortcut is 'G'. Under review mode, the program doesn't save the edits you made.

(2) Test saving edits under GUI mode
$ RFeditor oiinkrfsdemo testoutdb [-pf pf/RFeditor.pf]
[*] means optional since the default PF file is RFeditor.pf. Use -pf to specify alternative PF file.

* Explore with some editing.
* File -> Save & Go Next
* Quit the program: close the window or File -> Quit Without Save
* Open the saved db using dbe: $dbe testoutdb
* Open the editing table: tredit
* Explore the table to verify it includes the right editing statistics you made.

(3) Test saving edits under GUI-off mode
GUI-off means there will be not GUI window. This feature is designed for processing large database automatically using the parameters specified in the PF file.
$ RFeditor oiinkrfsdemo testoutdbGUIoff -d testoutGUIoffdata -go [-pf pf/RFeditor.pf]
* Here -d option specifies the directory to save the edited waveform data.
* -go: means GuiOff. This option is incompatible with review mode.
* When it finishes, check the output database.
* Go back to the PF file to explore by modifying the parameters used by GUI-off auto-mode. See the block auto_edit_parameters &Arr{ }
