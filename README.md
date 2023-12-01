**LowPtTracksProject**
======================

Generation and reconstruction of low-pT charged particles (pions, kaons, protons) from *pp* collisions at 13.6 TeV into a toy simulation of the **ITS2** and **TPC** of **ALICE**, under none or low magnetic field.

## **reconstruction**

  Modification of the example `B2a` of **GEANT4**.

* **Requirements:** GEANT4

* **Input:**

  * MC-generated particles file (as in `reconstruction/mc.csv`) with the following format (which corresponds to the output of `generator/` and `injector/`):

    ```
    PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
    ```

* **Output:**

  Three files will be output. Their information can be connected (or linked) via the variable `trackID`.

  * `eventXX_its.csv` contains the ITS2 hits information, where each line correspond to a different **hit** with the following format:

    ```
    trackID, layer number, x (cm), y (cm), z (cm), deposited energy (MeV), generation process
    ```

  * `eventXX_tpc.csv` contains the TPC hits information, where each line correspond to a different **hit** with the following format:

    ```
    trackID, x (cm), y (cm), z (cm), Px (MeV/c), Py (MeV/c), Pz (MeV/c), time (ns), deposited energy (MeV), generation process
    ```

  * `eventXX_traj.csv` contains the trajectories (or true particle) information, where each line correspond to a different **particle** with the following format:

    ```
    trackID, PDG Code, initial x (cm), initial y (cm), initial z (cm), initial Px (MeV/c), initial Py (MeV/c), initial Pz (MeV/c), parentID, charge
    ```

* **Interactive usage via Graphical User Interface:**

  First, make sure you already have both of your input files (check **Available commands** below). Then, inside the `reconstruction/` dir:

  1. Build with `bash build.sh`
  2. Execute with `bash run_gui.sh`
  3. In case you want to set different input/output files than the default ones, set the files as in **Available commands**
  4. Press the `/run/beamOn 1` button to generate a single event, propagating the injected particles

  * **Available commands**:

    * `/ALICE/signal_file <filename>` : set input signal file (default value: `reco/signal.csv`). To read only the background file, set `/ALICE/signal_file 0`
    * `/ALICE/output_file <filename>` : set output file  (default value: `reco/output_e<event_id>.csv`)
    * `/ALICE/output_file <filename>` : set output file  (default value: `reco/output_e<event_id>.csv`)
    * `/ALICE/output_file <filename>` : set output file  (default value: `reco/output_e<event_id>.csv`)
    * `/ALICE/bkg_pdg_code <pdg_code>` : PDG code of the injected background particle (default value: `-2112`)
    * `/ALICE/bkg_pdg_code <pdg_code>` : PDG code of the injected background particle (default value: `-2112`)

    (_pending_)

* **Command-line/terminal usage:**

  1. Enter `reconstruction/`
  2. Build with `bash build.sh`
  3. Enter `G4Setup_build/`
  4. Execute with:
     ```
     ./main
     ```
     (_pending_)

## **injector**

* **Requirements:** ROOT

* **Usage:**

  ```
  root -l -b -q 'injector/GenBox.C("<output_file>")'
  ```
  where `<output_file>` corresponds to the path of the output file.

* **Output:**

  A single CSV file called will be generated, where each line correspond to a different **particle** with the following format:
  ```
  PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
  ```

## **generator**

* **Requirements:** PYTHIA8 (can be easily installed with `generator/install_pythia.sh`)

* **Usage:**

  ```
  ./main_GenCollision --n <n> --config <config_file> --output <output_dir>
  ```
  where:
  * `--n <n>` : number of events to generate
  * `--config <config_file>` : input configuration file (mandatory), should look like `generator/config_pp.cmnd`
  * `--output <output_dir>` : output directory

* **Output:**

  A single file called `mc.csv` will be generated in case the number of events is 1, otherwise one `eventX_mc.csv` file will be generated per each event. In these files, each line correspond to a different **final-state particle** with the following format:
  ```
  PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
  ```

## **analysis**

* **Usage:**

## **tools**

* `send_production.sh` -- (_pending_)

* `merge_files.sh` -- merge every `eventX_mc.csv`, `eventX_traj.csv`, `eventX_its.csv` or `eventX_tpc.csv` within a run directory `<run_dir>` into a single `runX_mc.csv`, `runX_traj.csv`, `runX_its.csv`, or `runX_tpc.csv` (respectively) by adding the `eventID` as the first column.

  **Usage**: `./merge_files.sh <run_dir>`
