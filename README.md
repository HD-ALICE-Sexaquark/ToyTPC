**LowPtTracksProject**
======================

Generation and reconstruction of low-pT charged particles (pions, kaons, protons) from pp collisions at 13.6 TeV into a toy simulation of the **ITS2** and **TPC** of **ALICE**, under none or low magnetic field.

## **reconstruction**

  Modification of the example `B2a` of **GEANT4**.

* **Requirements:** GEANT4

* **Input:**

  * MC-generated particles file (as in `reconstruction/mc.csv`) with the following format:

    ```
    1,PDGCode,Px,Py,Pz
    ```

* **Output:**

  * Reconstruction file (`.csv`) with the following format:

     (_pending_)

* **Interactive usage via Graphical User Interface:**

  First, make sure you already have both of your input files (check **Available commands** below). Then, inside the `reco/` dir:

  1. Build with `bash build.sh`
  2. Execute with `bash run_gui.sh`
  3. In case you want to set different input/output files than the default ones, set the files as in **Available commands**
  4. Press the `/run/beamOn 1` button to generate a single event, propagating the injected particles

  * **Available commands**:

    * `/ALICE/signal_file <filename>` : set input signal file (default value: `reco/signal.csv`). To read only the background file, set `/ALICE/signal_file 0`
    * `/ALICE/bkg_file <filename>` : set input background file  (default value: `reco/bkg.csv`)
    * `/ALICE/output_file <filename>` : set output file  (default value: `reco/output_e<event_id>.csv`)
    * `/ALICE/bkg_pdg_code <pdg_code>` : PDG code of the injected background particle (default value: `-2112`) (**IMPORTANT**: it can be used to modify the trigger condition to keep relevant events, however, it will not update the inelastic scattering cross-section, for that, one has to restart the program with the proper code as first argument)
    * `/ALICE/mag_field <b_z>` : set the z-component of the magnetic field (in Tesla)

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

  `root -l -b -q 'injector/GenBox.C("<output_file>")'`

  where:
  * `<output_file>` : path of output file

## **generator**

  (_pending_)

## **send_production.sh**

  (_pending_)
