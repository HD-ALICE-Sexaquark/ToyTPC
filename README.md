# **LowPtTracksProject**

Generation and reconstruction of low-pT charged particles (pions, kaons, protons) into a toy simulation of the **ITS2** and **TPC** of **ALICE**, under none or low magnetic field.

## **reconstruction**

- **Description:** Modification of the example `B2a` of **GEANT4**. (_will be extended later_)

- **Requirements:**

  - GEANT4

- **Interactive usage via Graphical User Interface:**

  1. Enter the `reconstruction` dir and build with `bash build.sh`
  2. Execute with `bash run_gui.sh`
  3. Here, you have two options:
     1. Press the `/run/beamOn 1` button to inject and propagate particles from the standard input file `reconstruction/mc.csv` (input file can be changed with commands below). **IMPORTANT**: don't forget to set the magnetic field if needed (see commands below).
     2. Or, if you want to inject a single particle, you can use the commands `/ALICE/pdg_single_part` and `/ALICE/py_single_part` (see below) to set the particle's properties and to not input from file. Then, you can press the `/run/beamOn 1` button to inject and propagate it.

  - **Available menus**:

    - **Single Particle**: a single particle will be injected with just Py component, choose species "pion / kaon / proton" and momentum magnitude "Py = 25 / 50 / 75 / 100 / 150 / 200 / 225 / 250 / 300 MeV/c"

    - **Magnetic Field**: choose magnitude of z-component of magnetic field.

    - **Viewer** : viewer-related options.

    - **Style**: style-related options.

  - **Available commands**:

    - `/ALICE/input_file <filename>` : set input file (default value: `reconstruction/mc.csv`)
    - `/ALICE/traj_file <filename>` : set output file containing trajectories info (default value: `reconstruction/event<event_id>_traj.csv`)
    - `/ALICE/its_file <filename>` : set output file containing info of hits on the ITS2 (default value: `reconstruction/event<event_id>_its.csv`)
    - `/ALICE/tpc_file <filename>` : set output file containing info of hits on the TPC (default value: `reconstruction/event<event_id>_tpc.csv`)
    - `/ALICE/pdg_single_part <pdg_code>` : disables input from file, set a single particle instead with the chosen PDG code
    - `/ALICE/py_single_part <py> <units>` : y-component of the momentum of the single particle to be injected
    - `/globalField/setValue <mag_field_vector> <units>` : set magnetic field (default value `0 0 0`)
    - `/vis/ogl/export <filename>.<extension>` : export a screenshot of current view as an image file

- **Command-line/terminal usage:**

  1. Build with `bash build.sh`
  2. Enter `G4Setup_build/`
  3. Execute with:
     ```
     ./main <MC_CSV> <TRAJ_CSV> <ITS_CSV> <TPC_CSV>
     ```
     where:
     - `<MC_CSV>` : path of input file with MC particles (format described below)
     - `<TRAJ_CSV>` : path of output file containing trajectories info (format described below)
     - `<ITS_CSV>` : path of output file containing info of hits on the ITS2 (format described below)
     - `<TPC_CSV>` : path of output file containing info of hits on the TPC (format described below)

- **Input:**

  - MC-generated particles file (as in `reconstruction/mc.csv`) with the following format (which corresponds to the output of `generator/` and `injector/`):

  ```
  PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
  ```

- **Output:**

  Three files will be output. Their information can be connected (or linked) via the variable `trackID`.

  - `eventXX_its.csv` contains the ITS2 hits information, where each line correspond to a different **hit** with the following format:

  ```
  trackID, layer number, x (cm), y (cm), z (cm), Px (MeV/c), Py (MeV/c), Pz (MeV/c), deposited energy (MeV), generation process
  ```

  - `eventXX_tpc.csv` contains the TPC hits information, where each line correspond to a different **hit** with the following format:

  ```
  trackID, x (cm), y (cm), z (cm), Px (MeV/c), Py (MeV/c), Pz (MeV/c), time (ns), deposited energy (MeV), generation process
  ```

  - `eventXX_traj.csv` contains the trajectories (or true particle) information, where each line correspond to a different **particle** with the following format:

  ```
  trackID, PDG Code, initial x (cm), initial y (cm), initial z (cm), initial Px (MeV/c), initial Py (MeV/c), initial Pz (MeV/c), parentID, charge
  ```

## **injector**

- **Requirements:**

  - ROOT

- **Description:**

  ROOT macro that generates positive pions, positive kaons, protons, and their corresponding anti-particles, within an uniform distribution in transverse momentum, rapidity and azimuthal angle. It's currently set up to repeat `n_times`. Information is output in a `.csv` file.

- **Usage:**

  ```
  root -l -b -q 'injector/GenBox.C("<output_file>")'
  ```

- **Output:**

  A single CSV file called will be generated, where each line correspond to a different **particle** with the following format:

  ```
  PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
  ```

## **tools**

- `send_production.sh` -- bash script to send simulation jobs in batches, via **HTCondor**.

  - **Requirements**:

    - ROOT
    - GEANT4

  - **Usage**:

    ```
    ./send_production.sh --sim_dir <sim_dir> --run1 <run1> --run2 <run2> --serv <serv> --outsd <>
    ./send_production.sh --sim_dir <sim_dir> --rn <rn> --serv <serv> --outsd <out_sub_dir>
    where:
    <sim_dir>     = main directory (NOT output directory!)
                    (default: if ${PWD} == tools/, then sim_dir = "tools/.."
                              if not, then sim_dir = ${PWD})
    <rn>          = specific run numbers, separated by comma
                    for example: --rn 0,1,2
    <run1>        = run number of the first job (starting point of loop)
    <run2>        = run number of the last job (end of loop)
    <serv>        = select which machine to use: alice-serv<serv>
                    it can be: 10, 12 or 14
    <out_sub_dir> = subdirectory (it's a string!) within the output directory, for organizational purposes

    EXAMPLE:
    ./send_production.sh --run1 1 --run2 100 --serv 10 --outsd 001
    ```

  - **Hard-coded parameters**:

    - `OUTPUT_DIR` : top output directory
    - `N_EVENTS_PER_RUN` : number of events per run (recommended: 100)

- `merge_files.sh` -- merge every `eventX_mc.csv`, `eventX_traj.csv`, `eventX_its.csv` or `eventX_tpc.csv` within a run directory `<run_dir>` into a single `runX_mc.csv`, `runX_traj.csv`, `runX_its.csv`, or `runX_tpc.csv` (respectively) by adding the `eventID` as the first column.

  - **Usage**: `./merge_files.sh <run_dir>`

## **analysis**

(_description pending_)

- **Usage:** `root -l -b -q ParseCSVFiles(<run_n>)`, where `<run_n>` will read the respective `run<run_n>_traj.csv`, `run<run_n>_its.csv`, `run<run_n>_tpc.csv` files, and it will output a single `run<run_n>_ana.root` file.

## **generator**

(_deprecated_)

- **Requirements:** PYTHIA8 (can be easily installed with `generator/install_pythia.sh`)

- **Usage:**

  ```
  ./main_GenCollision --n <n> --config <config_file> --output <output_dir>
  ```

  where:

  - `--n <n>` : number of events to generate
  - `--config <config_file>` : input configuration file (mandatory), should look like `generator/config_pp.cmnd`
  - `--output <output_dir>` : output directory

- **Output:**

  A single file called `mc.csv` will be generated in case the number of events is 1, otherwise one `eventX_mc.csv` file will be generated per each event. In these files, each line correspond to a different **final-state particle** with the following format:

  ```
  PDG Code, Px (GeV/c), Py (GeV/c), Pz (GeV/c)
  ```
