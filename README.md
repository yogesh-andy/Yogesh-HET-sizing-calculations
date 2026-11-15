# Cathode sizing for a small argon Hall effect thruster

Calculation files for the MSc thesis "Testing and Thermionic Cathode
Development for an Argon Hall Effect Thruster" (Cranfield University, 2026,
Y. Andiyappan, supervisor Dr E. Grustan-Gutierrez).

## Files

| File | What it does | Thesis chapter |
|---|---|---|
| `cathode_sizing.m` | Sizing tool: all the calculations below in one MATLAB script | 4, 5, 6 |
| `Bare_W_Filament_Calculator.xlsx` | Bare tungsten filament: temperature, heater power, current and voltage | 4, 6 |
| `thoriated tungsten sizing and calculation.xlsx` | Phase-2 thoriated-tungsten filament sizing (sheet `Phase2_Calc`) | 4, 6 |
| `Langmuir_Sheet.xlsx` | Collector current and emission plateaus of the cathode test rig | 6 |
| `Hollow_Cathode_Design_Sizing_Tool_UPDATED.xlsx` | Hollow cathode sizing, validation cases and argon design | 5 |
| `components_list.xlsx` | Full list of bought components | Appendix C |

## What the tool is for

`cathode_sizing.m` answers the questions a student asks before building or
testing a cathode for this thruster:

- **Filament cathode.** For a wire of a given material, diameter and length,
  what heater current, voltage and power give the emission current the
  thruster needs? How hot does the wire get? Does the heater supply have
  enough current and voltage, and is one feedthrough pin enough?
- **Cathode test rig.** What collector voltage is needed to draw that
  current, and how much current should each heater setting give? These are
  the numbers to expect during the test.
- **Hollow cathode.** For a given thruster power, how long must the LaB6
  insert be, how hot must it run, what pressure is reached inside it, and
  what heater is needed?

It is useful when something changes: a different wire, a higher discharge
current, a new power supply, a different thruster power, or a new emitter
material. Change the input and the script gives the new operating point
and flags any limit that is exceeded.

## How to use it, step by step

1. **Get the file.** Download `cathode_sizing.m` and put it in a folder.
2. **Open MATLAB** (any recent version, no toolboxes needed). Set that folder
   as the current folder (the folder bar at the top, or `cd` in the command
   window).
3. **Open the script** (`edit cathode_sizing`) and go to **Section 1,
   INPUTS**. This is the only part to change. Every input has its unit and a
   short note beside it.
4. **Change what you need.** Examples:
   - More emission current: change `fil.I_em = 0.50;` to the new value, in A.
   - A thicker or longer wire: change `fil.d` (m) or `fil.L` (m).
   - A new emitter: add a row to `fil.cases` with a name, the material
     (`'W-ThO2'` or `'W'`) and the emittance at 2000 K and 3000 K.
   - A new heater supply or feedthrough: change `lim.heater_I`,
     `lim.heater_V` or `lim.pin_I`.
   - Different heater steps for the test rig: change `rig.heater_steps`.
   - A different thruster power: change `hc.P_d` and `hc.V_d`.
5. **Save and run.** Type `cathode_sizing` in the command window and press
   Enter. It takes a few seconds.
6. **Read the results.** They are printed in four blocks:
   - **A. Filament cathodes.** For each emitter, three columns: a uniform
     wire (quick hand estimate), with heat lost into the joints (the design
     value), and with the design margin. Lines starting `OK` passed a check.
     Lines starting `CHECK` need attention, for example "2 pins per lead
     needed".
   - **B. Cathode test rig.** The collector voltage for the required current,
     and a table of the expected emission at each heater current.
   - **C. Hollow cathode.** Insert length, temperatures, pressure and heater.
   - **D. Validation cases.** Three published cathodes, recalculated to show
     the method works. These are fixed and are not inputs.
7. **Look at the figures.** Two PNG files are saved in the folder:
   `collector_iv.png` (collector current against voltage for each heater
   step) and `emission_vs_heater.png` (emission current against heater
   current for each emitter).
8. **Keep the numbers (optional).** Type `r = cathode_sizing;` instead. All
   results are then in the structure `r`, for example
   `r.filaments{1}.I_h` is the heater current of the first emitter.

To skip the figures or the validation cases, set `opt.figures` or
`opt.validation` to `false` in Section 1.

## What the model assumes

- Emission follows the Richardson-Dushman law with published constants
  (Wolf 1995 for tungsten and thoriated tungsten, Lafferty 1951 for LaB6).
- The filament is solved along its length: Joule heating, radiation,
  electron cooling and conduction into joints held at `fil.T_joint`. The
  joint temperature is an estimate (400-800 K changes the heater current by
  under 1 %).
- Thoriated tungsten is assumed fully activated. A new wire may need a
  higher current until activation is complete.
- The collector current is limited by space charge in a coaxial gap.
- The results are a starting point for a test, not a substitute for one.
  Always raise the heater current slowly, at constant current, and watch
  the emission and pressure.
