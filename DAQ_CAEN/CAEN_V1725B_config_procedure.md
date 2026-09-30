# V1725B New Card Setup Runbook

Consolidated, tested procedure for bringing up a new V1725B digitizer on the existing V1718 (VME-USB) chain. Every command here was verified working during the original bring-up session for boards 426 / 488 / 421.

Replace `<VME_BASE_ADDR>` throughout with this new board's actual VME base address (see Step 1).

---

## Step 0 — Before powering anything on

- [ ] Confirm the crate is **off**.
- [ ] Seat the new V1725B in a free slot.
- [ ] **Physically check the on-board SW2 DIP switch is set to INT**, not EXT. This is the single most common failure mode encountered tonight — an EXT-set switch causes a "PLL not locked / Board Failure" error that no firmware or software step can fix, because it's a hardware clock-routing switch, not a firmware setting.
- [ ] Confirm no cable is connected to the front-panel CLK-IN connector (unless you specifically intend external-clock synchronization for this board).

## Step 1 — Set the VME base address (rotary switches)

- Read the board's front-panel rotary switches — four hex digits (0–F each).
- Append `0000` to get the full 32-bit VME base address. Example: dials reading `5A30` → base address `0x5A300000`.
- **Confirm this address doesn't collide with any board already in the crate** (426 = `0x04260000`, 488 = `0x32100000`, 421 = `0x32120000`).
- Record the board's serial number (visible on a label) alongside this address — you'll want this mapping permanently.

## Step 2 — Power on and confirm basic connectivity

Power the crate on. Then:

```bash
caen-toolbox dig1 check -c USB -l 0 -b <VME_BASE_ADDR>
```

Expected output: `Model: V1725B`, a serial number, current firmware info, and an **Unlock code** — record that unlock code now, against this board's serial number. It's what lets you restore DPP firmware later if this board ever needs to go back to DPP mode.

If this command fails to open the board at all (not a PLL error, a connection failure), stop here — that's an addressing, seating, or USB-link problem, not anything covered by the rest of this runbook.

## Step 3 — (Only if you need Waveform Recording / WaveDump support)

Skip this step if the board should stay on its current firmware (e.g. DPP-DAW for CoMPASS use).

### 3a. Flash the D-WAVE firmware

```bash
caen-toolbox dig1 upgrade -c USB -l 0 -b <VME_BASE_ADDR> ~/CAEN/firmware/x725_rev4.29_0.09.cfa
```

Confirm the tool logs every page programmed *and* verified (not just "success") — watch for something like `Verifying page: N/N` with N matching total page count.

You'll see: `Flash programmed successfully. Please power cycle your device now.`

**This message is not optional** — the new firmware sits in flash but the FPGA keeps running the old bitstream until a genuine power cycle.

### 3b. Power cycle (physically off, wait \~5s, back on)

### 3c. Verify the firmware switch took effect

```bash
caen-toolbox dig1 check -c USB -l 0 -b <VME_BASE_ADDR>
```

Confirm the AMC firmware line now shows `(Waveform Rec.)` instead of a DPP variant.

## Step 4 — Reload the PLL configuration

Even a fresh D-WAVE flash needs its PLL config reloaded — this is a separate flash region from the application firmware.

```bash
caen-toolbox dig1 pll -c USB -l 0 -b <VME_BASE_ADDR> -R /opt/caen-toolbox/_internal/_common/pll/v1725_vcxo500_ref50_pll_out0.rbf
```

This assumes: internal 50 MHz reference, no clock output propagated (`ref50`, `out0`) — correct for a standalone board with no CLK-IN cable and no daisy-chain sync to other boards yet.

Again: `Flash programmed successfully. Please power cycle your device now.`

**Power cycle again, genuinely.**

## Step 5 — Confirm PLL lock

Two independent checks — do both:

**Physical**: look at the board's front-panel LEDs. **PLL LOCK should be on. CLK IN should be off** (unless you intentionally wired external sync). If CLK IN is lit here despite no cable connected and despite loading `ref50`, go back to Step 0 — the SW2 switch is very likely still set to EXT.

**Register-level**:

```bash
caen-toolbox dig1 dump -c USB -l 0 -b <VME_BASE_ADDR> dump_check.csv
grep "0x8104" dump_check.csv
```

## Step 6 — First WaveDump test

```bash
cd ~/CAEN/wavedump-configs
sed 's/0x[0-9A-Fa-f]\{8\}/<VME_BASE_ADDR>/' pedestal_slot3.txt > pedestal_newcard.txt
wavedump pedestal_newcard.txt
```

If you see `Board error detected: PLL not locked. Board Failure` here despite Step 5 looking clean, don't re-run PLL/firmware steps blindly — re-check the physical LEDs first; this was the actual root cause both times it came up tonight.

Once it launches cleanly:

- Press **SPACE** to see the full key list for this WaveDump build before doing anything else (key bindings can differ from what you expect).
- `s` to start/stop acquisition.
- Fire software triggers to build up pedestal statistics (aim for 100+ events per channel for good statistics — 18 was a bit thin for the demo run).
- `q` to quit.

**Known quirk**: `OUTPUT_PATH` is not a valid config directive in this WaveDump build — it's silently rejected. Output files land in whatever directory you ran `wavedump` from. Plan your working directory accordingly rather than relying on the config file for this.

## Step 7 — Pedestal sanity check

Once you have `wave0.txt` through `wave15.txt`:

```bash
python3 analyze_pedestal.py
```

(Adjust the script's channel count / file paths if this board has a different channel count, or if you're storing this board's output somewhere other than the default location.)

Look for:

- All channels clustered around a consistent pedestal mean (not necessarily identical across channels, but no wild outliers).
- Small, consistent standard deviation (noise) per channel — the reference boards showed roughly 1.3–1.7 ADC counts.
- **Any channel reading flat zero across every sample** — this is not noise, it's a disconnected or dead channel (this is exactly what channel 13 showed on the first board tested).

---

## Quick reference — commands used throughout

| Purpose | Command |
| --- | --- |
| Check firmware/serial/unlock code | `caen-toolbox dig1 check -c USB -l 0 -b <ADDR>` |
| Flash new firmware | `caen-toolbox dig1 upgrade -c USB -l 0 -b <ADDR> <file>.cfa` |
| Reload PLL config | `caen-toolbox dig1 pll -c USB -l 0 -b <ADDR> -R <file>.rbf` |
| Dump raw registers | `caen-toolbox dig1 dump -c USB -l 0 -b <ADDR> <output>.csv` |
| Run WaveDump | `wavedump <config>.txt` |

## Board log (update this table as you bring up each new card)

| Slot | VME base address | Serial | Firmware | Unlock code | Notes |
| --- | --- | --- | --- | --- | --- |
| 3 | 0x04260000 | 426 | D-WAVE (was DPP-DAW) | 00048D74A9BF0400 | ch13 dead |
| 5 | 0x32100000 | 488 | D-WAVE (was DPP-DAW) | 00048D2C08BA0C00 |  |
| 7 | 0x32120000 | 421 | D-WAVE (was DPP-DAW) | 0004813E96AA0000 | needed PLL reload, CLK-IN was off |
|  |  |  |  |  |  |