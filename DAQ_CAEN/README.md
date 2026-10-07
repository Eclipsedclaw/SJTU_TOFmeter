# V1725B New Card Setup Runbook

Consolidated, tested procedure for bringing up a new V1725B digitizer on the existing V1718 (VME-USB) chain. Every command here was verified working during the original bring-up session for boards 426 / 488 / 421.

Replace `<VME_BASE_ADDR>` throughout with this new board's actual VME base address (see Step 1).

---

## Step 0 — Before powering anything on

- [ ] Confirm the crate is **off**.
- [ ] Seat the new V1725B in a free slot. Do not put over 4 V1725B inside the VME8002 unit at the same time, it will flow overcurrent.
- [ ] **Physically check the on-board SW2 DIP switch is set to INT**, not EXT. An EXT-set switch causes a "PLL not locked / Board Failure" error that no firmware or software step can fix, because it's a hardware clock-routing switch, not a firmware setting. If you see CLK IN LED lights ON, that means the card is in EXT mode and expect an external clock.
<img width="312" height="430" alt="CAEN_V1725_switch" src="https://github.com/user-attachments/assets/c3c445d0-5335-41e2-97be-61323775ab51" />

- [ ] Confirm no cable is connected to the front-panel CLK-IN connector (unless you specifically intend external-clock synchronization for this board).

## Step 1 — Set the VME base address (rotary switches)

- Read the board's rotary switches located on the back on the PCB — four hex digits (0–F each).
- Append `0000` to get the full 32-bit VME base address. Example: dials reading `5A30` → base address `0x5A300000`.
- **Confirm this address doesn't collide with any board already in the crate** (426 = `0x04260000`, 488 = `0x32100000`, 421 = `0x32120000`, these 3 were used for test).

## Step 2 — Power on and confirm basic connectivity

Power the crate on. Then:

```bash
caen-toolbox dig1 check -c USB -l 0 -b <VME_BASE_ADDR>
```

Expected output: `Model: V1725B`, a serial number, current firmware info, and an Unlock code.

If this command fails to open the board at all (not a PLL error, a connection failure), stop here — that's an addressing, seating, or USB-link problem, not anything covered by the rest of this runbook.

## Step 3 — (Only if you need Waveform Recording / WaveDump support)

Skip this step and beyond if the board should stay on its current firmware (e.g. DPP-DAW for CoMPASS use).

### 3a. Flash the D-WAVE firmware
Download the firmware from CAEN website if it is not already on your local: [https://caen.it/products/v1725-v1725s/](https://caen.it/products/v1725-v1725s/). This particular procedure introduces x725 Waveform Recording Firmware x725_rev4.29_0.09.cfa (V1725SB uses a different version of wavedump firmware)

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
Here is an example script `pedestal_slot3.txt` for taking data with slot3 CAEN card which has base address `0x04260000`. Change its variable respectively for your own card's base address.

```bash
[COMMON]
OPEN USB 0 0x04260000

RECORD_LENGTH  200
DECIMATION_FACTOR  1
POST_TRIGGER  80
PULSE_POLARITY  POSITIVE
EXTERNAL_TRIGGER   DISABLED
FPIO_LEVEL  NIM
OUTPUT_FILE_FORMAT  ASCII
OUTPUT_FILE_HEADER  NO
TEST_PATTERN   NO

ENABLE_INPUT          YES
BASELINE_LEVEL        50
TRIGGER_THRESHOLD     100
CHANNEL_TRIGGER       DISABLED

[0]
ENABLE_INPUT           YES
[1]
ENABLE_INPUT           YES
[2]
ENABLE_INPUT           YES
[3]
ENABLE_INPUT           YES
[4]
ENABLE_INPUT           YES
[5]
ENABLE_INPUT           YES
[6]
ENABLE_INPUT           YES
[7]
ENABLE_INPUT           YES
[8]
ENABLE_INPUT           YES
[9]
ENABLE_INPUT           YES
[10]
ENABLE_INPUT           YES
[11]
ENABLE_INPUT           YES
[12]
ENABLE_INPUT           YES
[13]
ENABLE_INPUT           YES
[14]
ENABLE_INPUT           YES
[15]
ENABLE_INPUT           YES
```


```bash
cd ~/CAEN/wavedump-configs
wavedump pedestal_slot3.txt
```

If you see `Board error detected: PLL not locked. Board Failure` here despite Step 5 looking clean, don't re-run PLL/firmware steps blindly — re-check the physical LEDs first; this was the actual root cause both times it came up tonight.

Once it launches cleanly:

- Press **SPACE** to see the full key list for this WaveDump build before doing anything else (key bindings can differ from what you expect).
- `s` to start/stop acquisition.
- Fire software triggers to build up pedestal statistics (aim for 100+ events per channel for good statistics — 18 was a bit thin for the demo run).
- `q` to quit.


## Step 7 — Pedestal sanity check

Once you have `wave0.txt` through `wave15.txt`:

```bash
python3 analyze_pedestal.py <data_folder>
```

Look for:

- All channels clustered around a consistent pedestal mean (not necessarily identical across channels, but no wild outliers).
- Small, consistent standard deviation (noise) per channel — the reference boards showed roughly 1.3–1.7 ADC counts.
- **Any channel reading flat zero across every sample** — this is not noise, it's a disconnected or dead channel.

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
| 7 | 0x32120000 | 421 | D-WAVE (was DPP-DAW) | 0004813E96AA0000 |  |
|  |  |  |  |  |  |
