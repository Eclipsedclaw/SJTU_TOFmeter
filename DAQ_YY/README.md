# Set up and ping a DAQ card on dbar

Follow this for each digitizer card that connects to dbar by fiber, one card at a time: connect it, give its server port an address, ping it, make the setting permanent, and record the result. Allow about ten minutes per card.

Server

dbar, a Dell PowerEdge R740 running Debian 12

Fiber ports

`ens1f0` to `ens1f3` (network card in riser slot 1) and `ens3f0` to `ens3f3` (riser slot 3)

Lab network

`eno4`, 192.168.35.100. Never change it.

Card settings

One small file per card in `/etc/network/interfaces.d/`

Already done

Card 1: 192.168.69.19 on `ens3f3`, server side 192.168.69.1/24

## Before you start

- Open a terminal on dbar and run this once, typing your password when asked:

  ```
  sudo -v
  ```

  After that, paste one block at a time and wait for the prompt to come back before pasting the next. If sudo stops to ask for a password in the middle of a pasted block, it swallows the lines that follow. Sudo remembers your password for 15 minutes; after a long break, run `sudo -v` again.
- Only work on the eight fiber ports. Never change `eno4`, and never edit `/etc/network/interfaces` itself.
- One card, one port, one address range. Two ports must never get addresses in the same 192.168.X range.
- Never add a `gateway` line to a card's settings file. The server reaches the internet through `eno4`, and a second gateway breaks that.
- If a step fails twice, stop and ask for help (see the end of this page).

## Your card

Fill these in first. Every command on this page updates to match, so you can copy them as they are.

Card number

Address on tag The tag shows the last two numbers: 69.19 means 192.168.69.19.

Server port Step 1 tells you which one.

dbar, port

fiber, 850 nm

Card

The server port takes the card's first three numbers and ends in .1. The /24 means the whole 192.168.X range is reached through that one port.

## Steps for one card

1. ### Connect the card and find its port

   Power the card on and plug its fiber into a free port on the back of dbar. Then switch on all eight fiber ports. This is safe to repeat; ports that are already set up don't change.

   ```
   for p in ens1f0 ens1f1 ens1f2 ens1f3 ens3f0 ens3f1 ens3f2 ens3f3; do sudo ip link set $p up; done
   ```

   This waits five seconds, then lists every fiber port with its link and address:

   ```
   sleep 5; for p in ens1f0 ens1f1 ens1f2 ens1f3 ens3f0 ens3f1 ens3f2 ens3f3; do printf "%-8s link=%-4s addr=%s\n" $p "$(sudo ethtool $p | awk '/Link detected/{print $3}')" "$(ip -br -4 addr show $p | awk '{print $3}')"; done
   ```

   It passes if exactly one port shows `link=yes` with nothing after `addr=`. That's your card's port: choose it under [Your card](#card). Ports that already show an address belong to cards set up earlier. If card 2 were on `ens1f0`, you'd see something like:

   ```
   ens1f0   link=yes  addr=
   ens1f1   link=no   addr=
   ens1f2   link=no   addr=
   ens1f3   link=no   addr=
   ens3f0   link=no   addr=
   ens3f1   link=no   addr=
   ens3f2   link=no   addr=
   ens3f3   link=yes  addr=192.168.69.1/24
   ```

   To be sure which socket that is, make its light blink for ten seconds, then label both ends of the fiber with the card number:

   ```

   ```

   [ ] Step 1 passed for card
2. ### Check the light

   ```

   ```

   It passes if it shows `850nm` and a receive power between about −3 and −10 dBm. Card 1 reads −5.31 dBm. Some modules don't report power; if that line is missing, carry on. Write the power in the record table.

   [ ] Step 2 passed for card
3. ### Give the server port its address

   ```

   ```

   The first line clears any old address on only. Never run it on `eno4`.

   It passes if the third line shows ` UP ` and the last line starts with ` dev src `. That last line is Linux confirming it will reach the card through the right port, from the right address.

   [ ] Step 3 passed for card
4. ### Ping the card

   ```

   ```

   It passes if you see `5 received, 0% packet loss` with times around 0.1 ms. The first reply can be a little slower: that's Linux asking the card for its hardware (MAC) address. The second command should show ` lladdr REACHABLE` (`STALE` is fine too). Write the MAC and the ping result in the record table.

   The expected MAC follows card 1's pattern: the card's own address written in hexadecimal. If yours is different but the ping works, just record what you see.

   [ ] Step 4 passed for card
5. ### Make it permanent

   Create the card's settings file. Paste this whole block at once; it's a single command:

   ```

   ```

   Hand the port over to the networking service:

   ```

   ```

   Check it again:

   ```

   ```

   It passes if is back on and the ping works. The file has no gateway line on purpose.

   [ ] Step 5 passed for card
6. ### Record and label

   Click [Add to record table](#card) under Your card, then fill in the light level, MAC, and ping result from steps 2 and 4, the settings file name, and your name and the date. Check that the fiber is labeled at both ends. Then pick the next card number under Your card and start again at step 1.

   [ ] Step 6 passed for card

## If something goes wrong

| What you see | What it means | What to do |
| --- | --- | --- |
| No port shows `link=yes` in step 1 | No light gets through: the card is off, a plug isn't seated, or the two fibers are swapped. | Re-seat both SFP modules and the fiber. Swap the two LC plugs at one end only. Then read the light with step 2. |
| Receive power is `-inf` or below about −17 dBm | Almost no light arrives. | Same as above. Also check that the card's module is 850 nm multimode, like the server's. |
| `RTNETLINK answers: File exists` | The address is already on the port. | Run the step 3 block again from its first line, which clears the old address. |
| `ip route get` names a different port, or `eno4` | Another port already uses this 192.168.X range. | Stop. Run `ip -br addr` and send the output to your supervisor. |
| `100% packet loss`, and `ip neigh` shows `FAILED` or `INCOMPLETE` | The card didn't answer at all: wrong address, card not running, or a different range. | Re-read the tag and power-cycle the card. Then list every device that answers in that range:      Or listen for 30 seconds while you reset the card: |
| Ping fails, but `ip neigh` shows a MAC | The card is reachable, but its firmware doesn't answer ping. | Record it and tell your supervisor. Data may still flow normally. |
| `ifup: interface … already configured` | The networking service already set up this port. | Restart it: |
| The address is gone after a reboot | The settings file is missing or has a typo. | Show it and compare with step 5: |

If `arp-scan` isn't found, install it first with `sudo apt install -y arp-scan`.

## Undo a card's setup

To remove a card's settings, for example after choosing the wrong port, run this, then start again from step 3:

```

```

If `ifdown` says the interface isn't configured, that's fine; the last line still clears the address.

## Record

Entries are saved in this browser only. Use Copy table to paste them into the lab log or a spreadsheet.

| Card | Card IP | Server port | Server IP | MAC | Light (dBm) | Ping | Settings file | Checked by, date |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Check all cards at once

Run this after a reboot or before a data run. It pings every card in the record table, prints OK or FAIL for each, then lists the address on every fiber port.

```

```

## How a ping reaches the card

1. **Light.** Each SFP module has to see light from the other end before the link comes up. Steps 1 and 2 check this.
2. **Address and route.** Giving the port tells Linux that every .x address is reached through that port. Step 3 checks this with `ip route get`.
3. **ARP.** Before talking to the card, Linux broadcasts "who has ?" and the card answers with its hardware (MAC) address. `ip neigh` in step 4 shows that answer.
4. **Ping.** An echo request goes to that MAC address and the card sends back an echo reply. On a direct fiber the round trip takes about 0.1 ms.

## If you get stuck

Stop after a step fails twice. Don't try other commands or edit other files. Copy the full terminal output of the failing step, including the command you typed, and send it to your supervisor with the card number and port.