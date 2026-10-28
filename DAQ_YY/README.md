# DAQ card link setup on dbar

How to connect a digitizer card to the dbar server over fiber, give the link an address, ping the card, and make the setting permanent. Do one card at a time; each takes about ten minutes.

## Setup at a glance

| Item | Value |
|---|---|
| Server | dbar (Dell PowerEdge R740, Debian 12) |
| Fiber ports | `ens1f0`–`ens1f3` (network card in riser slot 1) and `ens3f0`–`ens3f3` (riser slot 3) |
| Lab network | `eno4`, 192.168.35.100. Never change it. |
| Card settings | One file per card in `/etc/network/interfaces.d/` |
| Address rule | Tag `69.19` means the card is `192.168.69.19`; the server port gets `192.168.69.1/24` |
| Already done | Cards 1–3; see the [record](#record) |

## Before you start

> [!IMPORTANT]
> Run `sudo -v` once and type your password. Then paste **one code block at a time** and wait for the prompt to come back before pasting the next. If sudo asks for a password in the middle of a pasted block, it swallows the lines that follow. Sudo remembers your password for 15 minutes; after a long break, run `sudo -v` again.

> [!WARNING]
> Only work on the eight fiber ports. Never change `eno4`, never edit `/etc/network/interfaces` itself, and never add a `gateway` line to a card's settings file. The server reaches the internet through `eno4`, and a second gateway breaks that.

- One card, one port, one address range: two ports must never have addresses in the same `192.168.X` range.
- The commands use shell variables (`CARD`, `CARD_IP`, `PORT`, and so on). They only exist in the terminal where you set them. If you open a new terminal or log in again, run step 0 and the `PORT=` line from step 1 again.
- If a step fails twice, stop and see [If you get stuck](#if-you-get-stuck).

## Step 0: Set your card's values

Edit the first two lines for your card, then paste the block:

```bash
CARD=2                    # card number
CARD_IP=192.168.70.19     # from the tag: "70.19" means 192.168.70.19
NET=${CARD_IP%.*}         # the first three numbers, e.g. 192.168.70
SERVER=$NET.1             # server side of this link
[ "${CARD_IP##*.}" = 1 ] && SERVER=$NET.100   # only if the card's own address ends in .1
echo "card $CARD: $CARD_IP   server: $SERVER/24"
```

**Passes if** the printed line matches the card's tag.

## Step 1: Connect the card and find its port

Power the card on and plug its fiber into a free port on the back of dbar. Then switch on all eight fiber ports. This is safe to repeat; ports that are already set up don't change.

```bash
for p in ens1f0 ens1f1 ens1f2 ens1f3 ens3f0 ens3f1 ens3f2 ens3f3; do sudo ip link set $p up; done
```

Wait five seconds, then list every fiber port with its link and address:

```bash
sleep 5; for p in ens1f0 ens1f1 ens1f2 ens1f3 ens3f0 ens3f1 ens3f2 ens3f3; do printf "%-8s link=%-4s addr=%s\n" $p "$(sudo ethtool $p | awk '/Link detected/{print $3}')" "$(ip -br -4 addr show $p | awk '{print $3}')"; done
```

**Passes if** exactly one port shows `link=yes` with nothing after `addr=`. That's your card's port. Ports that already show an address belong to cards set up earlier. For example, with card 2 on `ens1f0`:

```text
ens1f0   link=yes  addr=
ens1f1   link=no   addr=
ens1f2   link=no   addr=
ens1f3   link=no   addr=
ens3f0   link=no   addr=
ens3f1   link=no   addr=
ens3f2   link=no   addr=
ens3f3   link=yes  addr=192.168.69.1/24
```

Set `PORT` to the port you found (edit the name), and check it:

```bash
PORT=ens1f0
[[ $PORT == ens[13]f[0-3] ]] && echo "OK: $PORT is a fiber port" || echo "STOP: $PORT is not a fiber port"
```

Make the port's LED blink for ten seconds so you can find the socket, then label both ends of the fiber with the card number and port:

```bash
sudo ethtool -p $PORT 10
```

## Step 2: Check the light

```bash
sudo ethtool -m $PORT | grep -iE 'wavelength|optical power'
```

**Passes if** it shows `850nm` and a receive power between about −3 and −10 dBm. Card 1 reads −5.31 dBm. If the power line is missing, the module doesn't report it; carry on.

## Step 3: Give the server port its address

> [!CAUTION]
> The first line below clears every address on `$PORT`. Never run it with `PORT=eno4`.

```bash
sudo ip addr flush dev $PORT
sudo ip addr add $SERVER/24 dev $PORT
ip -br addr show $PORT
ip route get $CARD_IP
```

**Passes if** the third line shows your port as `UP` with `$SERVER/24`, and the last line reads `<card IP> dev <your port> src <server IP>`. That last line is Linux confirming it will reach the card through the right port, from the right address.

## Step 4: Ping the card

```bash
ping -c 5 $CARD_IP
ip neigh show dev $PORT
printf 'expected MAC: 02:00:%02x:%02x:%02x:%02x\n' ${CARD_IP//./ }
```

**Passes if** you see `5 received, 0% packet loss` with times around 0.1 ms, and `ip neigh` shows the card's IP with `lladdr <MAC> REACHABLE` (`STALE` is fine too). The first reply can be a little slower: that's Linux asking the card for its hardware (MAC) address.

The last line prints the MAC expected from card 1's pattern: the card's own address written in hexadecimal. If yours is different but the ping works, just record what you see.

## Step 5: Make it permanent

Create the card's settings file and show it. Paste this whole block at once:

```bash
sudo tee /etc/network/interfaces.d/daq-card$CARD > /dev/null <<EOF
# DAQ card $CARD ($CARD_IP) on fiber port $PORT
auto $PORT
iface $PORT inet static
    address $SERVER/24
EOF
cat /etc/network/interfaces.d/daq-card$CARD
```

The file must show real values, such as `auto ens1f0`. If you see an empty name or a `$`, your variables aren't set: run step 0 and the `PORT=` line again, then repeat this step.

Hand the port over to the networking service:

```bash
sudo ip addr flush dev $PORT
sudo ifup $PORT
```

Check it again:

```bash
ip -br addr show $PORT
ping -c 3 $CARD_IP
```

**Passes if** `$SERVER/24` is back on your port and the ping works. The file has no gateway line on purpose.

## Step 6: Record and label

Add the card's row to the [record](#record) below, with the light level from step 2 and the MAC and ping result from step 4. Check that the fiber is labeled at both ends. Then go back to step 0 for the next card.

## Record

| Card | Card IP | Server port | Server IP | MAC | Light (dBm) | fiber port | Checked by, date |
|---|---|---|---|---|---|---|---|
| 1 | 192.168.63.19 | ens1f0 | 192.168.63.1/24 | 02:00:c0:a8:3f:13 | -6.86 | 1-a | 2026.10.22|
| 2 | 192.168.59.19 | ens1f1 | 192.168.59.1/24 | 02:00:c0:a8:3b:13 | -7.12 | 1-b | 2026.10.22 |
| 3 | 192.168.110.19  | ens1f2 | 192.168.110.1/24 | 02:00:c0:a8:6e:13 | -7.57 | 1-c | 2026.10.22 |
| 4 | 192.168.109.19 | ens1f3 | 192.168.109.1/24 | 02:00:c0:a8:6d:13 | -7.12 | 1-d | 2026.10.28 |
| 5 | 192.168.60.19 | ens3f0 | 192.168.60.1/24 | 02:00:c0:a8:3c:13 | -6.43 | 3-a | 2026.10.28 |
| 6 | 192.168.92.19 | ens3f1 | 192.168.92.1/24 | 02:00:c0:a8:5c:13 | -6.75 | 3-b | 2026.10.28 |
| 7 | 192.168.65.19 | ens3f2 | 192.168.65.1/24 | 02:00:c0:a8:41:13 | -9.69 | 3-c | 2026.10.28 |
| 8 | 192.168.107.19 | ens3f3 | 192.168.107.1/24 | 02:00:c0:a8:6b:13 | -6.48 | 3-d | 2026.10.28 |
<img width="374" height="261" alt="fiber_IO" src="https://github.com/user-attachments/assets/ea9696b1-be5b-4dd4-863c-a669c7b880aa" />

## Show the record from the server

Prints one row of the record table for every card that has a settings file, ready to paste above. It only reads what is already on dbar and sends nothing to the cards: the MAC comes from the server's ARP table and the light from the SFP module. A MAC shows `?` if the server hasn't talked to that card recently. The fiber port label is only on the physical label, so fill it in by hand.

```bash
sudo -v
for f in /etc/network/interfaces.d/daq-card*; do
  card=${f##*daq-card}
  port=$(awk '/^auto/{print $2}' "$f")
  srv=$(awk '/address/{print $2}' "$f")
  ip=$(sed -n 's/^# DAQ card .*(\(.*\)) on.*/\1/p' "$f")
  mac=$(ip neigh show "$ip" dev "$port" | awk '{for(i=1;i<NF;i++) if($i=="lladdr") print $(i+1)}')
  light=$(sudo ethtool -m "$port" 2>/dev/null | grep -iE 'receiver signal|rcvr signal|rx power' | grep -oE -- '-?[0-9.]+ dBm' | head -1)
  echo "| $card | $ip | $port | $srv | ${mac:-?} | ${light% dBm} | | |"
done | sort -t'|' -k2 -n
```


## Check all cards at once

Run this after a reboot or before a data run. It finds every card from its settings file, so a new card is included as soon as step 5 is done. It prints OK or FAIL for each card, then lists the address on every fiber port.

```bash
CARDS=$(sed -n 's/^# DAQ card .*(\(.*\)) on.*/\1/p' /etc/network/interfaces.d/daq-card*)
echo "cards:" $CARDS
for ip in $CARDS; do
  if ping -c 2 -W 1 $ip > /dev/null; then echo "$ip  OK"; else echo "$ip  FAIL"; fi
done
ip -br addr | grep ens
```

## If something goes wrong

| What you see | What it means | What to do |
|---|---|---|
| `STOP: ... is not a fiber port` | `PORT` is wrong or not set | Set `PORT=` to the port step 1 found |
| No port shows `link=yes` in step 1 | No light gets through: the card is off, a plug isn't seated, or the two fibers are swapped | Re-seat both SFP modules and the fiber. Swap the two LC plugs at one end only. Then read the light with step 2 |
| Receive power is `-inf` or below about −17 dBm | Almost no light arrives | Same as above. Also check that the card's module is 850 nm multimode, like the server's |
| `RTNETLINK answers: File exists` | The address is already on the port | Run the step 3 block again from its first line, which clears the old address |
| `ip route get` names a different port, or `eno4` | Another port already uses this `192.168.X` range | Stop. Run `ip -br addr` and send the output to your supervisor |
| `100% packet loss`, and `ip neigh` shows `FAILED` or `INCOMPLETE` | The card didn't answer at all: wrong address, card not running, or a different range | Re-read the tag and power-cycle the card. List every device that answers in the range with `sudo arp-scan -I $PORT $NET.0/24` (install it with `sudo apt install -y arp-scan`), or listen while you reset the card with `sudo timeout 30 tcpdump -i $PORT -n -e` |
| Ping fails, but `ip neigh` shows a MAC | The card is reachable, but its firmware doesn't answer ping | Record it and tell your supervisor. Data may still flow normally |
| `ifup: interface ... already configured` | The networking service already set up this port | Run `sudo ifdown $PORT && sudo ifup $PORT` |
| The address is gone after a reboot | The settings file is missing or has a typo | Run `cat /etc/network/interfaces.d/daq-card$CARD` and compare it with step 5 |

## Undo a card's setup

To remove a card's settings, for example after choosing the wrong port. `CARD` and `PORT` must be set to the values you used:

```bash
sudo ifdown $PORT
sudo rm /etc/network/interfaces.d/daq-card$CARD
sudo ip addr flush dev $PORT
```

If `ifdown` says the interface isn't configured, that's fine; the last line still clears the address. Then start again from step 3.

## How a ping reaches the card

1. **Light.** Each SFP module has to see light from the other end before the link comes up. Steps 1 and 2 check this.
2. **Address and route.** Giving the port `192.168.X.1/24` tells Linux that every `192.168.X` address is reached through that port. Step 3 checks this with `ip route get`.
3. **ARP.** Before talking to the card, Linux broadcasts "who has 192.168.X.Y?" and the card answers with its hardware (MAC) address. `ip neigh` in step 4 shows that answer.
4. **Ping.** An echo request goes to that MAC address and the card sends back an echo reply. On a direct fiber the round trip takes about 0.1 ms.

## If you get stuck

Stop after a step fails twice. Don't try other commands or edit other files. Copy the full terminal output of the failing step, including the command you typed, and send it to your supervisor with the card number and port.
