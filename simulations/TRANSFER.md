# O2 File Transfer -- global-scholars

Exact commands for moving code and results between your machine and Harvard O2.
Source: HMS RC "File Transfer" wiki
(https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1588662157/File+Transfer).

**Identity for this project**

| Item | Value |
|------|-------|
| Transfer host | `transfer.rc.hms.harvard.edu` (port 22) |
| Username | `dma12` (lowercase HMS ID) |
| Auth | **SSH key** (set up once, below) -- no passwords |
| Scratch root | `/n/scratch/users/d/dma12/` |
| Home quota | 100 GiB (final results only) |

## ⚠️ Rules (do not skip)

- **Never transfer while connected to the HMS VPN.** The VPN is not built for
  bulk data and RC will terminate the transfer. Disconnect VPN first.
- The transfer servers have **no modules and no job scheduler** -- run `scp` /
  `rsync` / `cp` directly after logging in; do not `sbatch` there.
- Always give the **full `/n/scratch/...` path**. Tools default to `$HOME`
  (100 GiB), which fills fast.
- Scratch is **temporary and not backed up** -- keep source code in git and pull
  only *final* results back to your machine.

## One-time SSH-key setup (no more passwords)

```bash
# 1. Create a key if you don't already have one (Ed25519 recommended).
ls ~/.ssh/id_ed25519.pub 2>/dev/null || ssh-keygen -t ed25519 -C "dma12@o2"

# 2. Install your public key on O2 (prompts for your HMS password ONCE).
#    Do this off the VPN.
ssh-copy-id dma12@transfer.rc.hms.harvard.edu

# 3. (Optional) add a short alias so later commands are terse.
cat >> ~/.ssh/config <<'EOF'

Host o2-transfer
    HostName transfer.rc.hms.harvard.edu
    User dma12
    Port 22
EOF

# Test: should log in with no password prompt.
ssh o2-transfer 'echo connected; hostname'
```

After this, every command below runs key-only. `o2-transfer` == the full host.

## Push code TO O2 (your machine -> O2 home)

Simulation code lives in git; the usual path is to `git pull` on O2. Use rsync
only for files not in git (e.g. large fixed inputs).

```bash
# Preferred: version-controlled code
#   (on O2)  cd ~/global-scholars && git pull

# rsync a directory of inputs to scratch (NOT home):
rsync -avz --progress \
  ./local-inputs/ \
  o2-transfer:/n/scratch/users/d/dma12/global-scholars/inputs/

# Single file to scratch:
scp ./bigfile.rds \
  o2-transfer:/n/scratch/users/d/dma12/global-scholars/inputs/
```

## Pull RESULTS FROM O2 (O2 home -> your machine)

Only the final combined result per run lives in home `results/`. Pull that, not
scratch intermediates.

```bash
# One study's final results (home dir on O2 -> local):
rsync -avz --progress \
  o2-transfer:~/global-scholars/simulations/<study-name>/results/ \
  ./simulations/<study-name>/results/

# A single run's result file:
scp o2-transfer:~/global-scholars/simulations/<study-name>/results/<run-id>.rds \
  ./simulations/<study-name>/results/
```

## Large / long transfers (keep alive)

A transfer dies if your laptop disconnects. Two options from the wiki:

```bash
# Option A: nohup + background. Output goes to nohup.out; safe to disconnect.
nohup rsync -avz \
  o2-transfer:~/global-scholars/simulations/ ./simulations/ &

# Option B: screen (resume the session later from anywhere).
#   (on O2)  screen           # start a session
#            rsync -av --remove-source-files SRC DST
#            Ctrl-A d          # detach; reattach later with:  screen -r
```

`rsync` is resumable: if it stops, re-run the same command and it continues from
the breakpoint.

## Batch copy on the cluster (large sets, as a job)

Do **not** run large copies on the O2 login nodes. Submit a copy job instead
(login nodes have the scheduler; transfer nodes do not):

```bash
# (on an O2 login node)
sbatch -p short -t 0-12:00 --wrap="rsync -a /n/scratch/users/d/dma12/global-scholars/ ~/global-scholars/archive/"
```

## Quick reference

| Task | Command |
|------|---------|
| Log in to transfer node | `ssh o2-transfer` |
| Push inputs to scratch | `rsync -avz ./in/ o2-transfer:/n/scratch/users/d/dma12/global-scholars/in/` |
| Pull one study's results | `rsync -avz o2-transfer:~/global-scholars/simulations/<study>/results/ ./simulations/<study>/results/` |
| Keep long transfer alive | prefix `nohup ... &` or use `screen` |
| Resume interrupted rsync | re-run the same `rsync` command |
