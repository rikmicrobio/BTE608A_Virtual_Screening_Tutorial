# 🧪 Virtual Screening Practical
## Installation Notes

Before starting the practical, please install the following tools:

- **Ubuntu/Linux or WSL** — for command-line work
- **OpenBabel** — chemical file conversion and 3D structure generation
- **DataWarrior** — molecular property and toxicity analysis

---

## 1. Set Up Ubuntu / WSL

### If you already have Ubuntu/Linux

You can directly use your existing **Ubuntu/Linux** installation.

### If you are using Windows

Install **WSL + Ubuntu**.

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

Restart your computer if prompted.

After restarting, open **Ubuntu** from the Windows Start Menu.

Check that Ubuntu is working:

```bash
uname -a
```

---

## 2. Install OpenBabel

Open your **Ubuntu/Linux terminal** and run:

```bash
sudo apt-get update
```

Then install OpenBabel:

```bash
sudo apt-get install -y openbabel
```

Check that OpenBabel is installed correctly:

```bash
obabel -V
```

You should see something similar to:

```text
Open Babel 3.x.x
```

### Quick Test

Run:

```bash
obabel -L formats
```

This should display a list of chemical file formats supported by OpenBabel.

### Additional Installation Guide

For more information:

[OpenBabel Installation on Ubuntu — Bioinformatics Review](https://bioinformaticsreview.com/20240926/tutorial-how-to-install-openbabel-on-ubuntu-linux/)

---

## 3. Install DataWarrior

DataWarrior is available for **both Windows and Linux**.

Choose the version appropriate for your operating system:

- 🪟 **Windows** → Download the Windows version
- 🐧 **Linux/Ubuntu** → Download the Linux version

[Download DataWarrior](https://openmolecules.org/datawarrior/download.html)

Install DataWarrior normally for your operating system and make sure that it opens successfully.

---

## 4. Final Check ✅

Before coming to the practical, make sure the following are working:

```text
☑ Ubuntu / WSL is working
☑ OpenBabel is installed
☑ DataWarrior is installed
```

Test OpenBabel one more time:

```bash
obabel -V
```

If you can see the OpenBabel version, **you are ready to start! 🚀**

---

## 🆘 Having Problems?

If something does not work:

1. **Do not randomly change the commands.**
2. Copy the complete error message from the terminal.
3. Take a screenshot if necessary.
4. Ask for help before proceeding.

We will troubleshoot it together.
