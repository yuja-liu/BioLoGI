# BioLoGI

![Biologic LOGO](img/logo.png)

BioLoGI (**Bio**logical **Lo**gic **G**ate **I**nference system) is a computational tool that automates the modeling
of logic gate-like interaction in syntheic genetic circuits and human transcriptome.

## Installation

Install by [pip](https://pip.pypa.io/en/stable/):

```sh
pip install https://github.com/RandolphLiu/BioLoGI/releases/download/v0.1/biologi-0.1-py3-none-any.whl
```

One of the dependencies, `theano` may require python development package depending on your platform.
If `theano` throws an error, try:

```sh
# Ubuntu
sudo apt install python3-dev

# CentOS, RHEL or Fedora
sudo dnf install python3-devel
```

## Quick Start

Modeling the 3OC6HSL receiver device, which is a synthetic circuit uses the 3OC6HSL inducible promoter pLux. 
See iGEM registry: [BBa\_T9002](https://parts.igem.org/Part:BBa_T9002)

After installing BioLoGI:

```sh
git clone https://github.com/RandolphLiu/BioLoGI
cd BioLoGI
biologi -c data/BBa_T9002/BBa_T9002_config.json
```

The output will include details of all canditate models and the best model, 
the dose-response curve and SBOL genetic circuit figure.
