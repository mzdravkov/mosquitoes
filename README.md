# mosquitos
A playground to test some research ideas for defeating mosquitos.


## Prerequisites:
You have to have python3 installed.

You have to have the **blat** executable in the your PATH.

You can download it from: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/
If you're using bash, you can make it accessible by: 
```bash
$ mkdir ~/bin
$ mv path/to/blat ~/bin
$ chmod +x ~/bin/blat

# then add this to your ~/.bashrc (or .zshrc if you're using zsh)
export PATH=$PATH:~/bin

# Make sure that the changes are sourced
$ source ~/.bashrc
# or source ~/.zshrc (respectively)
```

## Setup

```bash
$ python -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```