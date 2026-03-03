# HInt_entry.py
import multiprocessing as mp
from absl import app
from HInt.HInt import main as HInt_main

def main():
    mp.set_start_method("spawn", force=True)  # important pour propagation Ctrl+C
    app.run(HInt_main)

if __name__ == "__main__":
    main()
