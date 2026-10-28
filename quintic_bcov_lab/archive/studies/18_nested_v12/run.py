#!/usr/bin/env python3
"""Prediction-only orchestrator. Reads data/ and never imports validation/."""
from predict_terminal import run as terminals
from predict_terminal_values import run as reconstruction
from predict_diagonal import run as diagonals
from predict_r6 import run as r6

if __name__=='__main__':
 terminals();reconstruction();diagonals();r6()
