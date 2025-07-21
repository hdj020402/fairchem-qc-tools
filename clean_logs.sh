#!/bin/bash

cd "$(dirname "$0")" || exit 1

LOG_DIR="logs"
SERVICE_LOG="service.log"

find "$LOG_DIR" -maxdepth 1 -type f ! -name "$SERVICE_LOG" -exec rm -f {} \; && echo "Deleted old log files."

find "$LOG_DIR" -maxdepth 1 -type d -not -path "$LOG_DIR" -mtime +3 -exec rm -rf {} \; && echo "Deleted folders older than 3 days."
