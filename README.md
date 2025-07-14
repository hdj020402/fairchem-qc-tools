# QC Service README

## Overview
The QC Service is a FastAPI-based web service designed to run quantum chemistry calculations on molecular structures. It supports both CLI and web interfaces for submitting tasks, dynamically loads models for prediction, and logs execution details for debugging and tracking purposes.

## Features
- **FastAPI Framework**: Provides RESTful API endpoints for handling requests.
- **Quantum Chemistry Calculations**:
  - Accepts molecular structure files (e.g., `.xyz`, `.gjf`, etc.).
  - Supports configuration via request payloads or global config files.
  - Logs detailed task information per execution.
- **Model Management**:
  - Uses the UMA model from [FairChem](https://github.com/facebookresearch/fairchem) for molecular property prediction.
  - Loads predictor models based on configuration (`PredictorLoader`).
  - Clears CUDA cache gracefully during shutdown if GPU is used.
- **Logging System**:
  - Global logging for the service lifecycle.
  - Task-specific logging with unique identifiers.
- **Error Handling**:
  - Gracefully handles missing files, permission issues, and unexpected errors.
  - Returns appropriate HTTP status codes and error messages.

---

## Architecture

### Core Components
1. **Main Application (`service.py`)**:
   - Sets up FastAPI routes and middleware.
   - Initializes global configurations and model loader.
   - Manages application lifespan events (startup/shutdown).

2. **Core Modules**:
   - `core/loader.py`: Handles loading of machine learning models.
   - `core/calculator.py`: Contains the `QCCalculator` class responsible for executing quantum chemistry tasks.
   - `core/utils.py`: Utility functions for data processing, logging, and configuration management.

3. **Configuration**:
   - Uses a YAML file (`global_config.yml`) for default settings.
   - Configurable parameters include:
     - Model size and device (CPU/GPU).
     - Structure file path, charge, spin.
     - Task types and their arguments.

4. **Endpoints**
   - **POST `/run-cli`**:
     - Accepts a JSON body with task configuration.
     - Merges input with global config and runs the calculation.
   - **POST `/run-web`**:
     - Accepts an uploaded structure file from a web interface.
     - Stores it temporarily and processes it.
   - **GET `/logs/{task_id}`**:
     - Retrieves log content for a specific task by its ID.

---

## Usage

### Running the Service

#### With `uv`
Ensure all dependencies are installed:

```bash
uv sync
```

Start the service using Uvicorn:

```bash
uv run uvicorn service:app --host 0.0.0.0 --port 8000
```

The service will be available at `http://localhost:8000`.

#### With `conda`
Create a new environment and ensure all dependencies are installed:

```bash
conda create -n fairchem-qc-tools python=3.12
conda activate fairchem-qc-tools
pip install -r requirements.txt
```

Start the service using Uvicorn:

```bash
uvicorn service:app --host 0.0.0.0 --port 8000
```

The service will be available at `http://localhost:8000`.

---

### Example Requests

#### Running via CLI Interface

The CLI interface allows you to submit quantum chemistry tasks directly from the command line. Unlike running the full web service, using the CLI requires fewer dependencies and does not require starting a FastAPI server.

##### Minimum Required Files
To run the CLI, the following files are **required**:
- `cli.py`: The main CLI script.
- `task_config.yml`: Configuration file defining task parameters.
- `core/utils.py`: Contains utility functions for configuration management, logging, and common data operations.

These files can be placed in the same directory for simplicity.

##### Setup (Optional)
If you want to create a new environment, the following minimal dependencies are required:

```plain text
ase==3.25.0
numpy==2.3.1
PyYAML==6.0.2
Requests==2.32.4
```

These are significantly fewer than the full service dependencies, making the CLI a lightweight option for running quantum chemistry tasks without the overhead of the web server or machine learning model dependencies.

##### Execution
You can execute the CLI by specifying configuration options directly from the command line:

```bash
python cli.py --service.host <host_ip> --structure.file_path <structure_file>
```

Alternatively, you can modify the `task_config.yml` file in the same directory to predefine your settings.

##### Output and Logging
After submission, a `logs/` directory will be created within the structure file's directory. Each task generates a subdirectory named according to the format `<task_timestamp_uid>/`, and the log file is saved as `<task_timestamp_uid>.log`.

> ⚠️ Permissions Note: If the service runs under elevated privileges (e.g., root), generated logs may be read-only for regular users. Consider adjusting permissions or running the service under a standard user if this causes issues.

##### Summary
The CLI mode offers a lightweight way to perform quantum chemistry calculations without launching the full web service. It supports direct command-line arguments or configuration through `task_config.yml`, and it generates structured logs for tracking and debugging.

#### Upload via Web Interface
Use `/run-web` to upload a molecular structure file directly:

```bash
curl -X POST http://localhost:8000/run-web \
  -H "Content-Type: multipart/form-data" \
  -F "structure_file=@your_structure.xyz"
```

---

### Retrieve Task Logs
To retrieve logs for a completed task, use:

```bash
GET /logs/{task_id}
```

This returns the full log content for the specified task.

---

## Error Handling
The service includes robust error handling:
- **404 Not Found**: If a structure file is not found.
- **403 Forbidden**: If there's a permission issue.
- **500 Internal Server Error**: For any unexpected runtime exceptions.

---

## Logging
- All service-level logs are written to `logs/service.log`.
- Each task generates a dedicated log file under `logs/{task_id}/{task_id}.log`.

---

## Future Enhancements
- Add support for asynchronous task execution and result polling.
- Implement a frontend dashboard for monitoring and managing tasks.
- Expand supported file formats and calculation types.

---

## License
TBD

--- 

## Contact
For questions or contributions, please reach out to the project maintainers.