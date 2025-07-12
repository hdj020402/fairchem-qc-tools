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
Ensure all dependencies are installed:

```bash
uv sync
```

Start the service using Uvicorn:

```bash
uv run uvicorn service:app --host 0.0.0.0 --port 8000
```

The service will be available at `http://localhost:8000`.

---

### Example Requests

#### Run via CLI Interface
Send a POST request to `/run-cli` with a JSON payload containing your task configuration:

```json
{
  "structure": {
    "file_path": "path/to/structure.xyz",
    "charge": 0,
    "spin": 1
  },
  "task": {
    "optimize_geometry": {
      "run": true,
      "steps": 100
    }
  }
}
```

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