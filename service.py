import yaml, os, logging
import tempfile
from fastapi import FastAPI, UploadFile, File, Body, Form, Query, HTTPException
from fastapi.responses import StreamingResponse
from contextlib import asynccontextmanager
from core.loader import PredictorLoader
from core.calculator import QCCalculator
from core.utils import (gen_atoms,
                        gen_xyz,
                        setup_logger,
                        build_task_log_info,
                        get_config_value,
                        update_dict)
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request

CONFIG_PATH = "global_config.yml"
with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config: dict[str, dict] = yaml.safe_load(f)

LOG_DIR = "logs"
os.makedirs(LOG_DIR, exist_ok=True)

SERVICE_LOG_PATH = os.path.join(LOG_DIR, "service.log")

service_logger = setup_logger("ServiceLogger", SERVICE_LOG_PATH)
service_logger.info('Loading calculator...')
loader = PredictorLoader(
    model_size=config["model"]["size"],
    device=config["model"]["device"]
)
service_logger.info('Done!')

@asynccontextmanager
async def lifespan(app: FastAPI):
    service_logger.info("QC Service started.")
    service_logger.info(f'Model size: {config["model"]["size"]}, Device: {config["model"]["device"]}')
    yield
    service_logger.info("QC Service shutting down.")
    if config["model"].get("device") == "cuda":
        import torch
        torch.cuda.empty_cache()
        service_logger.info("CUDA cache cleared.")

app = FastAPI(lifespan=lifespan)

class LoggingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        service_logger.info(f"Request received: {request.method} {request.url}")
        response = await call_next(request)
        service_logger.info(f"Request completed with status code: {response.status_code}")
        return response

app.add_middleware(LoggingMiddleware)

@app.post("/run-cli")
async def run_qc_cli(task_config: dict = Body(...)):
    updated_config = update_dict(config, task_config)
    return await run_qc_public(updated_config)

@app.post("/run-web")
async def run_qc_web(
    structure_file: UploadFile = File(...),
    task_config: str = Form(...)
    ):
    task_config = yaml.safe_load(task_config)
    updated_config = update_dict(config, task_config)
    # TODO: edit config via front end
    # ?: where to save log files?
    # ?: when a structure is submitted from web, a temporary file is created at `./`,
    # ?: and thus the log file will be saved to `./logs/task_id/task_id.log`
    return await run_qc_public(updated_config, structure_file)

async def run_qc_public(
    updated_config: dict,
    structure_file: UploadFile | None = None,
    ) -> dict[str, str]:
    if structure_file:
        filename = structure_file.filename
    else:
        structure_file_path = get_config_value(updated_config, "structure", "file_path")
        filename = structure_file_path

    log_file, task_id = build_task_log_info(filename)
    current_task_logger = setup_logger(task_id, log_file)

    try:
        if structure_file:
            current_task_logger.info(f"Received uploaded file: {filename}")
            suffix = os.path.splitext(structure_file.filename)[1]
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmpfile:
                content = await structure_file.read()
                tmpfile.write(content)
                tmpfile_path = tmpfile.name
        else:
            current_task_logger.info(f"Using structure file from config: {filename}")
            structure_file_path = get_config_value(updated_config, "structure", "file_path")
            if not os.path.isfile(structure_file_path):
                raise FileNotFoundError(f"Structure file not found: {structure_file_path}")
            tmpfile_path = structure_file_path

        charge = get_config_value(updated_config, 'structure', 'charge', default=0)
        spin = get_config_value(updated_config, 'structure', 'spin', default=1)
        atoms = gen_atoms(tmpfile_path, charge, spin)
        qc = QCCalculator(atoms=atoms, calc=loader.calc, logger=current_task_logger)
        current_task_logger.info(
            f'Model size: {config["model"]["size"]}, Device: {config["model"]["device"]}')
        current_task_logger.info(f'Input orientation: \n{gen_xyz(atoms)}')
        current_task_logger.info(f"Charge: {atoms.info['charge']}, Spin: {atoms.info['spin']}")

        for task, args in updated_config['task'].items():
            args: dict
            if args['run']:
                kwargs = {k: v for k, v in args.items() if k != 'run'}
                getattr(qc, task)(**kwargs)

        if structure_file:
            os.unlink(tmpfile_path)

        return {
            "status": "success",
            "task_id": task_id,
            "task_log_path": log_file
        }

    except FileNotFoundError as e:
        error_msg = f"Structure file not found: {e}"
        service_logger.error(error_msg, exc_info=True)
        current_task_logger.error(error_msg, exc_info=True)
        raise HTTPException(status_code=404, detail=error_msg)

    except PermissionError as e:
        error_msg = f"Permission denied: {e}"
        service_logger.error(error_msg, exc_info=True)
        current_task_logger.error(error_msg, exc_info=True)
        raise HTTPException(status_code=403, detail=error_msg)

    except Exception as e:
        error_msg = f"Unexpected error: {e}"
        service_logger.error(error_msg, exc_info=True)
        current_task_logger.error(error_msg, exc_info=True)
        raise HTTPException(status_code=500, detail=error_msg)

    finally:
        if current_task_logger:
            for handler in current_task_logger.handlers[:]:
                if isinstance(handler, logging.FileHandler):
                    handler.close()
                current_task_logger.removeHandler(handler)

@app.get("/logs/{task_id}")
async def get_task_log(task_id: str) -> dict[str, str]:
    log_dir = os.path.join("logs", task_id)
    log_file = os.path.join(log_dir, f"{task_id}.log")

    if not os.path.exists(log_file):
        raise HTTPException(status_code=404, detail="Log file not found")

    with open(log_file, "r", encoding="utf-8") as f:
        content = f.read()

    return {
        "task_id": task_id,
        "log": content
    }

@app.get("/download/{task_id}")
async def download_task_folder(
    task_id: str,
    archive_format: str = Query("zip", regex="^(zip|tar.gz)$"),
    ) -> StreamingResponse:
    import io, tarfile
    from zipfile import ZipFile, ZIP_DEFLATED
    from core.utils import get_folder_size

    task_dir = os.path.join("logs", task_id)

    if not os.path.isdir(task_dir):
        raise HTTPException(status_code=404, detail="Task folder not found")

    max_size = 1024 * 1024 * 1024
    folder_size = get_folder_size(task_dir)
    if folder_size > max_size:
        raise HTTPException(status_code=413, detail="Task folder too large to download")

    try:
        if archive_format == "zip":
            zip_buffer = io.BytesIO()
            with ZipFile(zip_buffer, mode="w", compression=ZIP_DEFLATED) as zip_file:
                for root, dirs, files in os.walk(task_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.relpath(file_path, start="logs")
                        zip_file.write(file_path, arcname=arcname)

            zip_buffer.seek(0)
            content = zip_buffer.getvalue()
            media_type = "application/zip"
            headers = {
                "Content-Disposition": f"attachment; filename={task_id}.zip"
            }

        elif archive_format == "tar.gz":
            tar_buffer = io.BytesIO()
            with tarfile.open(fileobj=tar_buffer, mode="w:gz") as tar:
                tar.add(task_dir, arcname=os.path.basename(task_dir))
            tar_buffer.seek(0)
            content = tar_buffer.getvalue()
            media_type = "application/gzip"
            headers = {"Content-Disposition": f"attachment; filename={task_id}.tar.gz"}

        return StreamingResponse(io.BytesIO(content), media_type=media_type, headers=headers)

    except Exception as e:
        service_logger.error(f"Error creating {archive_format} archive for task {task_id}: {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Failed to create {archive_format} archive: {str(e)}")

