FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# Copy the project into the image

WORKDIR /app
COPY pyproject.toml /app
COPY uv.lock /app
RUN uv sync
COPY . /app
# RUN mkdir /model
CMD ["uv", "run", "main.py"]
