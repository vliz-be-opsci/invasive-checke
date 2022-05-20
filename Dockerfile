FROM tiangolo/uvicorn-gunicorn-fastapi:python3.9-slim
WORKDIR /code
COPY ./requirements.txt /code/requirements.txt

RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt
RUN apt-get update && apt-get install -y \
vim\
&& rm -rf /var/lib/apt/lists/*

COPY ./app /code/app
COPY ./invasive_checker /code/app/invasive_checker/
