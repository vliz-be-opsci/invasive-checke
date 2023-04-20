FROM mambaorg/micromamba:bullseye-slim
WORKDIR /code
COPY ./environment.yml /code/environment.yml
COPY ./app /code/app
COPY ./invasive_checker /code/app/invasive_checker/
RUN micromamba install -y -n base -f /code/environment.yml && micromamba clean --all --yes
