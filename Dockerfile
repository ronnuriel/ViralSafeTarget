FROM mambaorg/micromamba:1.5.10

WORKDIR /workspace
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace
RUN micromamba install -y -n base -f environment.yml && micromamba clean --all --yes

EXPOSE 8501 8888
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0"]
