# MDZen MCP Server Docker Image
# Provides all scientific dependencies (AmberTools, OpenMM, RDKit, etc.)
# with FastMCP servers exposed via Streamable HTTP transport

FROM condaforge/mambaforge:latest

LABEL maintainer="Matsunaga Lab"
LABEL description="MDZen MCP Server - Molecular Dynamics Setup Tools"
LABEL version="0.2.0"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV MDZEN_LOG_LEVEL=INFO

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install scientific dependencies via conda-forge
# These are difficult to install via pip and require specific builds
RUN mamba install -y -c conda-forge \
    python=3.11 \
    openmm \
    rdkit \
    ambertools \
    mdanalysis \
    biopython \
    pandas \
    numpy \
    scipy \
    openblas \
    pdbfixer \
    packmol \
    smina \
    && mamba clean -afy

# Fix packmol-memgen numpy compatibility (NumPy 1.24+ removed np.float)
RUN SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])") && \
    if [ -f "$SITE_PACKAGES/packmol_memgen/lib/pdbremix/v3numpy.py" ]; then \
        sed -i.bak "s/np\.float)/float)/g; s/np\.int)/int)/g" \
            "$SITE_PACKAGES/packmol_memgen/lib/pdbremix/v3numpy.py"; \
    fi

# Copy requirements for pip dependencies (MCP server only, no ADK)
COPY docker/requirements-mcp.txt /app/requirements.txt

# Install Python dependencies via pip
RUN pip install --no-cache-dir -r requirements.txt

# Copy MCP servers and common utilities
COPY servers/ /app/servers/
COPY common/ /app/common/
COPY src/mdzen/config.py /app/src/mdzen/config.py
COPY src/mdzen/schemas.py /app/src/mdzen/schemas.py
COPY src/mdzen/utils.py /app/src/mdzen/utils.py

# Set PYTHONPATH to include source directories
ENV PYTHONPATH=/app:/app/src

# Copy entrypoint script
COPY docker/entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

# Create workdir for user files
RUN mkdir -p /workdir
WORKDIR /workdir

# Expose MCP server port (Streamable HTTP)
EXPOSE 3000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:3000/health || exit 1

# Run MCP server
ENTRYPOINT ["/app/entrypoint.sh"]
