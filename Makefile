DOCKER_CONTEXT ?= colima-hypsipyle
IMAGE ?= hypsipyle:dev
DOCKERFILE ?= Dockerfile.dev
CONTAINER_NAME ?= hypsipyle-dev
PWD := $(shell pwd)
PORTS := --publish 0.0.0.0:80:80/tcp --publish 0.0.0.0:8000:8000/tcp

.PHONY: help context build rebuild run run-detach stop logs shell rm

help:
	@echo "Available targets:"
	@echo "  make context        # Set docker context"
	@echo "  make build          # Build the image"
	@echo "  make rebuild        # Build the image without cache"
	@echo "  make run            # Run container (interactive)"
	@echo "  make run-detach     # Run container detached"
	@echo "  make stop           # Stop and remove container"
	@echo "  make logs           # Follow container logs"
	@echo "  make shell          # Exec a shell in running container"

context:
	docker context use $(DOCKER_CONTEXT)

build:
	docker build -t $(IMAGE) -f $(DOCKERFILE) .

rebuild:
	docker build --no-cache -t $(IMAGE) -f $(DOCKERFILE) .

run: build
	docker container run $(PORTS) -ti -v $(PWD):/app --name $(CONTAINER_NAME) $(IMAGE)

run-detach: build
	docker container run $(PORTS) -d -v $(PWD):/app --name $(CONTAINER_NAME) $(IMAGE)

stop:
	docker container stop $(CONTAINER_NAME) || true
	docker container rm $(CONTAINER_NAME) || true

logs:
	docker logs -f $(CONTAINER_NAME)

shell:
	docker exec -it $(CONTAINER_NAME) /bin/bash

rm:
	docker container rm -f $(CONTAINER_NAME) || true
