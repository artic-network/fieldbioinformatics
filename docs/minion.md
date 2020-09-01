---
title: minion
summary: Outline of minion workflows.
authors:
  - Will Rowe
  - Nick Loman
date: 2020-09-01
---

# Core pipeline

## About

This page describes the core pipeline which is called using the `artic minion` command.

There are **2 workflows** baked into the core pipeline, one which uses signal data (via [nanopolish](https://github.com/jts/nanopolish)) and one that does not (via [medaka](https://github.com/nanoporetech/medaka)). As the workflows are identical in many ways, this page will describe the pipeline as whole and notify the reader when there is dfferent behaviour between the two workflows.
