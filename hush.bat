@echo off
rem HUSH launcher (Windows) - runs the CLI from the project's venv without
rem needing to activate it. Works from any working directory: %~dp0 is this
rem script's own folder (the repo root), where .venv lives.
"%~dp0.venv\Scripts\hush.exe" %*
