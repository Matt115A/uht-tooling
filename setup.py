from setuptools import setup, find_packages

# Read dependencies from requirements.txt
with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="uht-tooling",
    version="0.1.0",
    packages=find_packages(include=["scripts", "scripts.*"]),
    install_requires=requirements,
    author="Matt115A",
    description="Tooling for ultra-high throughput screening workflows",
    url="https://github.com/Matt115A/uht-tooling",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Adjust license if needed
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
