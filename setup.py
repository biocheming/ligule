from setuptools import setup

setup(
    name='ligule',
    version="0.1",
    py_modules=["ligule"],
    install_requires=[
        'Click',
        'Numpy',
        'Pandas',
        'MdTraj'
        ],
    entry_points="""
        [console_scripts]
        ligule=ligule:cli
    """,
)
