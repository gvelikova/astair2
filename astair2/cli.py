import logging
import warnings

import click

import astair2
import astair2.aligner as aligner
import astair2.caller as caller
import astair2.filter as filter
import astair2.finder as finder
import astair2.idbiaser as idbiaser
import astair2.mbias as mbias
import astair2.phred as phred
import astair2.simulator as simulator
import astair2.summary as summary

# TODO make this config properly configurable using command line options
# For example, we could use a global option -v to change the log level to
# DEBUG and produce verbose output for all commands
logging.basicConfig(level=logging.WARNING)
warnings.simplefilter(action="ignore", category=UserWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)


@click.group()
def cli():
    """
    asTair2 (tools for processing cytosine modification sequencing data)
    """
    pass


cli.epilog = """
__________________________________About__________________________________
asTair2 is developed by Gergana V. Velikova. It is based on asTair, which
was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.
This code is made available under the GNU General Public License v3, see
LICENSE.txt for more details.

                                                         Version: __version__
"""
cli.epilog = cli.epilog.replace("__version__", astair2.__version__)


cli.add_command(caller.call)
cli.add_command(phred.phred)
cli.add_command(mbias.mbias)
cli.add_command(finder.find)
cli.add_command(aligner.align)
cli.add_command(simulator.simulate)
cli.add_command(filter.filter)
cli.add_command(summary.summarise)
cli.add_command(idbiaser.idbias)


if __name__ == "__main__":
    cli()
