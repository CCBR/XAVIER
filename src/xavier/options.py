#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os
from __future__ import print_function
from ccbr_tools.pipeline.util import permissions


def genome_options(parser, user_option, prebuilt):
    """Dynamically checks if --genome option is a valid choice. Compares against a
    list of prebuilt or bundled genome reference genomes and accepts a custom reference
    JSON file.
    @param parser <argparse.ArgumentParser object>:
        Parser object from which an exception is raised not user_option is not valid
    @param user_option <str>:
        Provided value to the xavier run, --genome argument
    @param prebuilt list[<str>]:
        List of pre-built reference genomes
    return user_option <str>:
        Provided value to the xavier run, --genome argument
        If value is not valid or custom reference genome JSON file not readable,
        an exception is raised.
    """
    # Checks for custom built genomes using xavier build
    if user_option.endswith(".json"):
        # Check file is readable or accessible
        permissions(parser, user_option, os.R_OK)
    # Checks against valid pre-built options
    # TODO: makes this more dynamic in the future to have it check against
    # a list of genomes (files) in config/genomes/*.json
    elif not user_option in prebuilt:
        # User did NOT provide a valid choice
        parser.error(
            """provided invalid choice, '{}', to --genome argument!\n
        Choose from one of the following pre-built genome options: \n
        \t{}\n
        or supply a custom reference genome JSON file generated from xavier build.
        """.format(
                user_option, prebuilt
            )
        )

    return user_option
