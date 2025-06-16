#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import singlecellmultiomics
import argparse
from colorama import Fore, Style
import os
from shutil import copyfile
import importlib.resources

def get_workflow_list(include_hidden=False):
    workflows_dir = importlib.resources.files('singlecellmultiomics').joinpath('snakemake_workflows')
    return [workflow.name
            for workflow in workflows_dir.iterdir()
            if workflow.is_dir() and (include_hidden or not workflow.name.startswith('_'))]


def deploy_workflow_files(name, clean=False, directory='./'):
    """ Deploy workflow configuration

    Args:
        name (str) : Alias of the workflow files to deploy
        directory (str) : deploy to this directory
    """

    workflow_list = get_workflow_list(True)
    if not name in workflow_list:
        print(workflow_list)
        raise ValueError(f"Unkown workflow {name}")


    base_path = importlib.resources.files('singlecellmultiomics').joinpath(f'snakemake_workflows/{name}')
    for file_to_copy in base_path.iterdir():
        # Files starting with a _ are never copied
        if file_to_copy.name.startswith('_'):
            continue
        print(Fore.GREEN + f'Creating {file_to_copy.name}' + Style.RESET_ALL, end="\t")

        target = os.path.join(directory, file_to_copy.name)
        if os.path.exists(target) and not clean:
            print(Fore.RED + f'The file {target} already exists! Skipping!' + Style.RESET_ALL)
        else:
            with importlib.resources.path('singlecellmultiomics', f'snakemake_workflows/{name}/{file_to_copy.name}') as source:
                copyfile(source, target)
            print(Fore.GREEN + f'[ok]' + Style.RESET_ALL)

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain single cell multiomics workflow configuration files, deploys configuration in current directory')
    argparser.add_argument('alias', type=str, help= 'nla, chic', nargs='?')
    argparser.add_argument('--clean', action='store_true',help='clear all existing files')
    args = argparser.parse_args()

    if args.alias is None:
        # Show available :
        print(Style.BRIGHT + "The following workflows are available:" +Style.RESET_ALL)
        print('\n'.join(get_workflow_list()))

    else:
        deploy_workflow_files(args.alias, clean=args.clean)
        #deploy_workflow_files('_general', clean=args.clean)

        print(f"""\nEdit config.json and then run either;\n
        On local computer:\n
        {Fore.BLUE}snakemake{Style.RESET_ALL}\n
        On SGE cluster:\n
        {Fore.BLUE}snakemake --cluster sge_wrapper.py  --jobs 20 --restart-times 3{Style.RESET_ALL}\n
        On SLURM cluster:\n
        {Fore.BLUE}snakemake --cluster slurm_wrapper.py  --jobs 20 --restart-times 3{Style.RESET_ALL}\n

        To keep intermediate files, add {Fore.BLUE}--nt{Style.RESET_ALL}
        """)
