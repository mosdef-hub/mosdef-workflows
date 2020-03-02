#!/usr/bin/env python3
import glob
import signac
import os
import shutil

project = signac.init_project('water-models')

for model in glob.glob('files/*.xml'):
    basename = os.path.basename(model)
    model_name = os.path.splitext(basename)[0]
    job = project.open_job({'model_name': model_name})
    job.init()
    shutil.copyfile(model, job.fn('ff.xml'))
    print(model_name, 'initialized in job', job)
