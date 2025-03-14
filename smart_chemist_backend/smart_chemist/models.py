from django.db import models

from smart_chemist_backend.models import JobModel


class PatternMatchingInputModel(models.Model):
    """PatternMatching input model

    Contains information about the input to a pattern matching task.
    """
    # link to parent job
    parent_pattern_matching_job = models.OneToOneField('PatternMatchingJob', on_delete=models.CASCADE)
    # input string for molecule data. Can be file content from .sdf,.smi or just comma separated SMILES.
    input_string = models.TextField()
    # input format specification. Can be file extension.
    input_format = models.CharField(max_length=16)
    # max. number of molecule inputs to consider
    input_max_nof_molecules_allowed = models.IntegerField(default=100)


class PatternMatchingOutputModel(models.Model):
    """PatternMatching output model

    Contains information about the output to a pattern matching task.
    """
    # link to parent job
    parent_pattern_matching_job = models.OneToOneField('PatternMatchingJob', on_delete=models.CASCADE)
    # will hold results of the job in JSON format
    output_json = models.JSONField()


class PatternMatchingJob(JobModel):
    """Django Model for pattern matching job objects"""

    # model that holds the input data of the job
    input_info = models.OneToOneField(PatternMatchingInputModel, on_delete=models.CASCADE, null=True)
    # model that will hold the output data of the job when the job is finished
    output_info = models.OneToOneField(PatternMatchingOutputModel, on_delete=models.CASCADE, null=True)


class AnnotatedPattern(models.Model):
    """Model to store SMARTS expressions with associated annotations

    For example a pattern describing a trivial name, functional groups, etc.
    """
    smarts = models.CharField(max_length=500)
    trivial_name = models.CharField(max_length=100)
    group = models.CharField(max_length=256)
    hierarchy = models.CharField(max_length=500, default=None)
    index_file = models.IntegerField(default=0)
    heavy_atoms = models.IntegerField(default=0)
    num_rings = models.IntegerField(default=0)
    n_nitrogens = models.IntegerField(default=0)
    n_sulfur = models.IntegerField(default=0)
    n_oxygen = models.IntegerField(default=0)
    n_carbon = models.IntegerField(default=0)
    n_phosphor = models.IntegerField(default=0)
    n_halogens = models.IntegerField(default=0)
    n_other_atom = models.IntegerField(default=0)


class BugReport(models.Model):
    molecule_name = models.CharField(max_length=500)
    smiles = models.CharField(max_length=500)
    user_message = models.CharField(max_length=500)
    matches = models.JSONField()
    viewed = models.BooleanField(default=False)
    posting_time = models.DateTimeField()
