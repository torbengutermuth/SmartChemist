from django.db import models


class AnnotatedPattern(models.Model):
    """Model to store SMARTS expressions with associated annotations

    For example a pattern describing a trivial name, functional groups, etc.
    """
    smarts = models.CharField(max_length=500)
    trivial_name = models.CharField(max_length=100)
    group = models.CharField(max_length=256)

#
# class Molecule(models.Model):
#     chemical_molecule = models.CharField()
#     file_format = models