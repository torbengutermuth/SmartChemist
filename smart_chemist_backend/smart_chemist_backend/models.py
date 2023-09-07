"""Base models for objects"""
import uuid

from django.db import models

from .job_handler import Status


class JobModel(models.Model):
    """Abstract base model for job objects"""

    class Meta:
        abstract = True

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    status = models.CharField(max_length=1, choices=Status.choices, default=Status.PENDING)
    error = models.TextField(null=True)
    error_detailed = models.TextField(null=True)
    date_created = models.DateField(auto_now_add=True)
    date_last_accessed = models.DateField(auto_now=True)
