"""Common classes and functions for all apps"""
import logging
import traceback
from rest_framework import serializers
# from drf_spectacular.utils import extend_schema_field, OpenApiTypes

logger = logging.getLogger(__name__)


class Status:  # pylint: disable=too-few-public-methods
    """Class wrapping a status enum"""
    PENDING = 'p'
    RUNNING = 'r'
    SUCCESS = 's'
    FAILURE = 'f'

    DETAILED = {PENDING: 'pending', RUNNING: 'running', SUCCESS: 'success', FAILURE: 'failure'}

    choices = [
        (PENDING, DETAILED[PENDING]),
        (RUNNING, DETAILED[RUNNING]),
        (SUCCESS, DETAILED[SUCCESS]),
        (FAILURE, DETAILED[FAILURE]),
    ]

    @staticmethod
    def to_string(status):
        """Convert status to a readable string

        :param status: status abbreviation
        :type status: str
        :return: detailed status name
        :rtype: str
        """
        return Status.DETAILED[status]


# @extend_schema_field(OpenApiTypes.STR)
class StatusField(serializers.Field):  # pylint: disable=abstract-method
    """Custom serializer field for status"""

    def to_representation(self, value):
        return Status.to_string(value)


def execute_job(task, job_id, job_type, tool_name):
    """Execute Job as a celery task following a centralized workflow

    :param task: task to be executed
    :type task: function
    :param job_id: Id of the job object
    :type job_id: int
    :param job_type: Database table in which to find the job object
    :type job_type: django.db.models.Model
    :raises error: If an error occurs during job execution
    """

    logger.info('Started task. Executing %s on %s with id %s.', task, job_type, job_id)
    job = job_type.objects.get(id=job_id)
    try:
        job.status = Status.RUNNING
        job.save()

        task(job)
    except Exception as error:
        job.status = Status.FAILURE
        job.error = f'An error occurred during the execution of {tool_name}.'
        job.error_detailed = traceback.format_exc()
        job.save()

        logger.error(
            'Error occurred during execution of task %s on %s with id %s.\n'
            'Error: \n\t%s', task, job_type, job_id, traceback.format_exc()
        )
        raise error
    else:
        job.status = Status.SUCCESS
        job.save()
        logger.info('Successfully finished executing %s on %s with id %s.', task, job_type, job_id)


def submit_task(job, task, immediate=False):
    """Retrieve an existing job from cache or start its execution

    :param job: Job to be executed or retrieved from cache
    :type job: ProteinsPlusJob
    :param task: The task function to be called if no caching occurrs
    :type task: celery.local.celery.local (celery shared task)
    :param immediate: Whether to process the task immediately
    :type immediate: bool
    :return: the job id and whether the job was retrieved from cache
    :rtype: tuple(job_id, retrieved)
    """

    job.save()
    task(job.id) if immediate else task.delay(job.id)

    return job.id
