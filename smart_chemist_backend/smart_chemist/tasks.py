"""smart_chemist celery tasks"""
from celery import shared_task
from smart_chemist_backend.job_handler import execute_job
from .models import PatternMatchingJob
from .smart_chemist import SmartChemist


@shared_task
def pattern_matching_task(job_id):
    """Start job execution

    :param job_id: Database id of the job object to be executed
    :type job_id: uuid
    """
    execute_job(pattern_matching, job_id, PatternMatchingJob, 'PatternMatching')


def pattern_matching(job):
    """Run pattern matching job and store results into database

    :param job: PatternMatchingJob object containing the job data
    :type job: PatternMatchingJob
    """
    SmartChemist.handle_pattern_matching(job)
