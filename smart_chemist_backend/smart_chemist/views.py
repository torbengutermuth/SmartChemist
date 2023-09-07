"""smart chemist views"""
from django.http import JsonResponse
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.parsers import JSONParser, MultiPartParser, FormParser
from rest_framework.response import Response
from drf_spectacular.utils import extend_schema
from rest_framework.viewsets import ReadOnlyModelViewSet

from smart_chemist_backend.job_handler import submit_task
from smart_chemist_backend.serializers import JobResponseSerializer
from .models_utils import make_pattern_matching_job
from .serializers import SmartChemistSubmitSerializer, BugReportSubmitSerializer, PatternMatchingJobSerializer, \
    PatternMatchingInputModelSerializer, PatternMatchingOutputModelSerializer
from .models import BugReport, PatternMatchingJob, PatternMatchingInputModel, PatternMatchingOutputModel
import datetime

from .tasks import pattern_matching_task


class SmartChemistView(APIView):
    """View for running smart chemist"""

    parser_classes = (JSONParser, MultiPartParser, FormParser)

    @extend_schema(
        request=SmartChemistSubmitSerializer,
        # responses=ProteinsPlusJobResponseSerializer
    )
    def post(self, request):
        """Annotate molecule with SmartChemist.

        """

        # validate the input parameters
        serializer = SmartChemistSubmitSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        request_data = serializer.validated_data

        # create a new pattern matching job object
        job = make_pattern_matching_job(request_data, nof_molecules_allowed=100)

        # submit the job
        job_id = submit_task(job, pattern_matching_task)
        serializer = JobResponseSerializer({
            'job_id': job_id,
        })
        return Response(serializer.data, status=status.HTTP_202_ACCEPTED)


class SmartChemistBugReportView(APIView):
    parser_classes = (JSONParser, MultiPartParser, FormParser)

    @extend_schema(
        request=BugReportSubmitSerializer,
        # responses=ProteinsPlusJobResponseSerializer
    )
    def post(self, request):
        serializer = BugReportSubmitSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        request_data = serializer.validated_data
        thank_you_json_dict = {}
        report = BugReport(molecule_name=request_data["name"],
                           smiles=request_data["smiles"],
                           user_message=request_data["bug_report"],
                           matches=request_data["matches"],
                           posting_time=datetime.datetime.now(),
                           )
        report.save()
        return JsonResponse(thank_you_json_dict, status=status.HTTP_202_ACCEPTED, safe=False)


class PatternMatchingJobViewSet(ReadOnlyModelViewSet):  # pylint: disable=too-many-ancestors
    """Retrieve specific or list all PatternMatching job"""
    queryset = PatternMatchingJob.objects.all()
    serializer_class = PatternMatchingJobSerializer


class PatternMatchingInputModelViewSet(ReadOnlyModelViewSet):  # pylint: disable=too-many-ancestors
    """Retrieve specific or list all PatternMatchingInputModel info objects"""
    queryset = PatternMatchingInputModel.objects.all()
    serializer_class = PatternMatchingInputModelSerializer


class PatternMatchingOutputModelViewSet(ReadOnlyModelViewSet):  # pylint: disable=too-many-ancestors
    """Retrieve specific or list all PatternMatching result info objects"""
    queryset = PatternMatchingOutputModel.objects.all()
    serializer_class = PatternMatchingOutputModelSerializer
