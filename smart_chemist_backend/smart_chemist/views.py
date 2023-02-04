"""smart chemist views"""
from django.http import JsonResponse
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.viewsets import ReadOnlyModelViewSet
from rest_framework.parsers import JSONParser, MultiPartParser, FormParser
from rest_framework.response import Response
from drf_spectacular.utils import extend_schema

from .serializers import SmartChemistSubmitSerializer
from .smart_chemist import SmartChemist


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
        serializer = SmartChemistSubmitSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        request_data = serializer.validated_data

        result_json_dict = None
        if 'smiles' in request_data:
            result_json_dict = SmartChemist.handle_string_input(request_data["smiles"])
        elif 'molecule_file' in request_data:
            result_json_dict = SmartChemist.handle_file_input(request_data["molecule_file"])

        return JsonResponse(result_json_dict, status=status.HTTP_202_ACCEPTED, safe=False)
