"""SmartChemist model serializers for django rest framework"""

from rest_framework import serializers

from smart_chemist.models import PatternMatchingJob, PatternMatchingInputModel, PatternMatchingOutputModel
from smart_chemist_backend.serializers import JobSerializer


class PatternMatchingJobSerializer(JobSerializer):
    """Pattern matching job data"""

    class Meta(JobSerializer.Meta):
        model = PatternMatchingJob
        fields = JobSerializer.Meta.fields + [
            'input_info',
            'output_info',
        ]


class PatternMatchingInputModelSerializer(serializers.ModelSerializer):
    """Pattern matching input info data"""

    class Meta:
        model = PatternMatchingInputModel
        fields = ['id', 'input_max_nof_molecules_allowed', 'input_format', 'input_string',
                  'parent_pattern_matching_job']


class PatternMatchingOutputModelSerializer(serializers.ModelSerializer):
    """Pattern matching output info data"""

    class Meta:
        model = PatternMatchingOutputModel
        fields = ['id', 'output_json']


class SmartChemistSubmitSerializer(serializers.Serializer):  # pylint: disable=abstract-method
    """Smart chemist job submission data"""
    smiles = serializers.CharField(required=False)
    molecule_file = serializers.FileField(required=False)

    def validate(self, data):  # pylint: disable=arguments-renamed
        """Data validation
        :param data: Job submission data
        :type data: collections.OrderedDict
        :raises serializers.ValidationError: If molecule_str is empty.
        :return: Validated data
        :rtype data: collections.OrderedDict
        """

        if "smiles" in data.keys():
            if "molecule_file" in data.keys():
                raise serializers.ValidationError("Both molecule file and smiles given.")
            elif len(data["smiles"]) == 0:
                raise serializers.ValidationError("Molecule string is empty.")
        return data


class BugReportSubmitSerializer(serializers.Serializer):
    smiles = serializers.CharField(required=True)
    bug_report = serializers.CharField(required=True)
    name = serializers.CharField(required=False, default="")
    matches = serializers.CharField(required=False, default="")

    def validate(self, data):
        return data
