"""Base serializers for ProteinsPlus objects"""
from rest_framework import serializers
from smart_chemist_backend.job_handler import StatusField


class JobSerializer(serializers.ModelSerializer):
    """Base serializer for job models"""
    status = StatusField()

    class Meta:
        fields = ['id', 'status', 'date_created', 'date_last_accessed', 'error']


class JobSubmitSerializer(serializers.Serializer):  # pylint: disable=abstract-method
    """Base serializer for job submit data"""
    pass


class JobResponseSerializer(serializers.Serializer):  # pylint: disable=abstract-method
    """Job submission response data"""
    job_id = serializers.UUIDField(required=True)
