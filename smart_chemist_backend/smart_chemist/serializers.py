"""SmartChemist model serializers for django rest framework"""

from rest_framework import serializers


class SmartChemistSubmitSerializer(serializers.Serializer):  # pylint: disable=abstract-method
    """Smart chemist job submission data"""
    smiles = serializers.CharField(required=False)

    def validate(self, data):  # pylint: disable=arguments-renamed
        """Data validation
        :param data: Job submission data
        :type data: collections.OrderedDict
        :raises serializers.ValidationError: If molecule_str is empty.
        :return: Validated data
        :rtype data: collections.OrderedDict
        """

        if len(data['smiles']) == 0:
            raise serializers.ValidationError(
                'Molecule string is empty.')

        return data

#
# class SmartChemistResponseSerializer(serializers.Serializer):  # pylint: disable=abstract-method
#     """SmartChemist response data"""
#     job_id = serializers.UUIDField(required=True)
#     retrieved_from_cache = serializers.BooleanField(default=False)