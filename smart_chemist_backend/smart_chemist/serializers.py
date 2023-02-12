"""SmartChemist model serializers for django rest framework"""

from rest_framework import serializers


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
        elif "molecule_file" not in data.keys():
            raise serializers.ValidationError("Neither molecule string nor file given")
        return data
