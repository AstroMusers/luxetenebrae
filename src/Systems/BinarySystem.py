# from src.Properties.SystemParameters import SystemParameters
# from src.Properties.CommonEnvelopes import CommonEnvelopes
# from src.Properties.Supernovae import Supernovae
# from src.Properties.MassTransfers import MassTransfers
# from src.Properties.DoubleCompactObjects import DoubleCompactObjects

class BinarySystem:
    """
    A class to represent a binary star system in the luxetenebrae package.
    Public attributes: data, properties
    Internal (prefixed with _): _available_properties, _requested_properties, _build_properties, _set_property
    """

    def __init__(self, data, requested_properties=None):
        self.data = data
        self._available_properties = data.keys()
        self._requested_properties = requested_properties or []
        self.properties = self._build_properties()

    def _build_properties(self):
        """Build all requested properties at initialization."""
        props = {}
        for key in self._requested_properties:
            if key in self._available_properties:
                props[key] = self._set_property(key)
            else:
                print(f"Property {key} not available in data.")
        return props

    def get_property(self, key):
        """Retrieve a single property (compute if not cached)."""
        if key in self.properties:
            return self.properties[key]
        if key not in self._available_properties:
            raise KeyError(f"Property {key} not available in data.")
        value = self._set_property(key)
        self.properties[key] = value
        return value

    def _set_property(self, key):
        """Set property object/value."""
        if key == "SP":
            return SystemParameters(self.data['BSE_System_Parameters'])
        elif key == "CE":
            return CommonEnvelopes(self.data['BSE_Common_Envelopes'])
        elif key == "SN":
            return Supernovae(self.data['BSE_Supernovae'])
        elif key == "MT":
            return MassTransfers(self.data['BSE_RLOF'])
        elif key == "DCO":
            return DoubleCompactObjects(self.data['BSE_Double_Compact_Objects'])