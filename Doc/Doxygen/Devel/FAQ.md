Frequently Asked Questions    {#FAQPage}
==========================


## What is the color value range?

Colors are stored as R, G, B triplets. With meshes, each color value is a
floating point (single or double precision) in the range [0.0 ; 1.0]. With
point clouds, each color value is a 8 bits unsigned int in the range [0 ;
255].

The filters that apply to both meshes and point clouds may have to convert
between the two different ranges of color value. For example, the [Copy Graph
filter](https://projet.liris.cnrs.fr/mepp/doc/nightly/_filter_copy_graph.html)
filter can copy a mesh to a point cloud (and vice versa) with its vertex
colors. So it must convert colors to suit the range of color values supported
by the target data structure.

To make it easier to write such filters, some convenient functions that do
the conversion have been put into [the dedicated color_conversion.hpp
file](https://github.com/MEPP-team/MEPP2/blob/b2e15df7b47a7d58cca4977bc503519021220385/FEVV/Tools/Math/color_conversion.hpp).
The conversion functions use the syntax:

    OutputType convert_color_value< OutputType >(InputType value)

The conversion is based on the type of the color
value passed as a parameter, and the type of the value returned by the
function. For example, a conversion from type 'uint8_t' to type 'float' will
transform the color value from the range [0 ; 255] to the range [0.0 ; 1.0].
If the input type and the output type match the same range, then the color
value is not changed. In this case, only a cast is applied to produce the
right output type.

Examples:

    // convert color to range [0.0 ; 1.0]
    // if color is already in the target range, values are casted without change
    FPColor color_0_1(FEVV::Math::convert_color_value< float >(color[0]),
                      FEVV::Math::convert_color_value< float >(color[1]),
                      FEVV::Math::convert_color_value< float >(color[2]));

    // convert color to range [0 ; 255]
    // if color is already in the target range, values are casted without change
    UInt8Color color_0_255(FEVV::Math::convert_color_value< uint8_t >(color[0]),
                           FEVV::Math::convert_color_value< uint8_t >(color[1]),
                           FEVV::Math::convert_color_value< uint8_t >(color[2]));

