#pragma once

#include <LAS/Data/nVector.h>


namespace Symmetry {
    /// <summary>
    /// Namespace of colors used in visualization.
    /// </summary>
    namespace Color {
        const LAS::Data::Vector4d black =                   LAS::Data::Vector4d(0.0, 0.0, 0.0, 1.0);
        const LAS::Data::Vector4d gray =                    LAS::Data::Vector4d(0.4, 0.4, 0.4, 1.0);
        const LAS::Data::Vector4d red =                     LAS::Data::Vector4d(1.0, 0.0, 0.0, 1.0);
        const LAS::Data::Vector4d blue =                    LAS::Data::Vector4d(0.0, 0.0, 1.0, 1.0);
        const LAS::Data::Vector4d green =                   LAS::Data::Vector4d(0.0, 0.7, 0.0, 1.0);
        const LAS::Data::Vector4d magenta =                 LAS::Data::Vector4d(1.0, 0.0, 1.0, 1.0);
        const LAS::Data::Vector4d voxelLight =              LAS::Data::Vector4d(0.3, 0.3, 0.3, 0.4);
        const LAS::Data::Vector4d voxelDark =               LAS::Data::Vector4d(0.2, 0.2, 0.2, 0.4);
        const LAS::Data::Vector4d interestingVoxelLight =   LAS::Data::Vector4d(1.0, 0.8, 0.0, 0.4);
        const LAS::Data::Vector4d interestingVoxelDark =    LAS::Data::Vector4d(1.0, 0.6, 0.0, 0.4);
        const LAS::Data::Vector4d symmetryVoxelLightBlue =  LAS::Data::Vector4d(0.0, 0.0, 1.0, 0.4);
        const LAS::Data::Vector4d symmetryVoxelLightRed =   LAS::Data::Vector4d(1.0, 0.0, 0.0, 0.4);
        const LAS::Data::Vector4d symmetryVoxelLightGreen = LAS::Data::Vector4d(0.2, 0.9, 0.4, 0.4);
        const LAS::Data::Vector4d symmetryVoxelPurple =     LAS::Data::Vector4d(0.6, 0.2, 0.9, 0.4);
        const LAS::Data::Vector4d lineSegment =             LAS::Data::Vector4d(0.0, 0.5, 0.5, 1.0);
        const LAS::Data::Vector4d symmetryAxis =            LAS::Data::Vector4d(0.0, 0.8, 0.5, 1.0);
        const LAS::Data::Vector4d lightGreen =              LAS::Data::Vector4d(0.2, 0.9, 0.4, 1.0);
    };
};