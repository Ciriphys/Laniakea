-- TODO : Change Workspace Name, Choose 3-letter signature for macros based on name.
-- TODO : Add PCH
workspace "Laniakea" 
    architecture "x64"
    configurations { "Debug", "Release" }
    startproject "Sandbox"

    outDir = "%{cfg.buildcfg}_%{cfg.system}_%{cfg.architecture}"

    project "Sandbox"
        location "Build/Sandbox"
        kind "ConsoleApp"
        language "C++"
        cppdialect "C++17"

        targetdir   ( "Binaries/" .. outDir .. "/%{prj.name}" )
        objdir      ( "Binaries/Objects/" .. outDir .. "/%{prj.name}" )

        files { "Build/Sandbox/Source/**.cpp" }
        includedirs { "Build/Core", "Build/Core/Include" }

        links { "Library" }

        filter "system:Windows"
            staticruntime "On"
            systemversion "latest"
            system "windows"
            defines { "SBX_WIN", "LNK_WIN", "_CRT_SECURE_NO_WARNINGS" }

        filter "system:Macosx"
            system "macosx"
            defines { "SBX_MACOS", "LNK_MACOS" }

        filter "system:Linux"
            pic "On"
            system "Linux"
            defines { "SBX_LINUX", "LNK_LINUX" }
            buildoptions { "-Wno-unused-result" }

        filter { "configurations:Debug" }
            defines { "SBX_DEBUG", "DEBUG" }
            optimize "Debug"
            symbols "On"

        filter { "configurations:Release" }
            defines { "SBX_RELEASE", "NDEBUG" }
            optimize "Full"
            symbols "Off"

    project "Library"
        location "Build/Core"
        kind "SharedLib"
        language "C++"
        cppdialect "C++17"

        targetdir   ( "Binaries/" .. outDir .. "/%{prj.name}" )
        objdir      ( "Binaries/Objects/" .. outDir .. "/%{prj.name}" )

        files { "Build/Core/Source/**.cpp", "Build/Core/Include/**.h", "Build/Core/lnkpch.h"}
        includedirs { "Build/Core", "Build/Core/Include" }

        pchheader "lnkpch.h" 
        pchsource "Build/Core/Source/lnkpch.cpp" 

        defines { "LNK_BUILD_DLL" } 

        filter "system:Windows"
            staticruntime "On"
            systemversion "latest"
            system "windows"
            defines { "LNK_WIN", "_CRT_SECURE_NO_WARNINGS" } 

            postbuildcommands
            {
                "{COPY} ../../Binaries/" .. outDir .. "/Library/*.dll ../../Binaries/" .. outDir .. "/Sandbox",
                "{COPY} ../../Binaries/" .. outDir .. "/Library/*.lib ../../Binaries/" .. outDir .. "/Sandbox"
            }

        filter "system:Macosx"
            system "macosx"
            defines { "LNK_MACOS" } 

            prelinkcommands
            {
               ("mkdir -p ../../Binaries/" .. outDir .. "/Sandbox")
            }
    
            postbuildcommands
            {
                ("{COPY} %{cfg.buildtarget.relpath} ../../Binaries/" .. outDir .. "/Sandbox")
            }

        filter "system:Linux"
            pic "On"
            system "Linux"
            defines { "LNK_LINUX" }
            buildoptions { "-Wno-unused-result" }

            prelinkcommands
            {
               ("mkdir -p ../../Binaries/" .. outDir .. "/Sandbox")
            }
    
            postbuildcommands
            {
                ("{COPY} %{cfg.buildtarget.relpath} ../../Binaries/" .. outDir .. "/Sandbox")
            }

        filter { "configurations:Debug" }
            defines { "LNK_DEBUG", "DEBUG" }
            optimize "Debug"
            symbols "On"

        filter { "configurations:Release" }
            defines { "LNK_RELEASE", "NDEBUG" }
            optimize "Full"
            symbols "Off"
