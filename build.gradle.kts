plugins {
    kotlin("multiplatform") version "1.9.0"
    `maven-publish`
}

group = "com.stevenvdb"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
    mavenLocal()
}

kotlin {
    jvm {
        compilations.all {
            kotlinOptions.jvmTarget = "1.8"
        }
        withJava()
        testRuns["test"].executionTask.configure {
            useJUnitPlatform()
        }
    }
    js {
        browser()
        nodejs()
    }

    val hostOs = System.getProperty("os.name")
    val isMingwX64 = hostOs.startsWith("Windows")
    val nativeTarget = when {
        hostOs == "Mac OS X" -> macosX64("native")
        hostOs == "Linux" -> linuxX64("native")
        isMingwX64 -> mingwX64("native")
        else -> throw GradleException("Host OS is not supported in Kotlin/Native.")
    }

    sourceSets {
        val commonMain by getting {
            dependencies {
                implementation("com.stevenvdb:koptimize:1.0")
                implementation("com.stevenvdb:kroot:1.0")
            }
        }
        val commonTest by getting {
            dependencies {
                implementation(kotlin("test"))
            }
        }
    }
}



publishing {
    publications {
        create<MavenPublication>("maven") {
            groupId = "com.stevenvdb"
            artifactId = "kopenvoronoi"
            version = "1.0-SNAPSHOT"
        }
    }
}