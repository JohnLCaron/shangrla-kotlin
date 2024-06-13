plugins {
    kotlin("jvm") version "1.9.23"
    alias(libs.plugins.ktor)
    alias(libs.plugins.serialization)
}

group = "org.cryptobiotic"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
}

dependencies {
    implementation(libs.bull.result)
    implementation(libs.bundles.xmlutil)
    implementation(libs.ktor.serialization.kotlinx.json.jvm )
    implementation(libs.bundles.logging)
    testImplementation(kotlin("test"))
}

tasks.test {
    useJUnitPlatform()
}
kotlin {
    jvmToolchain(21)
}
