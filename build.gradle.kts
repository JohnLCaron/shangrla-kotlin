plugins {
    kotlin("jvm") version "1.9.23"
    alias(libs.plugins.ktor)
    alias(libs.plugins.serialization)
}

group = "org.cryptobiotic"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
    maven("https://packages.jetbrains.team/maven/p/kds/kotlin-ds-maven")
}

dependencies {
    implementation(libs.bull.result)
    implementation(libs.bundles.xmlutil)
    implementation(libs.ktor.serialization.kotlinx.json.jvm )
    implementation(libs.kotlinx.cli)
    implementation(libs.bundles.logging)

    implementation("org.jetbrains.kotlinx:kandy-lets-plot:0.6.0")
    implementation("org.jetbrains.kotlinx:kotlin-statistics-jvm:0.2.1")

    testImplementation(kotlin("test"))
}

tasks.test {
    useJUnitPlatform()
}
kotlin {
    jvmToolchain(21)
}

application {
    mainClass.set("org.cryptobiotic.rla.RunRlaTest")
}
