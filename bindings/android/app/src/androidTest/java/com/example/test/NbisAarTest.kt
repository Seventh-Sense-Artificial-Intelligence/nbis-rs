package com.example.test

import androidx.test.platform.app.InstrumentationRegistry
import androidx.test.ext.junit.runners.AndroidJUnit4
import android.content.Context
import androidx.test.core.app.ApplicationProvider
import java.io.InputStream
import org.junit.Test
import org.junit.runner.RunWith
import org.junit.Assert.*

// Import your class from the AAR
import ai.seventhsense.sdk.nbis.NbisExtractor
import ai.seventhsense.sdk.nbis.NbisExtractorSettings
import ai.seventhsense.sdk.nbis.Minutiae
import ai.seventhsense.sdk.nbis.newNbisExtractor
/**
 * Instrumented test, which will execute on an Android device.
 *
 * See [testing documentation](http://d.android.com/tools/testing).

@RunWith(AndroidJUnit4::class)
class ExampleInstrumentedTest {
    @Test
    fun useAppContext() {
        // Context of the app under test.
        val appContext = InstrumentationRegistry.getInstrumentation().targetContext
        assertEquals("ai.seventhsense.sdk.nbis.test", appContext.packageName)
    }
}
 */

@RunWith(AndroidJUnit4::class)
class MinutiaeInstrumentationTest {

    // Helper to load a fingerprint image from test assets (you should place test image in `app/src/androidTest/assets/`)
    private fun loadFingerprintBytes(fileName: String): ByteArray {
        val context = ApplicationProvider.getApplicationContext<Context>()
        val inputStream: InputStream = context.assets.open(fileName)
        return inputStream.readBytes()
    }

    @Test
    fun testExtractMinutiaeAndCompare() {
        // Load sample image(s)
        val image1 = loadFingerprintBytes("p1_1.png")
        val image2 = loadFingerprintBytes("p1_2.png")

        val settings = NbisExtractorSettings(
            minQuality = 0.2,         // Example value
            getCenter = true,         // Whether to extract core point
            checkFingerprint = true,  // Whether to verify valid fingerprint
            computeNfiq2 = true,      // Whether to compute NFIQ2 quality score
            ppi = null                // Pixels per inch of sensor (standard is 500)
        )

        val extractor: NbisExtractor = newNbisExtractor(settings)

        val minutiae1: Minutiae = extractor.extractMinutiae(image1)
        val minutiae2: Minutiae = extractor.extractMinutiae(image2)

        assertTrue("Expected a non-zero quality score", minutiae1.quality().score > 0u)
        assertTrue("Expected a non-zero quality score", minutiae2.quality().score > 0u)

        assertNotNull(minutiae1)
        assertNotNull(minutiae2)

        // Compare the extracted minutiae
        val similarityScore = minutiae1.compare(minutiae2)
        println("Similarity score: $similarityScore")

        // Check if match is likely
        assertTrue("Expected some similarity", similarityScore >= 0)
    }
}