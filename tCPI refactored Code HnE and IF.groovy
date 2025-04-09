/**
Script to determine the tCPI of an enclosed lesion. Script will detect cell objects within the defined region selected before calculating the median signed distance from lesion border and segmenting the area of the lesion into a central and peripheral median. The quotient of the peripheral area and the total area is defined as the total Central Preference Index (tCPI) by Sawyer et al. (2023).

Code to create generate the infiltration annotation is inspired from code produced on the Image.sc Forum, particularly from Pete Bankhead. Regions of code where adapted script appears have been indicated in a comment block.

The script requires a lesion annotation with the "Lesion" class and an overall tissue annotation with class "Tissue".

Script will by default clear all detections when run. Ensure no detections are needed before running the script.

To preserve annotations change their class temporarily from "Lesion" so that they are not included in calculation of another lesion.
**/

import qupath.ext.stardist.StarDist2D
import qupath.lib.scripting.QP

import org.locationtech.jts.geom.Geometry
import qupath.lib.common.GeneralTools
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.ROIs

import java.awt.Rectangle
import java.awt.geom.Area
import org.locationtech.jts.precision.GeometryPrecisionReducer
import org.locationtech.jts.geom.PrecisionModel
import qupath.lib.gui.measure.ObservableMeasurementTableData




// Do you want to keep inner annotation at end of script? Y/N
keepInnerAnnotation = "N"
// Is it Brightfield or IF
scanType = "IF"
// What are the exact names of the marker's used?
def MarkerNames = ["CD20","CD8","CD4","MPO","CD68","Ki67"];
// Do you want ImmCPI (Y/N)
def IncludeImmCPI= "Y"

// Define ImmCPI Calculation
def ImmCpiCalc(ImmuneCellName, anno) {
    
    def ImmuneCellDistance = [];
    def RegularCellDistance = [];

    getDetectionObjects().each{
        if (it.getPathClass()!=null && it.getPathClass().toString().contains(ImmuneCellName)) {
            ImmuneCellDistance.push(it.measurements.get("Signed distance to annotation with Lesion µm"))
        }
    }
    
    if(ImmuneCellDistance.size()==0) {
       ImmuneCellDistance.push(0) 
    }

    getDetectionObjects().each{
        RegularCellDistance.push(it.measurements.get("Signed distance to annotation with Lesion µm"))
    }


    double ImmuneMeanDistance = ImmuneCellDistance.average().toDouble()
    double RegularMeanDistance = RegularCellDistance.average().toDouble()

    double ImmCpi = ImmuneMeanDistance/RegularMeanDistance

    def GranulomaObject = anno
    GranulomaObject.measurements.put("ImmCPI("+ImmuneCellName+")",ImmCpi)
    print("ImmuneCPI for "+ImmuneCellName+" is: "+ImmCpi)
}
// Define median function not in groovy by standard
def median(data) {
    def copy = data.toSorted()
    def middle = data.size().intdiv(2)
    data.size() %2 ? copy[middle] : (copy[middle-1] + copy[middle])/2
}


// Clear Detections
clearDetections()


selectObjectsByClassification("Lesion")

// Replace the following code block with settings for cell detection determined for particular tissue or with alternative models such as StarDist


if(scanType=="Brightfield"){
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield":"Hematoxylin OD","requestedPixelSizeMicrons":0.4994,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.4,"minAreaMicrons":5.0,"maxAreaMicrons":400.0,"threshold":0.35,"maxBackground":2.0,"watershedPostProcess":true,"cellExpansionMicrons":0.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')
}

if(scanType=="IF"){
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"DAPI","requestedPixelSizeMicrons":0.50,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.5,"minAreaMicrons":8.0,"maxAreaMicrons":500.0,"threshold":8,"watershedPostProcess":true,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}');
if(IncludeImmCPI=="Y"){
    runObjectClassifier("95.2.AllMarkers");
}
}

//Measure distance from individual cell to lesion border
detectionToAnnotationDistancesSigned(true)


def signedDistances = [];
getDetectionObjects().each {
    signedDistances.push(it.measurements.get("Signed distance to annotation with Lesion µm"))
}

medCellDistance = Math.abs(median(signedDistances)).toDouble()
print("Median distance of cells to lesion border is "+medCellDistance.toString()+" µm")


/**
Following code block adapted from Pete Bankhead, Mike Nelson, and MicroscopyRA’s script at https://gist.github.com/Svidro/5829ba53f927e79bb6e370a6a6747cfd as Tumor invasion areas 0.3.0.groovy.
**/

PrecisionModel PM = new PrecisionModel(PrecisionModel.FIXED)

//Establish base env Variables
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
if (!cal.hasPixelSizeMicrons()){
    print 'No pixel size calibration detected'
    return
}

double expandPixels = medCellDistance / cal.getAveragedPixelSizeMicrons()
def initLesion = getAnnotationObjects().find{it.getPathClass() == getPathClass("Lesion")}
def lesionGeom = getAnnotationObjects().find{it.getPathClass()==getPathClass("Lesion")}.getROI().getGeometry()
def plane = ImagePlane.getDefaultPlane()
def tissueGeom = getAnnotationObjects().find{it.getPathClass() == getPathClass("Tissue")}.getROI().getGeometry()

/**
End code block
**/

// Ensure lesion geometry is bound by tissue annotation. This covers cases where the lesion occurs on the edge of the tissue

lesionIntersectGeom = tissueGeom.intersection(lesionGeom)
lesionROIClean = GeometryTools.geometryToROI(lesionIntersectGeom, plane)

lesionIntersect = PathObjects.createAnnotationObject(lesionROIClean, getPathClass("Lesion"))
lesionIntersect.setName("Intersected Lesion")

generatedAnnotations = []
generatedAnnotations << lesionIntersect



/**
Following code block adapted from Pete Bankhead, Mike Nelson, and MicroscopyRA’s script at https://gist.github.com/Svidro/5829ba53f927e79bb6e370a6a6747cfd as Tumor invasion areas 0.3.0.groovy.
**/


// Get the central area
def geomCentral = lesionIntersectGeom.buffer(-expandPixels)
geomCentral = geomCentral.intersection(tissueGeom)
def roiCentral = GeometryTools.geometryToROI(geomCentral, plane)
def annotationCentral = PathObjects.createAnnotationObject(roiCentral)
annotationCentral.setName("Center")

// Get the inner margin area
def geomInner = lesionIntersectGeom
geomInner = geomInner.difference(geomCentral)
geomInner = geomInner.intersection(tissueGeom)
def roiInner = GeometryTools.geometryToROI(geomInner, plane)
def annotationInner = PathObjects.createAnnotationObject(roiInner)
annotationInner.setName("Inner margin")
periph = getPathClass("Periphery")
annotationInner.setPathClass(periph)


addObjects(annotationInner)
/**
End code block
**/

def regions = []
regions << initLesion
regions << annotationInner

def area = "Area µm^2"

def ob = new ObservableMeasurementTableData()
ob.setImageData(imageData,regions)

double lesionTotalArea = ob.getStringValue(initLesion, area).toDouble()
double peripheralArea = ob.getStringValue(annotationInner, area).toDouble()
double tCPI = peripheralArea/lesionTotalArea

initLesion.measurements.put("tCPI",tCPI)
print("tCPI is "+tCPI)


if(scanType=="IF" && IncludeImmCPI=="Y") {
    for (marker in MarkerNames) { 
    ImmCpiCalc(marker,initLesion)
    }
}

if(keepInnerAnnotation=="N"){
    removeObject(annotationInner, true)
}