# trace generated using paraview version 5.6.2
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
output000000 = FindSource('output000000*')

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=output000000)
tableToPoints1.XColumn = '1'
tableToPoints1.YColumn = '2'
tableToPoints1.ZColumn = '3'

# get active view
spreadSheetView1 = GetActiveViewOrCreate('SpreadSheetView')
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(output000000, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [233, 551]

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display_1 = Show(tableToPoints1, renderView1)

# trace defaults for the display properties.
tableToPoints1Display_1.Representation = 'Surface'
tableToPoints1Display_1.ColorArrayName = [None, '']
tableToPoints1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display_1.ScaleFactor = 19.5
tableToPoints1Display_1.GlyphType = 'Arrow'
tableToPoints1Display_1.GaussianRadius = 0.975
tableToPoints1Display_1.SetScaleArray = [None, '']
tableToPoints1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.OpacityArray = [None, '']
tableToPoints1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display_1.SelectionCellLabelFontFile = ''
tableToPoints1Display_1.SelectionPointLabelFontFile = ''
tableToPoints1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tableToPoints1Display_1.DataAxesGrid.XTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.YTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.ZTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.XLabelFontFile = ''
tableToPoints1Display_1.DataAxesGrid.YLabelFontFile = ''
tableToPoints1Display_1.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tableToPoints1Display_1.PolarAxes.PolarAxisTitleFontFile = ''
tableToPoints1Display_1.PolarAxes.PolarAxisLabelFontFile = ''
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextFontFile = ''
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Glyph'
glyph1 = Glyph(Input=tableToPoints1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = [None, '0']
glyph1.ScaleFactor = 19.5
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.ScaleArray = ['POINTS', '0']

# show data in view
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.ScaleFactor = 2228.2000000000003
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GaussianRadius = 111.41
glyph1Display.SetScaleArray = [None, '']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = [None, '']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# Properties modified on glyph1
glyph1.ScaleArray = ['POINTS', '4']

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 0.007818765036086608

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip1 = Clip(Input=glyph1)
clip1.ClipType = 'Plane'
clip1.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [9.75, -1.0, 3.0]

# Properties modified on clip1
clip1.Scalars = ['POINTS', '']

# show data in view
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.OSPRayScaleArray = '0'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = '0'
clip1Display.ScaleFactor = 19.78999938964844
clip1Display.SelectScaleArray = '0'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = '0'
clip1Display.GaussianRadius = 0.9894999694824219
clip1Display.SetScaleArray = ['POINTS', '0']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', '0']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 12.34241667806186

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display.DataAxesGrid.XTitleFontFile = ''
clip1Display.DataAxesGrid.YTitleFontFile = ''
clip1Display.DataAxesGrid.ZTitleFontFile = ''
clip1Display.DataAxesGrid.XLabelFontFile = ''
clip1Display.DataAxesGrid.YLabelFontFile = ''
clip1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(glyph1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(tableToPoints1, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(clip1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, ('POINTS', '7'))

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'a7'
a7LUT = GetColorTransferFunction('a7')
a7LUT.RGBPoints = [3.0, 0.231373, 0.298039, 0.752941, 3.000244140625, 0.865003, 0.865003, 0.865003, 3.00048828125, 0.705882, 0.0156863, 0.14902]
a7LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'a7'
a7PWF = GetOpacityTransferFunction('a7')
a7PWF.Points = [3.0, 0.0, 0.5, 0.0, 3.00048828125, 1.0, 0.5, 0.0]
a7PWF.ScalarRangeInitialized = 1

# Properties modified on a7LUT
a7LUT.InterpretValuesAsCategories = 1
a7LUT.AnnotationsInitialized = 1

# Properties modified on a7LUT
a7LUT.Annotations = ['3', '3']
a7LUT.IndexedColors = [1.0, 1.0, 1.0]
a7LUT.IndexedOpacities = [1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['3', 'premitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['3', 'premitotic', '', '']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['3', 'premitotic', '1', '']

# Properties modified on a7LUT
a7LUT.Annotations = ['3', 'premitotic', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '', '', '1', 'postmitotic']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', '', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '', '', '1', 'postmitotic']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', '', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '', '', '1', 'postmitotic']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', '', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic ', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '', '', '1', 'postmitotic']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', '', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', 'necrotic lysing', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', 'necrotic lysis', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', 'necrotic lysis', '', '', '1', 'postmitotic']
a7LUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]
a7LUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', 'necrotic lysis', '103', '', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.Annotations = ['0', 'premitotic', '3', 'Ki67 negative', '100', 'apoptotic', '101', 'necrotic swelling', '102', 'necrotic lysis', '103', 'necrotic', '1', 'postmitotic']

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 0.0, 0.0, 0.3333333333333333, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 1.0, 0.0, 0.3333333333333333, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 1.0, 1.0, 0.0, 0.0, 0.3176470588235294, 0.3176470588235294, 0.3176470588235294, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.0]

# Properties modified on a7LUT
a7LUT.IndexedColors = [0.0, 0.6666666666666666, 0.4980392156862745, 0.0, 0.3333333333333333, 1.0, 1.0, 0.0, 0.0, 0.3176470588235294, 0.3176470588235294, 0.3176470588235294, 0.5803921568627451, 0.5803921568627451, 0.5803921568627451, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.0]

# get layout
layout1 = GetLayout()

# resize frame
layout1.SetSplitFraction(0, 1.0)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [11.207098553248258, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [11.207098553248258, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [15.642223400734785, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [15.642223400734785, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [23.97358258817356, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [23.97358258817356, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [32.304941775612335, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [32.304941775612335, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [38.67222649613955, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [38.67222649613955, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [37.722103856554696, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [37.722103856554696, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [32.971490658630024, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [32.971490658630024, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [31.989453425174098, -1.0, 3.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [31.989453425174098, -1.0, 3.0]

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [500.95281098001414, 153.5687667766207, -650.4656551555004]
renderView1.CameraFocalPoint = [-0.9999999999999994, -0.4999999999999989, 2.9999999999999973]
renderView1.CameraViewUp = [-0.5476323895348909, 0.8041845465956594, -0.2310540649913766]
renderView1.CameraParallelScale = 395.47643624550597

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).