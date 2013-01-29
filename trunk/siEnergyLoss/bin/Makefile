EVENTS = 10000

all : .realLife .processPixel .processStrip

.realLife : realLife
	$< 0 &
	$< 1
	$< 2 &
	$< 3
	touch $@

.processPixel : clusterFit .generatePixel
	$< -process  $(EVENTS) -pixel
	touch $@

.processStrip : clusterFit .generateStrip
	$< -process  $(EVENTS) -strip
	touch $@

.generatePixel : generate
	$< -generate $(EVENTS) -pixel
	touch $@

.generateStrip : generate
	$< -generate $(EVENTS) -strip
	touch $@
