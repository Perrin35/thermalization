OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0639609) q[0];
sx q[0];
rz(-1.504717) q[0];
sx q[0];
rz(1.328701) q[0];
rz(0.98923367) q[2];
sx q[2];
rz(-2.6510112) q[2];
sx q[2];
rz(1.7218334) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0393385) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(-1.1182055) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47031109) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(0.51629984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(0.84709644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0711489) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(1.1308934) q[0];
rz(0.77913021) q[2];
sx q[2];
rz(-1.3574294) q[2];
sx q[2];
rz(-1.9386292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-0.54547711) q[1];
rz(-pi) q[2];
rz(0.91244016) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(-0.57953366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2505328) q[0];
sx q[0];
rz(-0.30936229) q[0];
sx q[0];
rz(0.022457794) q[0];
rz(1.8965365) q[2];
sx q[2];
rz(-2.2176761) q[2];
sx q[2];
rz(-1.5145472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4434102) q[1];
sx q[1];
rz(-1.3122307) q[1];
sx q[1];
rz(-2.3727388) q[1];
rz(1.2445883) q[3];
sx q[3];
rz(-2.2766114) q[3];
sx q[3];
rz(2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(0.8262659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4429312) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(-2.4311964) q[0];
x q[1];
rz(-2.4414908) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(-2.5259279) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(-3.0613042) q[1];
rz(-0.63185933) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3591101) q[0];
sx q[0];
rz(-1.6164653) q[0];
sx q[0];
rz(0.029687667) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(0.15304676) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51296556) q[1];
sx q[1];
rz(-1.291853) q[1];
sx q[1];
rz(-0.97475027) q[1];
rz(-pi) q[2];
rz(0.072237416) q[3];
sx q[3];
rz(-2.9655955) q[3];
sx q[3];
rz(-2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(-2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8343617) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(-2.2063072) q[0];
x q[1];
rz(1.8413576) q[2];
sx q[2];
rz(-0.62108835) q[2];
sx q[2];
rz(-1.5121216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.084701531) q[1];
sx q[1];
rz(-1.4086205) q[1];
sx q[1];
rz(-0.2918891) q[1];
rz(-0.018499231) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95773762) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366074) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.8100912) q[0];
x q[1];
rz(-3.0954103) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(1.6398167) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.096958728) q[1];
sx q[1];
rz(-0.63475906) q[1];
sx q[1];
rz(1.8083014) q[1];
rz(-pi) q[2];
rz(-0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(-1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.404495) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.6361902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.226107) q[0];
rz(0.6263528) q[2];
sx q[2];
rz(-0.92712958) q[2];
sx q[2];
rz(-2.6477637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(3.0909096) q[1];
rz(-pi) q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.8441895) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18683534) q[0];
sx q[0];
rz(-1.9592013) q[0];
sx q[0];
rz(-0.79083058) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55391295) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(1.5253138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39667323) q[1];
sx q[1];
rz(-1.6066237) q[1];
sx q[1];
rz(0.14573914) q[1];
rz(1.7725138) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6181347) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20477644) q[0];
sx q[0];
rz(-1.6813155) q[0];
sx q[0];
rz(-0.21893455) q[0];
x q[1];
rz(-0.98430888) q[2];
sx q[2];
rz(-1.7585635) q[2];
sx q[2];
rz(-0.078660065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-0.34983695) q[1];
x q[2];
rz(2.0189832) q[3];
sx q[3];
rz(-0.19482329) q[3];
sx q[3];
rz(-2.4217055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-2.0586661) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(-0.54995723) q[3];
sx q[3];
rz(-0.64823845) q[3];
sx q[3];
rz(1.3241495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
