OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8528557) q[0];
sx q[0];
rz(5.5878162) q[0];
sx q[0];
rz(6.5509808) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(4.0623436) q[1];
sx q[1];
rz(10.698591) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201694) q[0];
sx q[0];
rz(-1.5258938) q[0];
sx q[0];
rz(1.3292759) q[0];
rz(-pi) q[1];
rz(-0.8795514) q[2];
sx q[2];
rz(-2.3143907) q[2];
sx q[2];
rz(2.4776221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3706412) q[1];
sx q[1];
rz(-2.2302366) q[1];
sx q[1];
rz(1.6573011) q[1];
rz(-pi) q[2];
rz(3.1056728) q[3];
sx q[3];
rz(-1.714141) q[3];
sx q[3];
rz(-0.28456056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6196809) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(3.0989975) q[2];
rz(-0.28111449) q[3];
sx q[3];
rz(-1.576985) q[3];
sx q[3];
rz(0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6302781) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(-2.0654772) q[0];
rz(1.0307182) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(-1.3105185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56633184) q[0];
sx q[0];
rz(-1.5499513) q[0];
sx q[0];
rz(2.2290609) q[0];
rz(2.3064447) q[2];
sx q[2];
rz(-1.877376) q[2];
sx q[2];
rz(0.94294244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2952463) q[1];
sx q[1];
rz(-0.2057067) q[1];
sx q[1];
rz(0.61726112) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6077529) q[3];
sx q[3];
rz(-1.9263785) q[3];
sx q[3];
rz(-1.3188286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.217546) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(-2.1853866) q[2];
rz(-1.8337967) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(0.1213049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8420551) q[0];
sx q[0];
rz(-2.6326023) q[0];
sx q[0];
rz(3.1053542) q[0];
rz(0.80129519) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42388445) q[0];
sx q[0];
rz(-0.1926271) q[0];
sx q[0];
rz(2.1602551) q[0];
rz(-pi) q[1];
rz(1.6341798) q[2];
sx q[2];
rz(-2.5434539) q[2];
sx q[2];
rz(1.5352406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50967783) q[1];
sx q[1];
rz(-1.011112) q[1];
sx q[1];
rz(2.6259929) q[1];
x q[2];
rz(0.68978975) q[3];
sx q[3];
rz(-1.8588603) q[3];
sx q[3];
rz(-0.40460872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7654045) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(2.9236531) q[2];
rz(2.229522) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(-0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4172149) q[0];
sx q[0];
rz(-1.0715002) q[0];
sx q[0];
rz(2.3260314) q[0];
rz(1.0527481) q[1];
sx q[1];
rz(-0.88200724) q[1];
sx q[1];
rz(-1.7900593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4824351) q[0];
sx q[0];
rz(-2.0391984) q[0];
sx q[0];
rz(-1.6574366) q[0];
x q[1];
rz(-2.4250373) q[2];
sx q[2];
rz(-1.5045696) q[2];
sx q[2];
rz(2.9481681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4554418) q[1];
sx q[1];
rz(-1.0956005) q[1];
sx q[1];
rz(3.1093756) q[1];
rz(-pi) q[2];
rz(2.4249486) q[3];
sx q[3];
rz(-1.9633246) q[3];
sx q[3];
rz(1.5115304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.012933) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(-2.7093757) q[2];
rz(0.84248078) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(0.55111849) q[0];
rz(2.5235858) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(-2.0726223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38902125) q[0];
sx q[0];
rz(-2.4271936) q[0];
sx q[0];
rz(0.32984372) q[0];
rz(-2.3927116) q[2];
sx q[2];
rz(-0.62070337) q[2];
sx q[2];
rz(1.4399897) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0615017) q[1];
sx q[1];
rz(-1.6964165) q[1];
sx q[1];
rz(0.83275034) q[1];
rz(2.8884747) q[3];
sx q[3];
rz(-0.52682823) q[3];
sx q[3];
rz(1.4344858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3507639) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(-0.9551777) q[2];
rz(0.59349924) q[3];
sx q[3];
rz(-0.94624001) q[3];
sx q[3];
rz(1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3747568) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(1.7497077) q[0];
rz(1.1426686) q[1];
sx q[1];
rz(-1.7997768) q[1];
sx q[1];
rz(1.8291738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29762938) q[0];
sx q[0];
rz(-2.2052551) q[0];
sx q[0];
rz(2.5353801) q[0];
rz(-pi) q[1];
rz(-2.2696804) q[2];
sx q[2];
rz(-2.12247) q[2];
sx q[2];
rz(2.0786503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9095972) q[1];
sx q[1];
rz(-1.0568406) q[1];
sx q[1];
rz(-1.6327052) q[1];
x q[2];
rz(1.8862861) q[3];
sx q[3];
rz(-1.5157561) q[3];
sx q[3];
rz(2.6590527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.352508) q[2];
sx q[2];
rz(-2.3859873) q[2];
sx q[2];
rz(1.4136723) q[2];
rz(-0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(0.55268923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308052) q[0];
sx q[0];
rz(-3.0785705) q[0];
sx q[0];
rz(-0.68156534) q[0];
rz(-1.956578) q[1];
sx q[1];
rz(-2.3763035) q[1];
sx q[1];
rz(-2.3550745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44030066) q[0];
sx q[0];
rz(-1.1806187) q[0];
sx q[0];
rz(1.2083645) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1952536) q[2];
sx q[2];
rz(-1.1558487) q[2];
sx q[2];
rz(-3.0111661) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6326846) q[1];
sx q[1];
rz(-0.68636319) q[1];
sx q[1];
rz(-1.4014161) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83966484) q[3];
sx q[3];
rz(-1.6014301) q[3];
sx q[3];
rz(0.015686663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6098392) q[2];
sx q[2];
rz(-0.24007758) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(-1.2872559) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(-0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59931961) q[0];
sx q[0];
rz(-0.71745187) q[0];
sx q[0];
rz(1.188311) q[0];
rz(-0.020847281) q[1];
sx q[1];
rz(-0.7531082) q[1];
sx q[1];
rz(-2.5426224) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1513903) q[0];
sx q[0];
rz(-1.4722451) q[0];
sx q[0];
rz(-2.0771785) q[0];
rz(0.85592593) q[2];
sx q[2];
rz(-1.2675261) q[2];
sx q[2];
rz(-1.0452458) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5934893) q[1];
sx q[1];
rz(-0.91750359) q[1];
sx q[1];
rz(-2.1849225) q[1];
rz(-2.0557559) q[3];
sx q[3];
rz(-2.3502878) q[3];
sx q[3];
rz(-0.51998653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.282436) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(-0.51521987) q[2];
rz(0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(0.93713078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027997967) q[0];
sx q[0];
rz(-0.6830712) q[0];
sx q[0];
rz(-0.37044507) q[0];
rz(-1.0796374) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(-3.0830141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48029362) q[0];
sx q[0];
rz(-0.96503557) q[0];
sx q[0];
rz(-1.3783216) q[0];
x q[1];
rz(3.09715) q[2];
sx q[2];
rz(-1.7979597) q[2];
sx q[2];
rz(2.4054766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9875659) q[1];
sx q[1];
rz(-2.3515764) q[1];
sx q[1];
rz(1.7055737) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83567668) q[3];
sx q[3];
rz(-1.1181076) q[3];
sx q[3];
rz(-1.2293466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8687245) q[2];
sx q[2];
rz(-2.3385907) q[2];
sx q[2];
rz(-0.33099428) q[2];
rz(-0.013414772) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(-2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5633504) q[0];
sx q[0];
rz(-2.712482) q[0];
sx q[0];
rz(2.6934534) q[0];
rz(3.0478802) q[1];
sx q[1];
rz(-2.8401076) q[1];
sx q[1];
rz(-1.9261446) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5543723) q[0];
sx q[0];
rz(-1.1379114) q[0];
sx q[0];
rz(0.032399633) q[0];
rz(-pi) q[1];
rz(0.81268572) q[2];
sx q[2];
rz(-1.865109) q[2];
sx q[2];
rz(-0.049969604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3443904) q[1];
sx q[1];
rz(-1.0420402) q[1];
sx q[1];
rz(-2.8307298) q[1];
rz(-pi) q[2];
rz(3.0057109) q[3];
sx q[3];
rz(-1.965842) q[3];
sx q[3];
rz(-2.6108133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.965968) q[2];
sx q[2];
rz(-1.5886687) q[2];
sx q[2];
rz(-0.69941163) q[2];
rz(-1.8661963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(-1.8705961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70915837) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(-1.3954096) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(-3.0396661) q[2];
sx q[2];
rz(-2.473255) q[2];
sx q[2];
rz(-0.77108939) q[2];
rz(2.9085085) q[3];
sx q[3];
rz(-2.6641416) q[3];
sx q[3];
rz(-2.5269846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
