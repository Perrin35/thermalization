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
rz(-0.28873697) q[0];
sx q[0];
rz(-2.4462235) q[0];
sx q[0];
rz(-0.26779548) q[0];
rz(-2.7195622) q[1];
sx q[1];
rz(-0.92075092) q[1];
sx q[1];
rz(-1.2738127) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201694) q[0];
sx q[0];
rz(-1.5258938) q[0];
sx q[0];
rz(-1.8123167) q[0];
x q[1];
rz(0.8795514) q[2];
sx q[2];
rz(-2.3143907) q[2];
sx q[2];
rz(0.66397053) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2886658) q[1];
sx q[1];
rz(-1.5024606) q[1];
sx q[1];
rz(2.4803376) q[1];
rz(-pi) q[2];
rz(1.7142322) q[3];
sx q[3];
rz(-1.6063476) q[3];
sx q[3];
rz(1.2811023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6196809) q[2];
sx q[2];
rz(-0.64138594) q[2];
sx q[2];
rz(-0.042595159) q[2];
rz(2.8604782) q[3];
sx q[3];
rz(-1.576985) q[3];
sx q[3];
rz(-2.5458096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6302781) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(1.0761155) q[0];
rz(1.0307182) q[1];
sx q[1];
rz(-2.0298256) q[1];
sx q[1];
rz(-1.8310742) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1210107) q[0];
sx q[0];
rz(-2.2288929) q[0];
sx q[0];
rz(-0.026348635) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83514799) q[2];
sx q[2];
rz(-1.2642167) q[2];
sx q[2];
rz(-0.94294244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0243822) q[1];
sx q[1];
rz(-1.6892994) q[1];
sx q[1];
rz(-2.973053) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53383975) q[3];
sx q[3];
rz(-1.2152142) q[3];
sx q[3];
rz(-1.8227641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92404667) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(0.95620608) q[2];
rz(1.8337967) q[3];
sx q[3];
rz(-0.81495133) q[3];
sx q[3];
rz(0.1213049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8420551) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(3.1053542) q[0];
rz(0.80129519) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(1.3005728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0220003) q[0];
sx q[0];
rz(-1.4109857) q[0];
sx q[0];
rz(3.0335866) q[0];
rz(-pi) q[1];
rz(-2.168) q[2];
sx q[2];
rz(-1.6064715) q[2];
sx q[2];
rz(3.1247471) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6319148) q[1];
sx q[1];
rz(-2.1304806) q[1];
sx q[1];
rz(2.6259929) q[1];
rz(-0.68978975) q[3];
sx q[3];
rz(-1.8588603) q[3];
sx q[3];
rz(0.40460872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3761882) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-2.9236531) q[2];
rz(0.91207063) q[3];
sx q[3];
rz(-1.6488766) q[3];
sx q[3];
rz(-0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4172149) q[0];
sx q[0];
rz(-1.0715002) q[0];
sx q[0];
rz(0.81556129) q[0];
rz(-2.0888445) q[1];
sx q[1];
rz(-0.88200724) q[1];
sx q[1];
rz(-1.7900593) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4824351) q[0];
sx q[0];
rz(-1.1023942) q[0];
sx q[0];
rz(1.4841561) q[0];
rz(-pi) q[1];
rz(1.4830681) q[2];
sx q[2];
rz(-2.2854439) q[2];
sx q[2];
rz(1.8217979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.4554418) q[1];
sx q[1];
rz(-2.0459922) q[1];
sx q[1];
rz(-3.1093756) q[1];
rz(-pi) q[2];
rz(-0.71664401) q[3];
sx q[3];
rz(-1.1782681) q[3];
sx q[3];
rz(-1.5115304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.012933) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(0.43221691) q[2];
rz(-0.84248078) q[3];
sx q[3];
rz(-1.0336927) q[3];
sx q[3];
rz(-2.1459818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1603482) q[0];
sx q[0];
rz(-2.6762185) q[0];
sx q[0];
rz(-0.55111849) q[0];
rz(0.61800686) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(2.0726223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3270996) q[0];
sx q[0];
rz(-2.2393423) q[0];
sx q[0];
rz(-1.8446246) q[0];
rz(-pi) q[1];
rz(2.3927116) q[2];
sx q[2];
rz(-0.62070337) q[2];
sx q[2];
rz(1.701603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6275989) q[1];
sx q[1];
rz(-0.74666903) q[1];
sx q[1];
rz(-1.385266) q[1];
x q[2];
rz(0.25311796) q[3];
sx q[3];
rz(-2.6147644) q[3];
sx q[3];
rz(1.4344858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7908287) q[2];
sx q[2];
rz(-2.5574234) q[2];
sx q[2];
rz(-2.186415) q[2];
rz(0.59349924) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(-1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3747568) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(1.391885) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(1.8291738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8439633) q[0];
sx q[0];
rz(-2.2052551) q[0];
sx q[0];
rz(-2.5353801) q[0];
rz(-pi) q[1];
rz(2.2696804) q[2];
sx q[2];
rz(-1.0191227) q[2];
sx q[2];
rz(-1.0629423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3692664) q[1];
sx q[1];
rz(-1.516894) q[1];
sx q[1];
rz(0.51477706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3950609) q[3];
sx q[3];
rz(-2.8214957) q[3];
sx q[3];
rz(-1.2552346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.352508) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(-1.7279203) q[2];
rz(0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(-0.55268923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308052) q[0];
sx q[0];
rz(-3.0785705) q[0];
sx q[0];
rz(-2.4600273) q[0];
rz(-1.1850146) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(0.78651816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98726051) q[0];
sx q[0];
rz(-1.2367147) q[0];
sx q[0];
rz(2.7272237) q[0];
rz(-2.2165197) q[2];
sx q[2];
rz(-0.73410119) q[2];
sx q[2];
rz(0.9303329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1934172) q[1];
sx q[1];
rz(-1.4637636) q[1];
sx q[1];
rz(2.2501037) q[1];
rz(-1.6166572) q[3];
sx q[3];
rz(-0.73165441) q[3];
sx q[3];
rz(-1.589244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6098392) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(1.5717724) q[2];
rz(1.8543367) q[3];
sx q[3];
rz(-1.5232892) q[3];
sx q[3];
rz(0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.542273) q[0];
sx q[0];
rz(-0.71745187) q[0];
sx q[0];
rz(1.188311) q[0];
rz(-3.1207454) q[1];
sx q[1];
rz(-0.7531082) q[1];
sx q[1];
rz(2.5426224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99020236) q[0];
sx q[0];
rz(-1.4722451) q[0];
sx q[0];
rz(2.0771785) q[0];
x q[1];
rz(-0.85592593) q[2];
sx q[2];
rz(-1.2675261) q[2];
sx q[2];
rz(1.0452458) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61781949) q[1];
sx q[1];
rz(-1.0955155) q[1];
sx q[1];
rz(2.3889524) q[1];
x q[2];
rz(1.0858367) q[3];
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
rz(2.6263728) q[2];
rz(-2.6089) q[3];
sx q[3];
rz(-2.0566514) q[3];
sx q[3];
rz(-0.93713078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135947) q[0];
sx q[0];
rz(-0.6830712) q[0];
sx q[0];
rz(-2.7711476) q[0];
rz(-2.0619552) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(3.0830141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616042) q[0];
sx q[0];
rz(-1.4128886) q[0];
sx q[0];
rz(2.5270793) q[0];
rz(-pi) q[1];
rz(-1.3809105) q[2];
sx q[2];
rz(-2.9101964) q[2];
sx q[2];
rz(-0.93107241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96370041) q[1];
sx q[1];
rz(-2.3517015) q[1];
sx q[1];
rz(-0.13479418) q[1];
rz(0.58038099) q[3];
sx q[3];
rz(-0.92338744) q[3];
sx q[3];
rz(-3.1064432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8687245) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(0.33099428) q[2];
rz(-3.1281779) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(-1.0564055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57824221) q[0];
sx q[0];
rz(-2.712482) q[0];
sx q[0];
rz(-2.6934534) q[0];
rz(0.093712417) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.2154481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171612) q[0];
sx q[0];
rz(-1.6002065) q[0];
sx q[0];
rz(-2.0038811) q[0];
rz(0.39544659) q[2];
sx q[2];
rz(-0.85265358) q[2];
sx q[2];
rz(-1.2527695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3642973) q[1];
sx q[1];
rz(-2.5357995) q[1];
sx q[1];
rz(-2.0531027) q[1];
x q[2];
rz(1.8849202) q[3];
sx q[3];
rz(-2.7249955) q[3];
sx q[3];
rz(-2.9521717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1756246) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(0.69941163) q[2];
rz(-1.8661963) q[3];
sx q[3];
rz(-1.1558775) q[3];
sx q[3];
rz(1.8705961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324343) q[0];
sx q[0];
rz(-2.2812738) q[0];
sx q[0];
rz(-1.9914837) q[0];
rz(1.746183) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(1.6509612) q[2];
sx q[2];
rz(-0.90654984) q[2];
sx q[2];
rz(-0.64150099) q[2];
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
