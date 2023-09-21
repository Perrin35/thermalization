OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(6.05655) q[1];
sx q[1];
rz(1.5645138) q[1];
sx q[1];
rz(5.984879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053781833) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(1.5009297) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0478893) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(1.1790438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5391985) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(-0.56166517) q[1];
x q[2];
rz(2.2536623) q[3];
sx q[3];
rz(-0.13266064) q[3];
sx q[3];
rz(-0.56548972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-0.99386627) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(0.018218536) q[0];
rz(0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(-0.47168628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52662151) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(1.660166) q[0];
rz(-pi) q[1];
rz(-2.6623146) q[2];
sx q[2];
rz(-2.5844378) q[2];
sx q[2];
rz(0.011475871) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7029611) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(2.5665934) q[1];
rz(-0.96305965) q[3];
sx q[3];
rz(-1.7000323) q[3];
sx q[3];
rz(0.78101633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91711125) q[0];
sx q[0];
rz(-2.3492976) q[0];
sx q[0];
rz(-0.86865058) q[0];
rz(-pi) q[1];
rz(0.45708926) q[2];
sx q[2];
rz(-2.3466913) q[2];
sx q[2];
rz(2.7283816) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1959343) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(-0.083421589) q[1];
rz(-pi) q[2];
x q[2];
rz(-6/(13*pi)) q[3];
sx q[3];
rz(-2.5111755) q[3];
sx q[3];
rz(-2.279225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-0.96763119) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-0.63823429) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(-1.9086054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15105948) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(0.14453669) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4013176) q[2];
sx q[2];
rz(-2.4057655) q[2];
sx q[2];
rz(2.1528113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9150881) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(1.7712797) q[1];
x q[2];
rz(2.9144985) q[3];
sx q[3];
rz(-2.2507651) q[3];
sx q[3];
rz(2.1122776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(-2.9290501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034886995) q[0];
sx q[0];
rz(-1.5910774) q[0];
sx q[0];
rz(-3.1253392) q[0];
rz(-pi) q[1];
rz(-0.96617713) q[2];
sx q[2];
rz(-1.539131) q[2];
sx q[2];
rz(-2.4633173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37413874) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(-2.2351082) q[1];
rz(-0.46194525) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(-0.40601054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-2.6203716) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(0.22512063) q[0];
rz(-1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-0.37757847) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29992732) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(-1.6542692) q[0];
rz(0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-2.6968323) q[1];
sx q[1];
rz(-1.8163535) q[1];
rz(-pi) q[2];
x q[2];
rz(0.029929786) q[3];
sx q[3];
rz(-2.3464591) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(-0.48103508) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(-0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2339904) q[0];
sx q[0];
rz(-1.4617209) q[0];
sx q[0];
rz(-1.2489737) q[0];
rz(-0.9585602) q[2];
sx q[2];
rz(-1.6974291) q[2];
sx q[2];
rz(2.3938092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2563045) q[1];
sx q[1];
rz(-1.8418152) q[1];
sx q[1];
rz(-2.6726252) q[1];
x q[2];
rz(2.3354704) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(-0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1640132) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(-2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5751858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9109089) q[0];
sx q[0];
rz(-0.77373234) q[0];
sx q[0];
rz(0.45549972) q[0];
rz(-pi) q[1];
rz(0.17922108) q[2];
sx q[2];
rz(-2.1041098) q[2];
sx q[2];
rz(-0.78782493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52289256) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(-0.79843847) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67993645) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87551293) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(1.6850527) q[0];
rz(-2.5121571) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(-1.1368407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658543) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(-1.0435186) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37461899) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(1.3499201) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5053619) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(1.865571) q[1];
x q[2];
rz(-2.4056899) q[3];
sx q[3];
rz(-1.8527485) q[3];
sx q[3];
rz(-1.499093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(1.0409522) q[2];
rz(-0.0020290931) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-0.61202234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529022) q[0];
sx q[0];
rz(-0.86581794) q[0];
sx q[0];
rz(-0.61625723) q[0];
rz(-pi) q[1];
rz(-1.3304747) q[2];
sx q[2];
rz(-1.2346621) q[2];
sx q[2];
rz(2.0517595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8814197) q[1];
sx q[1];
rz(-2.2248785) q[1];
sx q[1];
rz(1.4152848) q[1];
rz(-pi) q[2];
rz(1.9480115) q[3];
sx q[3];
rz(-2.0743138) q[3];
sx q[3];
rz(-0.52770381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(0.65336147) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(0.75469771) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(2.4738612) q[2];
sx q[2];
rz(-2.3420391) q[2];
sx q[2];
rz(2.5426368) q[2];
rz(1.2470506) q[3];
sx q[3];
rz(-1.4907881) q[3];
sx q[3];
rz(-1.1992906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
