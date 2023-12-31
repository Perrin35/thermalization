OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(0.0013874887) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35933094) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(2.9666535) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6944461) q[2];
sx q[2];
rz(-1.6086173) q[2];
sx q[2];
rz(2.56074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4692987) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(1.6507571) q[1];
rz(-pi) q[2];
rz(-0.9790768) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-2.0548267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(0.81545365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52112752) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(2.492766) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4120861) q[2];
sx q[2];
rz(-1.7765877) q[2];
sx q[2];
rz(-1.6934998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7533469) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(-1.5335598) q[1];
rz(-pi) q[2];
rz(-1.3482413) q[3];
sx q[3];
rz(-0.9369623) q[3];
sx q[3];
rz(-2.8702877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(-1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(0.99951807) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16386579) q[0];
sx q[0];
rz(-0.3361055) q[0];
sx q[0];
rz(-2.0396114) q[0];
rz(-pi) q[1];
rz(2.2314084) q[2];
sx q[2];
rz(-2.0816457) q[2];
sx q[2];
rz(2.3597033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1516583) q[1];
sx q[1];
rz(-2.5807568) q[1];
sx q[1];
rz(-2.6480617) q[1];
rz(1.0619034) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(0.27483637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66514689) q[0];
sx q[0];
rz(-0.11867141) q[0];
sx q[0];
rz(-0.84400405) q[0];
rz(-pi) q[1];
rz(1.1103815) q[2];
sx q[2];
rz(-1.5722256) q[2];
sx q[2];
rz(1.187385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3596613) q[1];
sx q[1];
rz(-2.0205824) q[1];
sx q[1];
rz(0.91314544) q[1];
x q[2];
rz(-0.12675385) q[3];
sx q[3];
rz(-0.8781913) q[3];
sx q[3];
rz(0.17514378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.5083195) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(2.9798853) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(-0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6277916) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(2.6898726) q[0];
rz(-pi) q[1];
rz(2.6078569) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(-1.5460207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21302528) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-0.43069559) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77002854) q[3];
sx q[3];
rz(-2.3838245) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(-0.3240164) q[2];
rz(1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(0.34067672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621988) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(-1.6367153) q[0];
rz(-pi) q[1];
rz(0.75603007) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(1.1083958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6319879) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.8297086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86990279) q[3];
sx q[3];
rz(-0.60855908) q[3];
sx q[3];
rz(-0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(-2.3699956) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(0.46494928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3195254) q[0];
sx q[0];
rz(-1.7626581) q[0];
sx q[0];
rz(2.0454387) q[0];
x q[1];
rz(2.2642235) q[2];
sx q[2];
rz(-2.4420028) q[2];
sx q[2];
rz(0.96044651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.834224) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(-2.6878396) q[1];
rz(-1.0602337) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(-0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6358444) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(-1.4710674) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3381091) q[2];
sx q[2];
rz(-2.8065971) q[2];
sx q[2];
rz(1.4575046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8790508) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(2.2294728) q[1];
rz(-pi) q[2];
rz(-1.7844291) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(-0.3233288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7408961) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-2.267568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066814518) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(-0.76807036) q[0];
rz(-pi) q[1];
rz(-1.1633515) q[2];
sx q[2];
rz(-1.904084) q[2];
sx q[2];
rz(2.8826706) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8928788) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(1.0700657) q[1];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(-1.1788517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261242) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(1.8490851) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5095021) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(2.0096411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19503838) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(-0.22919319) q[1];
rz(1.3445271) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.4964676) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(-0.27579565) q[3];
sx q[3];
rz(-1.7134568) q[3];
sx q[3];
rz(-1.4207763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
