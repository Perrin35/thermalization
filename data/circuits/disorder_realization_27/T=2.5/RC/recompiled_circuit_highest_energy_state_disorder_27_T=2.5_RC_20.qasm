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
rz(-1.6838411) q[0];
sx q[0];
rz(-0.10310752) q[0];
sx q[0];
rz(0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(-2.4060251) q[1];
sx q[1];
rz(-0.29677376) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4838801) q[0];
sx q[0];
rz(-1.1917416) q[0];
sx q[0];
rz(-1.1227648) q[0];
rz(-0.71798433) q[2];
sx q[2];
rz(-1.7545926) q[2];
sx q[2];
rz(1.2280457) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.709503) q[1];
sx q[1];
rz(-2.5910834) q[1];
sx q[1];
rz(0.014660346) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9624309) q[3];
sx q[3];
rz(-1.1350313) q[3];
sx q[3];
rz(0.023438862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.890581) q[2];
sx q[2];
rz(-3.0457532) q[2];
sx q[2];
rz(-1.7354234) q[2];
rz(-2.3928394) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(-0.21193084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9894079) q[0];
sx q[0];
rz(-0.87225544) q[0];
sx q[0];
rz(2.8023791) q[0];
rz(2.5665414) q[1];
sx q[1];
rz(-1.9564068) q[1];
sx q[1];
rz(0.38160479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763827) q[0];
sx q[0];
rz(-1.8032224) q[0];
sx q[0];
rz(-1.1430955) q[0];
rz(-pi) q[1];
rz(-1.9100045) q[2];
sx q[2];
rz(-1.8198065) q[2];
sx q[2];
rz(-0.23819645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84880616) q[1];
sx q[1];
rz(-1.5611575) q[1];
sx q[1];
rz(1.9697492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6032002) q[3];
sx q[3];
rz(-0.88148553) q[3];
sx q[3];
rz(-1.4581084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7316458) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(0.48639578) q[2];
rz(-0.54388034) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(2.9774408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51089066) q[0];
sx q[0];
rz(-2.7257901) q[0];
sx q[0];
rz(-2.1422332) q[0];
rz(0.73127812) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(2.9670002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098487735) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(-1.5652324) q[0];
x q[1];
rz(-0.83589696) q[2];
sx q[2];
rz(-0.23698254) q[2];
sx q[2];
rz(-1.3224755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0831956) q[1];
sx q[1];
rz(-0.10279142) q[1];
sx q[1];
rz(-2.1071069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1483795) q[3];
sx q[3];
rz(-2.6452521) q[3];
sx q[3];
rz(-1.6272735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.752855) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(1.8641776) q[2];
rz(-2.3169005) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326062) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(-1.9933568) q[0];
rz(2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(1.0848328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6713632) q[0];
sx q[0];
rz(-1.4309023) q[0];
sx q[0];
rz(-1.2752007) q[0];
x q[1];
rz(2.8021028) q[2];
sx q[2];
rz(-0.84622806) q[2];
sx q[2];
rz(0.17865114) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12382896) q[1];
sx q[1];
rz(-2.4014086) q[1];
sx q[1];
rz(-3.0637118) q[1];
rz(-pi) q[2];
rz(0.46369073) q[3];
sx q[3];
rz(-1.9640018) q[3];
sx q[3];
rz(-0.56305199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3499902) q[2];
sx q[2];
rz(-2.377066) q[2];
sx q[2];
rz(0.0066268607) q[2];
rz(2.918112) q[3];
sx q[3];
rz(-2.1036744) q[3];
sx q[3];
rz(2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70581907) q[0];
sx q[0];
rz(-0.75278246) q[0];
sx q[0];
rz(-1.9594132) q[0];
rz(1.301282) q[1];
sx q[1];
rz(-2.2186406) q[1];
sx q[1];
rz(-0.15929793) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8790588) q[0];
sx q[0];
rz(-0.28994432) q[0];
sx q[0];
rz(3.0109809) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1012405) q[2];
sx q[2];
rz(-2.7795876) q[2];
sx q[2];
rz(-1.9792045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9827798) q[1];
sx q[1];
rz(-0.69852622) q[1];
sx q[1];
rz(2.34312) q[1];
x q[2];
rz(2.716655) q[3];
sx q[3];
rz(-2.1930088) q[3];
sx q[3];
rz(3.034301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2948239) q[2];
sx q[2];
rz(-2.939665) q[2];
sx q[2];
rz(-1.798604) q[2];
rz(-0.35074562) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(2.8158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2419234) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(0.1668461) q[0];
rz(-2.9150561) q[1];
sx q[1];
rz(-1.8140565) q[1];
sx q[1];
rz(-2.66364) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.307823) q[0];
sx q[0];
rz(-1.9109028) q[0];
sx q[0];
rz(-2.1653963) q[0];
x q[1];
rz(-1.9327546) q[2];
sx q[2];
rz(-2.2565292) q[2];
sx q[2];
rz(-2.509523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9354269) q[1];
sx q[1];
rz(-0.7784673) q[1];
sx q[1];
rz(-2.4099518) q[1];
rz(-2.7109954) q[3];
sx q[3];
rz(-1.1599419) q[3];
sx q[3];
rz(-1.915286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53092521) q[2];
sx q[2];
rz(-1.764856) q[2];
sx q[2];
rz(1.2923856) q[2];
rz(-2.2409706) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6677299) q[0];
sx q[0];
rz(-2.168499) q[0];
sx q[0];
rz(1.5245755) q[0];
rz(1.698311) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(2.6406094) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9568142) q[0];
sx q[0];
rz(-1.2010472) q[0];
sx q[0];
rz(2.2201204) q[0];
x q[1];
rz(0.47246859) q[2];
sx q[2];
rz(-0.61206619) q[2];
sx q[2];
rz(-1.2304359) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1107378) q[1];
sx q[1];
rz(-2.1424865) q[1];
sx q[1];
rz(2.5377889) q[1];
x q[2];
rz(0.30876183) q[3];
sx q[3];
rz(-1.4523066) q[3];
sx q[3];
rz(-2.091696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3524126) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(0.46335709) q[2];
rz(-3.0739259) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25245923) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(2.4976318) q[1];
sx q[1];
rz(-2.0168346) q[1];
sx q[1];
rz(0.28800979) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.889042) q[0];
sx q[0];
rz(-1.6091378) q[0];
sx q[0];
rz(-0.19186963) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3013822) q[2];
sx q[2];
rz(-2.3019493) q[2];
sx q[2];
rz(-2.3498442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0916413) q[1];
sx q[1];
rz(-2.5041951) q[1];
sx q[1];
rz(0.47224381) q[1];
x q[2];
rz(-2.8657416) q[3];
sx q[3];
rz(-2.7497661) q[3];
sx q[3];
rz(2.0656757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94706941) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(-0.21491773) q[2];
rz(1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(-2.9609093) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27338481) q[0];
sx q[0];
rz(-0.93623638) q[0];
sx q[0];
rz(2.0599763) q[0];
rz(-2.1514905) q[1];
sx q[1];
rz(-1.5568045) q[1];
sx q[1];
rz(0.33499151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7107527) q[0];
sx q[0];
rz(-2.8169605) q[0];
sx q[0];
rz(-1.6765094) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1352409) q[2];
sx q[2];
rz(-1.1489023) q[2];
sx q[2];
rz(-0.61074257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3051863) q[1];
sx q[1];
rz(-1.6560153) q[1];
sx q[1];
rz(-3.042074) q[1];
x q[2];
rz(1.1698059) q[3];
sx q[3];
rz(-1.5785517) q[3];
sx q[3];
rz(1.4435651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0830903) q[2];
sx q[2];
rz(-2.6164656) q[2];
sx q[2];
rz(-0.38880175) q[2];
rz(-2.7723516) q[3];
sx q[3];
rz(-0.25596127) q[3];
sx q[3];
rz(-0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-2.9487172) q[0];
sx q[0];
rz(-2.188864) q[0];
sx q[0];
rz(0.70190758) q[0];
rz(-1.6096055) q[1];
sx q[1];
rz(-0.67370266) q[1];
sx q[1];
rz(2.7598377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696357) q[0];
sx q[0];
rz(-1.5077295) q[0];
sx q[0];
rz(-0.081916787) q[0];
rz(-pi) q[1];
rz(1.4456621) q[2];
sx q[2];
rz(-1.9788673) q[2];
sx q[2];
rz(3.0681075) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3249141) q[1];
sx q[1];
rz(-0.35129282) q[1];
sx q[1];
rz(-0.80609806) q[1];
rz(-0.89729805) q[3];
sx q[3];
rz(-0.89400154) q[3];
sx q[3];
rz(-0.38161665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4664885) q[2];
sx q[2];
rz(-2.1110822) q[2];
sx q[2];
rz(0.64811903) q[2];
rz(-1.0068007) q[3];
sx q[3];
rz(-1.9961793) q[3];
sx q[3];
rz(-3.0692611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61160144) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(-2.8499659) q[1];
sx q[1];
rz(-1.3529774) q[1];
sx q[1];
rz(-1.1381961) q[1];
rz(-1.2164581) q[2];
sx q[2];
rz(-1.0783429) q[2];
sx q[2];
rz(1.1781296) q[2];
rz(1.9907822) q[3];
sx q[3];
rz(-1.715015) q[3];
sx q[3];
rz(0.7072995) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
