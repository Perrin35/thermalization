OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(-0.85715357) q[0];
sx q[0];
rz(0.13248086) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1231175) q[0];
sx q[0];
rz(-1.0536061) q[0];
sx q[0];
rz(3.0267107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1039256) q[2];
sx q[2];
rz(-2.1618927) q[2];
sx q[2];
rz(0.95962722) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9308656) q[1];
sx q[1];
rz(-0.14666808) q[1];
sx q[1];
rz(-1.055483) q[1];
rz(-pi) q[2];
rz(2.0308706) q[3];
sx q[3];
rz(-0.51957909) q[3];
sx q[3];
rz(0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(1.8923627) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55727977) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(0.31038196) q[0];
rz(-pi) q[1];
rz(-1.3524019) q[2];
sx q[2];
rz(-2.1973655) q[2];
sx q[2];
rz(0.44109694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7932574) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(2.4130915) q[1];
x q[2];
rz(-0.35238102) q[3];
sx q[3];
rz(-1.1412732) q[3];
sx q[3];
rz(-2.9126715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(0.87835971) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(-0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(0.032827854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093124495) q[0];
sx q[0];
rz(-1.801352) q[0];
sx q[0];
rz(0.31974098) q[0];
rz(-2.8509154) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(-0.79613396) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(-2.3418531) q[1];
x q[2];
rz(-1.0333943) q[3];
sx q[3];
rz(-0.61363797) q[3];
sx q[3];
rz(2.9463241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34439987) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(-2.8642505) q[2];
rz(0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(0.77082005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(-0.016419134) q[0];
x q[1];
rz(-1.7647676) q[2];
sx q[2];
rz(-0.49805222) q[2];
sx q[2];
rz(1.1454524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3168287) q[1];
sx q[1];
rz(-1.0526122) q[1];
sx q[1];
rz(-1.9903241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7951489) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(-2.2546774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(0.17000155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572896) q[0];
sx q[0];
rz(-2.347725) q[0];
sx q[0];
rz(-2.8110709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0427809) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(-0.99265487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0712142) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(2.1945164) q[1];
x q[2];
rz(2.4835028) q[3];
sx q[3];
rz(-2.1562025) q[3];
sx q[3];
rz(0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(0.68516723) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2049094) q[0];
sx q[0];
rz(-0.65403599) q[0];
sx q[0];
rz(1.685164) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77722705) q[2];
sx q[2];
rz(-1.4677375) q[2];
sx q[2];
rz(1.7412141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8243858) q[1];
sx q[1];
rz(-0.82815352) q[1];
sx q[1];
rz(-0.87888996) q[1];
x q[2];
rz(-1.5983726) q[3];
sx q[3];
rz(-0.58610361) q[3];
sx q[3];
rz(-0.54857777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(-1.7123327) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(-0.72934735) q[0];
rz(-0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.9940631) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.731819) q[0];
sx q[0];
rz(-1.496135) q[0];
sx q[0];
rz(0.41418196) q[0];
rz(-2.2044264) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(-2.5123951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8116248) q[1];
sx q[1];
rz(-2.6553272) q[1];
sx q[1];
rz(-3.0537513) q[1];
rz(1.958193) q[3];
sx q[3];
rz(-2.0347188) q[3];
sx q[3];
rz(1.0138318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(-2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(0.83129445) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-2.7430699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(-2.0429862) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4648828) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(-2.523409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9989222) q[1];
sx q[1];
rz(-0.97375662) q[1];
sx q[1];
rz(-0.02762694) q[1];
rz(2.1274444) q[3];
sx q[3];
rz(-2.017052) q[3];
sx q[3];
rz(-2.2198912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(2.0166345) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(-1.7522316) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-1.0983889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84868816) q[0];
sx q[0];
rz(-1.0875889) q[0];
sx q[0];
rz(-0.26341652) q[0];
rz(-pi) q[1];
rz(-2.6838052) q[2];
sx q[2];
rz(-2.2676003) q[2];
sx q[2];
rz(-0.84504499) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92139447) q[1];
sx q[1];
rz(-2.1500593) q[1];
sx q[1];
rz(1.5701576) q[1];
rz(2.8875651) q[3];
sx q[3];
rz(-1.4048178) q[3];
sx q[3];
rz(-2.0167054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(2.675132) q[0];
rz(2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(0.62896532) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53101978) q[0];
sx q[0];
rz(-1.5218698) q[0];
sx q[0];
rz(1.1958836) q[0];
x q[1];
rz(-2.2132657) q[2];
sx q[2];
rz(-1.7648342) q[2];
sx q[2];
rz(-2.5460668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91937376) q[1];
sx q[1];
rz(-2.0128951) q[1];
sx q[1];
rz(1.5685146) q[1];
x q[2];
rz(1.2260776) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(-2.7310839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8978867) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(-1.0591327) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(-1.5268486) q[2];
sx q[2];
rz(-2.2938001) q[2];
sx q[2];
rz(-3.0707035) q[2];
rz(-2.1577088) q[3];
sx q[3];
rz(-1.513653) q[3];
sx q[3];
rz(-0.37099864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];