OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3381391) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(1.8344318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88959496) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(-1.3537458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-0.70322137) q[1];
sx q[1];
rz(-2.1232848) q[1];
rz(-pi) q[2];
rz(2.6775042) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(-0.97809631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933724) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63576525) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(-2.2155227) q[0];
rz(0.26308665) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(0.69743246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.038139548) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(-1.000688) q[1];
rz(-0.52069943) q[3];
sx q[3];
rz(-1.4200913) q[3];
sx q[3];
rz(-3.0050788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-0.48167357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96268481) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(0.90895598) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1842314) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(2.7797109) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.958365) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(2.4135713) q[1];
rz(-pi) q[2];
rz(0.62044575) q[3];
sx q[3];
rz(-2.1944322) q[3];
sx q[3];
rz(2.7321531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(0.552185) q[0];
rz(1.5532956) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.8968556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23403215) q[0];
sx q[0];
rz(-1.2313156) q[0];
sx q[0];
rz(-1.5690804) q[0];
x q[1];
rz(-0.99889465) q[2];
sx q[2];
rz(-0.23795393) q[2];
sx q[2];
rz(1.1513125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6914122) q[1];
sx q[1];
rz(-1.0346864) q[1];
sx q[1];
rz(0.24713534) q[1];
rz(-pi) q[2];
rz(1.3454868) q[3];
sx q[3];
rz(-2.1300737) q[3];
sx q[3];
rz(0.51636458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(-1.6385471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25795709) q[0];
sx q[0];
rz(-0.99146087) q[0];
sx q[0];
rz(-0.15276076) q[0];
x q[1];
rz(-0.083085255) q[2];
sx q[2];
rz(-1.462888) q[2];
sx q[2];
rz(-1.8748869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(-1.2574408) q[1];
x q[2];
rz(-1.8893858) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6333255) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(-0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5313107) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(2.3430603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079433867) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(0.14607231) q[0];
rz(1.919826) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(0.40086056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27998566) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(3.0793889) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1875601) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(-1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-2.7748761) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-2.5352056) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-0.95373017) q[0];
sx q[0];
rz(-2.6130136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8115933) q[2];
sx q[2];
rz(-1.5903928) q[2];
sx q[2];
rz(-0.400825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53941599) q[1];
sx q[1];
rz(-2.5578941) q[1];
sx q[1];
rz(-1.0023487) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(-1.9627278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-2.8835473) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(-2.7602957) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.4415178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9863319) q[0];
sx q[0];
rz(-3.0763456) q[0];
sx q[0];
rz(2.9078527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5324138) q[2];
sx q[2];
rz(-0.95853171) q[2];
sx q[2];
rz(-1.4505163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4691094) q[1];
sx q[1];
rz(-1.6961349) q[1];
sx q[1];
rz(-0.13087665) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5278682) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(-1.7121401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.3990078) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-2.1549966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50842677) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(1.7500061) q[0];
x q[1];
rz(-2.2676342) q[2];
sx q[2];
rz(-1.1468337) q[2];
sx q[2];
rz(-1.6800113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66227312) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(2.8981266) q[1];
rz(-pi) q[2];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.4126561) q[3];
sx q[3];
rz(1.4827673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9916519) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.8206235) q[0];
rz(-pi) q[1];
rz(0.24869463) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(-2.8949708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0701329) q[1];
sx q[1];
rz(-1.3576641) q[1];
sx q[1];
rz(-0.15503426) q[1];
rz(-pi) q[2];
rz(-0.16898445) q[3];
sx q[3];
rz(-2.0770693) q[3];
sx q[3];
rz(0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(-0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-2.5406191) q[2];
sx q[2];
rz(-0.31243159) q[2];
sx q[2];
rz(2.5675788) q[2];
rz(-2.9028805) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
