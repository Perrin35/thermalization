OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8747099) q[0];
sx q[0];
rz(-2.0800135) q[0];
sx q[0];
rz(0.51577407) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(-0.20144784) q[1];
sx q[1];
rz(-1.0010866) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37487113) q[0];
sx q[0];
rz(-2.3618638) q[0];
sx q[0];
rz(-2.4089902) q[0];
x q[1];
rz(0.41462173) q[2];
sx q[2];
rz(-1.5650563) q[2];
sx q[2];
rz(1.8513377) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6382003) q[1];
sx q[1];
rz(-0.93825785) q[1];
sx q[1];
rz(2.1293541) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9941213) q[3];
sx q[3];
rz(-1.4502954) q[3];
sx q[3];
rz(-1.7013753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4150245) q[2];
sx q[2];
rz(-0.51076204) q[2];
sx q[2];
rz(1.4600352) q[2];
rz(0.36871746) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(-1.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52861315) q[0];
sx q[0];
rz(-0.11358417) q[0];
sx q[0];
rz(3.0533277) q[0];
rz(2.2093692) q[1];
sx q[1];
rz(-1.1270707) q[1];
sx q[1];
rz(0.50311911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4012484) q[0];
sx q[0];
rz(-0.49763864) q[0];
sx q[0];
rz(-2.9236682) q[0];
x q[1];
rz(0.90089779) q[2];
sx q[2];
rz(-1.42282) q[2];
sx q[2];
rz(2.072538) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4331665) q[1];
sx q[1];
rz(-1.462877) q[1];
sx q[1];
rz(0.34924653) q[1];
x q[2];
rz(2.3130619) q[3];
sx q[3];
rz(-1.9779543) q[3];
sx q[3];
rz(0.26155868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3162389) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(-2.6884354) q[2];
rz(-3.1406) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(-0.7411952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93333018) q[0];
sx q[0];
rz(-2.9587511) q[0];
sx q[0];
rz(0.46874794) q[0];
rz(2.5543429) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(2.8900878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5875781) q[0];
sx q[0];
rz(-1.9791326) q[0];
sx q[0];
rz(0.46550444) q[0];
rz(-2.7526546) q[2];
sx q[2];
rz(-1.374482) q[2];
sx q[2];
rz(2.0816572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1826102) q[1];
sx q[1];
rz(-1.512456) q[1];
sx q[1];
rz(-1.6946409) q[1];
rz(-2.4819314) q[3];
sx q[3];
rz(-0.53852496) q[3];
sx q[3];
rz(-0.3769359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59715366) q[2];
sx q[2];
rz(-0.73828283) q[2];
sx q[2];
rz(-0.047133751) q[2];
rz(-0.86826396) q[3];
sx q[3];
rz(-1.9009813) q[3];
sx q[3];
rz(2.6422083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.353001) q[0];
sx q[0];
rz(-2.7356739) q[0];
sx q[0];
rz(1.8274008) q[0];
rz(-0.81002533) q[1];
sx q[1];
rz(-1.6255197) q[1];
sx q[1];
rz(-0.31585082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0266036) q[0];
sx q[0];
rz(-2.3427377) q[0];
sx q[0];
rz(-2.3197915) q[0];
rz(-0.48607488) q[2];
sx q[2];
rz(-0.3559843) q[2];
sx q[2];
rz(0.69185585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3923287) q[1];
sx q[1];
rz(-2.3468319) q[1];
sx q[1];
rz(-0.41458798) q[1];
x q[2];
rz(-1.4203912) q[3];
sx q[3];
rz(-0.97109761) q[3];
sx q[3];
rz(0.57460845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4874068) q[2];
sx q[2];
rz(-2.5155289) q[2];
sx q[2];
rz(0.89228863) q[2];
rz(-1.4687294) q[3];
sx q[3];
rz(-2.0817616) q[3];
sx q[3];
rz(1.0631801) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6650218) q[0];
sx q[0];
rz(-2.5044818) q[0];
sx q[0];
rz(-2.2549905) q[0];
rz(0.84238148) q[1];
sx q[1];
rz(-1.8319172) q[1];
sx q[1];
rz(1.1096035) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1177301) q[0];
sx q[0];
rz(-0.95265022) q[0];
sx q[0];
rz(-2.6075075) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9363052) q[2];
sx q[2];
rz(-2.4078806) q[2];
sx q[2];
rz(3.044341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2354616) q[1];
sx q[1];
rz(-2.0812199) q[1];
sx q[1];
rz(2.4861369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97115406) q[3];
sx q[3];
rz(-2.1622938) q[3];
sx q[3];
rz(-0.20918555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8699009) q[2];
sx q[2];
rz(-0.39125189) q[2];
sx q[2];
rz(0.58763495) q[2];
rz(0.21688004) q[3];
sx q[3];
rz(-1.2809332) q[3];
sx q[3];
rz(1.3542401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13123913) q[0];
sx q[0];
rz(-2.6709747) q[0];
sx q[0];
rz(1.3838029) q[0];
rz(-0.17419392) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(-1.1879638) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0452703) q[0];
sx q[0];
rz(-2.850527) q[0];
sx q[0];
rz(-1.2805268) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30004487) q[2];
sx q[2];
rz(-0.52562537) q[2];
sx q[2];
rz(0.92074652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9932258) q[1];
sx q[1];
rz(-1.6794551) q[1];
sx q[1];
rz(-3.1106366) q[1];
rz(0.76455595) q[3];
sx q[3];
rz(-1.3455794) q[3];
sx q[3];
rz(1.5080687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69372988) q[2];
sx q[2];
rz(-1.8665946) q[2];
sx q[2];
rz(1.0542144) q[2];
rz(-0.55274719) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(-2.0235846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0746821) q[0];
sx q[0];
rz(-0.56147611) q[0];
sx q[0];
rz(-3.1232324) q[0];
rz(-1.8439937) q[1];
sx q[1];
rz(-1.3341787) q[1];
sx q[1];
rz(-2.0265354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92662382) q[0];
sx q[0];
rz(-1.4415662) q[0];
sx q[0];
rz(1.5885167) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1586559) q[2];
sx q[2];
rz(-0.52861428) q[2];
sx q[2];
rz(1.7129218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9950986) q[1];
sx q[1];
rz(-2.0420688) q[1];
sx q[1];
rz(-2.5033094) q[1];
rz(-pi) q[2];
rz(-2.2998718) q[3];
sx q[3];
rz(-2.5091565) q[3];
sx q[3];
rz(3.1068592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6899507) q[2];
sx q[2];
rz(-0.74863282) q[2];
sx q[2];
rz(-2.3484223) q[2];
rz(-2.4801109) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(-2.49559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20872214) q[0];
sx q[0];
rz(-2.0168309) q[0];
sx q[0];
rz(-2.9710508) q[0];
rz(1.0199245) q[1];
sx q[1];
rz(-0.41671697) q[1];
sx q[1];
rz(0.65933093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.764544) q[0];
sx q[0];
rz(-1.555849) q[0];
sx q[0];
rz(0.055146761) q[0];
rz(-pi) q[1];
rz(2.6213367) q[2];
sx q[2];
rz(-1.1197436) q[2];
sx q[2];
rz(0.28869155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48261343) q[1];
sx q[1];
rz(-0.70090997) q[1];
sx q[1];
rz(2.6644876) q[1];
rz(2.5468656) q[3];
sx q[3];
rz(-0.85053125) q[3];
sx q[3];
rz(1.8320465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3063804) q[2];
sx q[2];
rz(-0.39432085) q[2];
sx q[2];
rz(2.1739056) q[2];
rz(-0.59099284) q[3];
sx q[3];
rz(-1.7715745) q[3];
sx q[3];
rz(1.7041357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10385926) q[0];
sx q[0];
rz(-0.93872207) q[0];
sx q[0];
rz(-2.9528604) q[0];
rz(-2.0813148) q[1];
sx q[1];
rz(-0.66245586) q[1];
sx q[1];
rz(0.8998543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8322139) q[0];
sx q[0];
rz(-0.77118528) q[0];
sx q[0];
rz(0.82892771) q[0];
rz(-pi) q[1];
rz(1.5737888) q[2];
sx q[2];
rz(-2.4140777) q[2];
sx q[2];
rz(-1.2902009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65727216) q[1];
sx q[1];
rz(-1.5725699) q[1];
sx q[1];
rz(-2.491105) q[1];
rz(2.9505355) q[3];
sx q[3];
rz(-0.85874288) q[3];
sx q[3];
rz(-2.4993757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4047644) q[2];
sx q[2];
rz(-1.9177723) q[2];
sx q[2];
rz(-1.169211) q[2];
rz(0.87559187) q[3];
sx q[3];
rz(-0.67684042) q[3];
sx q[3];
rz(0.08918795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.5065696) q[0];
sx q[0];
rz(-1.4590141) q[0];
sx q[0];
rz(2.7154679) q[0];
rz(-1.9020853) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(0.32136163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023055619) q[0];
sx q[0];
rz(-1.5154953) q[0];
sx q[0];
rz(0.99744101) q[0];
rz(-pi) q[1];
rz(-0.05409492) q[2];
sx q[2];
rz(-2.0129497) q[2];
sx q[2];
rz(-1.4070321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9656941) q[1];
sx q[1];
rz(-1.7135317) q[1];
sx q[1];
rz(1.229415) q[1];
x q[2];
rz(-2.5697702) q[3];
sx q[3];
rz(-0.79872978) q[3];
sx q[3];
rz(-0.031599061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1539803) q[2];
sx q[2];
rz(-0.93903792) q[2];
sx q[2];
rz(-0.067848094) q[2];
rz(2.1765354) q[3];
sx q[3];
rz(-2.8128251) q[3];
sx q[3];
rz(-1.4062101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6168552) q[0];
sx q[0];
rz(-2.6829834) q[0];
sx q[0];
rz(0.99880698) q[0];
rz(-0.85159341) q[1];
sx q[1];
rz(-2.0709745) q[1];
sx q[1];
rz(-2.4877683) q[1];
rz(2.5851696) q[2];
sx q[2];
rz(-1.5172554) q[2];
sx q[2];
rz(0.42419626) q[2];
rz(1.5953993) q[3];
sx q[3];
rz(-2.2886124) q[3];
sx q[3];
rz(-0.99289865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
