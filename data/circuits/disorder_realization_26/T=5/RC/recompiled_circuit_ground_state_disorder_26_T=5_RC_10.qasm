OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3684664) q[0];
sx q[0];
rz(-2.3946895) q[0];
sx q[0];
rz(-0.84063831) q[0];
rz(-3.0199938) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(2.8425541) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5344319) q[0];
sx q[0];
rz(-2.2272155) q[0];
sx q[0];
rz(-2.5975511) q[0];
rz(1.7173355) q[2];
sx q[2];
rz(-0.49078178) q[2];
sx q[2];
rz(-2.7637568) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92235293) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(-0.24210614) q[1];
x q[2];
rz(-1.2848008) q[3];
sx q[3];
rz(-2.6702849) q[3];
sx q[3];
rz(-2.907674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3806939) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(0.32763457) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(0.49427858) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7647917) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-0.85533992) q[0];
rz(-0.0080464706) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(-2.6904552) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2349885) q[0];
sx q[0];
rz(-1.4140109) q[0];
sx q[0];
rz(-2.2995728) q[0];
rz(-2.3694384) q[2];
sx q[2];
rz(-1.4322965) q[2];
sx q[2];
rz(-2.0318299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0476183) q[1];
sx q[1];
rz(-2.6674358) q[1];
sx q[1];
rz(-2.7535901) q[1];
rz(0.47685949) q[3];
sx q[3];
rz(-1.9369159) q[3];
sx q[3];
rz(-1.9204634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96192876) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(0.12953225) q[2];
rz(2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.391908) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(2.5823197) q[0];
rz(0.094712146) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(-0.47725484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6420631) q[0];
sx q[0];
rz(-1.1746049) q[0];
sx q[0];
rz(2.8642333) q[0];
rz(-pi) q[1];
rz(2.7626286) q[2];
sx q[2];
rz(-2.7550089) q[2];
sx q[2];
rz(-1.6340337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86384799) q[1];
sx q[1];
rz(-1.7949634) q[1];
sx q[1];
rz(0.78296354) q[1];
rz(-pi) q[2];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.1616544) q[3];
sx q[3];
rz(-2.9144998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(-0.9300119) q[2];
rz(-1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(-2.1029396) q[0];
rz(2.5301798) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(-2.2183653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8002565) q[0];
sx q[0];
rz(-2.1390599) q[0];
sx q[0];
rz(0.41109127) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5973041) q[2];
sx q[2];
rz(-2.7237281) q[2];
sx q[2];
rz(-2.1919427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7103017) q[1];
sx q[1];
rz(-1.6552306) q[1];
sx q[1];
rz(-0.84655098) q[1];
x q[2];
rz(2.8037854) q[3];
sx q[3];
rz(-0.67854133) q[3];
sx q[3];
rz(1.1480918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-2.5941217) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(-2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9861458) q[0];
sx q[0];
rz(-2.2387945) q[0];
sx q[0];
rz(-1.6814394) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(0.11988457) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6836707) q[0];
sx q[0];
rz(-0.73212762) q[0];
sx q[0];
rz(-1.2326272) q[0];
rz(-pi) q[1];
rz(0.91953711) q[2];
sx q[2];
rz(-2.0555758) q[2];
sx q[2];
rz(-0.19249053) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6189055) q[1];
sx q[1];
rz(-2.9123016) q[1];
sx q[1];
rz(2.7875336) q[1];
rz(1.0644887) q[3];
sx q[3];
rz(-1.0303921) q[3];
sx q[3];
rz(-3.0960954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(-0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6108625) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(-2.7643381) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(2.196905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7101038) q[0];
sx q[0];
rz(-2.1444107) q[0];
sx q[0];
rz(1.3106034) q[0];
rz(0.60010054) q[2];
sx q[2];
rz(-2.1515111) q[2];
sx q[2];
rz(2.3490259) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8807043) q[1];
sx q[1];
rz(-2.7847538) q[1];
sx q[1];
rz(-2.6964705) q[1];
x q[2];
rz(1.3768436) q[3];
sx q[3];
rz(-2.5911281) q[3];
sx q[3];
rz(2.6725519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3433596) q[2];
sx q[2];
rz(-1.6075906) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-2.01912) q[3];
sx q[3];
rz(1.6759492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632904) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(-2.5569051) q[0];
rz(-1.0514642) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(0.47031602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10775751) q[0];
sx q[0];
rz(-1.4246539) q[0];
sx q[0];
rz(0.58476292) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1675179) q[2];
sx q[2];
rz(-1.2263311) q[2];
sx q[2];
rz(-0.044805077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85519771) q[1];
sx q[1];
rz(-0.87955399) q[1];
sx q[1];
rz(-1.0950086) q[1];
rz(-pi) q[2];
rz(-1.5585654) q[3];
sx q[3];
rz(-1.6837956) q[3];
sx q[3];
rz(1.6390273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8768002) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(-2.8301767) q[2];
rz(-0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(2.8581207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51157057) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(2.2054963) q[0];
rz(0.49631897) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(2.671303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966489) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(-2.6281283) q[0];
rz(-pi) q[1];
rz(1.5580721) q[2];
sx q[2];
rz(-0.1977405) q[2];
sx q[2];
rz(1.5997353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5230323) q[1];
sx q[1];
rz(-1.2292687) q[1];
sx q[1];
rz(0.073445436) q[1];
rz(-0.55996446) q[3];
sx q[3];
rz(-0.32295152) q[3];
sx q[3];
rz(1.1895869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5101667) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(-1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(0.9789595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706021) q[0];
sx q[0];
rz(-0.63828642) q[0];
sx q[0];
rz(0.19564512) q[0];
x q[1];
rz(-2.8674815) q[2];
sx q[2];
rz(-0.90702552) q[2];
sx q[2];
rz(1.889515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0461523) q[1];
sx q[1];
rz(-0.80712026) q[1];
sx q[1];
rz(0.67428204) q[1];
rz(-pi) q[2];
rz(-0.03421182) q[3];
sx q[3];
rz(-2.1497823) q[3];
sx q[3];
rz(1.2466696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0673361) q[2];
sx q[2];
rz(-0.29198519) q[2];
sx q[2];
rz(-3.0687029) q[2];
rz(-2.5439751) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(3.1387175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6947967) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(-2.6126675) q[0];
rz(-2.9534598) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(-2.5573152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408509) q[0];
sx q[0];
rz(-2.7817431) q[0];
sx q[0];
rz(0.76143439) q[0];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.9263679) q[2];
sx q[2];
rz(-0.3027161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1288695) q[1];
sx q[1];
rz(-1.1616316) q[1];
sx q[1];
rz(1.6175458) q[1];
rz(-pi) q[2];
rz(-0.27354555) q[3];
sx q[3];
rz(-2.1714032) q[3];
sx q[3];
rz(2.724444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(0.071852597) q[2];
rz(2.1182649) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(-2.6527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95579424) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(2.4664948) q[1];
sx q[1];
rz(-1.5470807) q[1];
sx q[1];
rz(3.0463228) q[1];
rz(-2.6266392) q[2];
sx q[2];
rz(-1.5025768) q[2];
sx q[2];
rz(1.2160355) q[2];
rz(-2.974398) q[3];
sx q[3];
rz(-1.7291369) q[3];
sx q[3];
rz(-0.72469934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
