OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(-1.9089729) q[1];
sx q[1];
rz(0.90484172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.605699) q[0];
sx q[0];
rz(-0.68471691) q[0];
sx q[0];
rz(0.84795714) q[0];
rz(-0.22613871) q[2];
sx q[2];
rz(-1.3790501) q[2];
sx q[2];
rz(-2.6297671) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20476725) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(2.8863465) q[1];
x q[2];
rz(2.1078029) q[3];
sx q[3];
rz(-0.35764965) q[3];
sx q[3];
rz(-1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-0.15371418) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(2.4786425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6843296) q[0];
sx q[0];
rz(-1.508679) q[0];
sx q[0];
rz(-1.6371884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7402678) q[2];
sx q[2];
rz(-2.1839645) q[2];
sx q[2];
rz(-0.80712986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6389097) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(1.6619976) q[1];
rz(1.4088267) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(-1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0654304) q[0];
sx q[0];
rz(-2.2632474) q[0];
sx q[0];
rz(-2.6718219) q[0];
x q[1];
rz(-2.8754183) q[2];
sx q[2];
rz(-1.4783579) q[2];
sx q[2];
rz(-3.0150974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7171214) q[1];
sx q[1];
rz(-0.79454225) q[1];
sx q[1];
rz(-2.8984927) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(-0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(3.0730491) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-3.0025735) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67447353) q[0];
sx q[0];
rz(-1.6498483) q[0];
sx q[0];
rz(2.5208958) q[0];
rz(0.94167534) q[2];
sx q[2];
rz(-0.86949124) q[2];
sx q[2];
rz(-0.47052449) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9070223) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(2.6837818) q[1];
x q[2];
rz(-1.7771878) q[3];
sx q[3];
rz(-1.5406113) q[3];
sx q[3];
rz(-0.15011945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-2.9117057) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(2.5255323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992633) q[0];
sx q[0];
rz(-3.0780601) q[0];
sx q[0];
rz(-0.97139831) q[0];
rz(-2.8370503) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(-2.9074557) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7459813) q[1];
sx q[1];
rz(-0.99318722) q[1];
sx q[1];
rz(1.5773768) q[1];
rz(2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(1.5622996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(0.88476673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038267604) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(0.079770712) q[0];
rz(-pi) q[1];
rz(-0.21005819) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(-2.7163598) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6073508) q[1];
sx q[1];
rz(-0.51041767) q[1];
sx q[1];
rz(0.95153248) q[1];
rz(2.8956036) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-0.3113783) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-2.4087002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830025) q[0];
sx q[0];
rz(-0.68651474) q[0];
sx q[0];
rz(0.65450432) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1204719) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(2.6598425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1454754) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(-1.5596418) q[1];
rz(-pi) q[2];
rz(1.6781428) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(2.6220208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(2.627009) q[2];
rz(-1.9040646) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(0.27871305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7473135) q[0];
sx q[0];
rz(-0.25082591) q[0];
sx q[0];
rz(-0.5745116) q[0];
x q[1];
rz(-0.5586241) q[2];
sx q[2];
rz(-1.363029) q[2];
sx q[2];
rz(-0.34780207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7292494) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(1.8885683) q[1];
rz(-3.0635733) q[3];
sx q[3];
rz(-1.3481513) q[3];
sx q[3];
rz(1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.1466325) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.1425346) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.8737428) q[0];
x q[1];
rz(-2.4382486) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(1.6714931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7361703) q[1];
sx q[1];
rz(-1.0338963) q[1];
sx q[1];
rz(-1.7169397) q[1];
rz(-pi) q[2];
rz(1.3167131) q[3];
sx q[3];
rz(-0.68613201) q[3];
sx q[3];
rz(0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.9030301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68574821) q[0];
sx q[0];
rz(-2.0861174) q[0];
sx q[0];
rz(2.8932163) q[0];
rz(-pi) q[1];
rz(0.752991) q[2];
sx q[2];
rz(-1.6082616) q[2];
sx q[2];
rz(2.3245036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8348332) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(1.761761) q[1];
rz(0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(-0.30944165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(3.0977541) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6702406) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(2.2123443) q[2];
sx q[2];
rz(-0.48454787) q[2];
sx q[2];
rz(1.5018644) q[2];
rz(2.9712215) q[3];
sx q[3];
rz(-0.80633612) q[3];
sx q[3];
rz(2.7663305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
