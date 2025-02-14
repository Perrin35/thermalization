OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6629163) q[0];
sx q[0];
rz(-1.5530246) q[0];
sx q[0];
rz(-1.8832062) q[0];
rz(-5.0212669) q[1];
sx q[1];
rz(6.8016383) q[1];
sx q[1];
rz(6.3432884) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4661048) q[0];
sx q[0];
rz(-1.6610613) q[0];
sx q[0];
rz(1.5399) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8608039) q[2];
sx q[2];
rz(-0.5092237) q[2];
sx q[2];
rz(-1.4851242) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9277674) q[1];
sx q[1];
rz(-2.7672184) q[1];
sx q[1];
rz(1.8273749) q[1];
rz(-pi) q[2];
rz(-2.1677371) q[3];
sx q[3];
rz(-1.1789448) q[3];
sx q[3];
rz(2.4490956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1897159) q[2];
sx q[2];
rz(-0.78236255) q[2];
sx q[2];
rz(2.961109) q[2];
rz(-2.8372724) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75905281) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(-0.10391129) q[0];
rz(-0.78481627) q[1];
sx q[1];
rz(-1.8321313) q[1];
sx q[1];
rz(-2.5596502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3510785) q[0];
sx q[0];
rz(-1.0389555) q[0];
sx q[0];
rz(1.5449338) q[0];
x q[1];
rz(2.5594219) q[2];
sx q[2];
rz(-2.5522531) q[2];
sx q[2];
rz(-1.8631528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84060639) q[1];
sx q[1];
rz(-2.0475629) q[1];
sx q[1];
rz(-0.029164) q[1];
x q[2];
rz(-0.41116233) q[3];
sx q[3];
rz(-1.9389279) q[3];
sx q[3];
rz(2.1187834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.71501422) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(-1.0283872) q[2];
rz(2.5905124) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-2.6330131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.5904163) q[0];
sx q[0];
rz(-1.506378) q[0];
sx q[0];
rz(-1.8835541) q[0];
rz(-2.538077) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(0.18361941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1574137) q[0];
sx q[0];
rz(-1.3065728) q[0];
sx q[0];
rz(1.9963032) q[0];
rz(2.8894561) q[2];
sx q[2];
rz(-2.647433) q[2];
sx q[2];
rz(1.8273938) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10079846) q[1];
sx q[1];
rz(-1.9996627) q[1];
sx q[1];
rz(-0.78227346) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9345064) q[3];
sx q[3];
rz(-1.6112176) q[3];
sx q[3];
rz(0.8881027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6151578) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(0.60205013) q[2];
rz(0.43271068) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(-1.6656779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3848307) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(2.5153644) q[0];
rz(1.2681883) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(0.38571206) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.884366) q[0];
sx q[0];
rz(-0.8943087) q[0];
sx q[0];
rz(-2.1348597) q[0];
rz(0.63718225) q[2];
sx q[2];
rz(-0.68693012) q[2];
sx q[2];
rz(-1.0174583) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0929132) q[1];
sx q[1];
rz(-1.6374705) q[1];
sx q[1];
rz(0.20079048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1139051) q[3];
sx q[3];
rz(-1.466179) q[3];
sx q[3];
rz(-2.0547158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6993774) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(-2.8507612) q[2];
rz(1.95131) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-1.0440089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6749343) q[0];
sx q[0];
rz(-3.0789154) q[0];
sx q[0];
rz(1.689893) q[0];
rz(2.2845204) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(2.7136386) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65954094) q[0];
sx q[0];
rz(-0.45705802) q[0];
sx q[0];
rz(-1.2808475) q[0];
rz(-pi) q[1];
rz(0.73889795) q[2];
sx q[2];
rz(-1.9618926) q[2];
sx q[2];
rz(2.1192239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8331205) q[1];
sx q[1];
rz(-0.80128925) q[1];
sx q[1];
rz(0.50442969) q[1];
rz(1.0048423) q[3];
sx q[3];
rz(-2.1814846) q[3];
sx q[3];
rz(2.7683059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2032623) q[2];
sx q[2];
rz(-1.0469971) q[2];
sx q[2];
rz(0.14673512) q[2];
rz(0.19715582) q[3];
sx q[3];
rz(-1.4898172) q[3];
sx q[3];
rz(0.76876918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305304) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(-0.16832571) q[0];
rz(2.7492145) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.6900774) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096286721) q[0];
sx q[0];
rz(-1.1340965) q[0];
sx q[0];
rz(-0.023604579) q[0];
rz(-pi) q[1];
rz(-0.2608725) q[2];
sx q[2];
rz(-2.2745273) q[2];
sx q[2];
rz(1.0359302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2830303) q[1];
sx q[1];
rz(-1.7114189) q[1];
sx q[1];
rz(-0.750641) q[1];
x q[2];
rz(2.0198972) q[3];
sx q[3];
rz(-2.3074352) q[3];
sx q[3];
rz(0.68461217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3671941) q[2];
sx q[2];
rz(-0.44782475) q[2];
sx q[2];
rz(1.9642584) q[2];
rz(-2.2199953) q[3];
sx q[3];
rz(-0.84601837) q[3];
sx q[3];
rz(-2.6109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717936) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(1.0569093) q[0];
rz(-0.22843703) q[1];
sx q[1];
rz(-2.3616796) q[1];
sx q[1];
rz(2.8905919) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91755501) q[0];
sx q[0];
rz(-2.7537808) q[0];
sx q[0];
rz(2.4075131) q[0];
rz(2.2784581) q[2];
sx q[2];
rz(-0.8388817) q[2];
sx q[2];
rz(-1.9307435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7556762) q[1];
sx q[1];
rz(-1.7337013) q[1];
sx q[1];
rz(-1.7079123) q[1];
rz(1.1993221) q[3];
sx q[3];
rz(-2.0544996) q[3];
sx q[3];
rz(2.7886645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2061657) q[2];
sx q[2];
rz(-0.82038227) q[2];
sx q[2];
rz(0.41681918) q[2];
rz(2.5721512) q[3];
sx q[3];
rz(-1.1007525) q[3];
sx q[3];
rz(2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5702629) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(2.6369693) q[0];
rz(2.8847671) q[1];
sx q[1];
rz(-2.3378614) q[1];
sx q[1];
rz(1.1508734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1543836) q[0];
sx q[0];
rz(-1.0797653) q[0];
sx q[0];
rz(2.9084297) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0382284) q[2];
sx q[2];
rz(-0.69004493) q[2];
sx q[2];
rz(1.104081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.098365) q[1];
sx q[1];
rz(-2.4154764) q[1];
sx q[1];
rz(2.4642706) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34567771) q[3];
sx q[3];
rz(-0.44000235) q[3];
sx q[3];
rz(2.1328762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5095832) q[2];
sx q[2];
rz(-0.42319599) q[2];
sx q[2];
rz(0.43935856) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(-1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87089649) q[0];
sx q[0];
rz(-2.3211711) q[0];
sx q[0];
rz(2.6743555) q[0];
rz(-2.1029419) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(1.4899303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56541967) q[0];
sx q[0];
rz(-1.7459946) q[0];
sx q[0];
rz(2.6433667) q[0];
rz(0.14665276) q[2];
sx q[2];
rz(-1.4770368) q[2];
sx q[2];
rz(-0.74945005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0831775) q[1];
sx q[1];
rz(-1.5442368) q[1];
sx q[1];
rz(-0.81574622) q[1];
rz(-pi) q[2];
rz(-0.64405264) q[3];
sx q[3];
rz(-0.67094147) q[3];
sx q[3];
rz(1.80598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6963639) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(2.1571531) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-0.84243542) q[0];
sx q[0];
rz(0.62826759) q[0];
rz(2.9382622) q[1];
sx q[1];
rz(-2.0275828) q[1];
sx q[1];
rz(-0.59439739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3056476) q[0];
sx q[0];
rz(-1.3653269) q[0];
sx q[0];
rz(2.1688658) q[0];
x q[1];
rz(-0.35697414) q[2];
sx q[2];
rz(-0.56783119) q[2];
sx q[2];
rz(-2.6326376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17972936) q[1];
sx q[1];
rz(-2.7309901) q[1];
sx q[1];
rz(-0.73335464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6420308) q[3];
sx q[3];
rz(-2.1826577) q[3];
sx q[3];
rz(-2.9737986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9336885) q[2];
sx q[2];
rz(-2.4173357) q[2];
sx q[2];
rz(-2.9653463) q[2];
rz(-1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(-0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9888196) q[0];
sx q[0];
rz(-1.4403227) q[0];
sx q[0];
rz(1.3829917) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(2.9043612) q[2];
sx q[2];
rz(-1.6085515) q[2];
sx q[2];
rz(-1.5534437) q[2];
rz(-0.38269855) q[3];
sx q[3];
rz(-2.0321587) q[3];
sx q[3];
rz(1.7176499) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
