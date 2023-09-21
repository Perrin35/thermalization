OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21539772) q[0];
sx q[0];
rz(-2.8370259) q[0];
sx q[0];
rz(2.347441) q[0];
rz(0.60249451) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(0.22533016) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029763873) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(-1.8087216) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5570251) q[3];
sx q[3];
rz(-2.3574986) q[3];
sx q[3];
rz(-2.9053094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(0.72201133) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(-0.66501578) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-0.87759334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7769593) q[0];
sx q[0];
rz(-1.7517125) q[0];
sx q[0];
rz(-2.4448256) q[0];
x q[1];
rz(-1.7052824) q[2];
sx q[2];
rz(-1.209139) q[2];
sx q[2];
rz(2.550617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3493537) q[1];
sx q[1];
rz(-1.6041479) q[1];
sx q[1];
rz(-1.6778212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9545994) q[3];
sx q[3];
rz(-1.882453) q[3];
sx q[3];
rz(-0.66031885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.6931504) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(3.0677632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89798112) q[0];
sx q[0];
rz(-1.5092761) q[0];
sx q[0];
rz(3.0520526) q[0];
rz(-pi) q[1];
rz(-0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.8287303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.044588305) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(-1.3929277) q[1];
rz(2.79014) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6510216) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(1.1086639) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(-1.1289319) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(-1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98722092) q[0];
sx q[0];
rz(-2.4043596) q[0];
sx q[0];
rz(3.0043976) q[0];
rz(-pi) q[1];
rz(0.75886274) q[2];
sx q[2];
rz(-0.43735158) q[2];
sx q[2];
rz(-2.8749089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97836271) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(-2.7510838) q[1];
x q[2];
rz(3.0152263) q[3];
sx q[3];
rz(-1.4411981) q[3];
sx q[3];
rz(1.6882997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554768) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(-1.0771846) q[0];
rz(2.8793648) q[2];
sx q[2];
rz(-1.3272459) q[2];
sx q[2];
rz(-1.7521996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74354467) q[1];
sx q[1];
rz(-3.0882356) q[1];
sx q[1];
rz(1.3392901) q[1];
rz(2.3718194) q[3];
sx q[3];
rz(-2.4379745) q[3];
sx q[3];
rz(2.9911656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(2.3357847) q[2];
rz(-2.6082883) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-2.1300952) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(-1.3202753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4754007) q[0];
sx q[0];
rz(-1.1328508) q[0];
sx q[0];
rz(-1.9802718) q[0];
rz(-pi) q[1];
rz(-2.5713021) q[2];
sx q[2];
rz(-1.7208793) q[2];
sx q[2];
rz(-0.97638408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4016061) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(2.2120038) q[1];
rz(-1.9404066) q[3];
sx q[3];
rz(-0.69023057) q[3];
sx q[3];
rz(1.6604916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(2.5406204) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.8619974) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(2.5860508) q[0];
rz(3.1069966) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.3909891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23404183) q[0];
sx q[0];
rz(-0.64767917) q[0];
sx q[0];
rz(-0.56869047) q[0];
rz(0.58629845) q[2];
sx q[2];
rz(-2.7099897) q[2];
sx q[2];
rz(-1.4177711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(0.89819737) q[1];
rz(-0.075918003) q[3];
sx q[3];
rz(-2.3644991) q[3];
sx q[3];
rz(2.0457552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(2.18816) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.7038201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9490818) q[0];
sx q[0];
rz(-1.7532945) q[0];
sx q[0];
rz(-2.0927246) q[0];
x q[1];
rz(-1.4346801) q[2];
sx q[2];
rz(-1.4021177) q[2];
sx q[2];
rz(-0.31751925) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2807323) q[1];
sx q[1];
rz(-0.8287462) q[1];
sx q[1];
rz(-0.6764722) q[1];
rz(-pi) q[2];
rz(1.1233166) q[3];
sx q[3];
rz(-1.5527225) q[3];
sx q[3];
rz(-1.7759089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20539595) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(-2.4049092) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(-2.0827983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1595488) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(0.3190785) q[0];
rz(0.40760298) q[2];
sx q[2];
rz(-2.0098445) q[2];
sx q[2];
rz(1.1328732) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.345929) q[1];
sx q[1];
rz(-1.9469896) q[1];
sx q[1];
rz(-0.13161195) q[1];
x q[2];
rz(1.9731673) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8686691) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(-1.6607364) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2925443) q[0];
sx q[0];
rz(-1.2546854) q[0];
sx q[0];
rz(-2.6228117) q[0];
rz(-pi) q[1];
rz(-2.6860793) q[2];
sx q[2];
rz(-1.7835155) q[2];
sx q[2];
rz(0.46609512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4869625) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(-2.589588) q[1];
rz(1.3833952) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(2.8965829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6802406) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(1.6067778) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(-0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5948982) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(2.7728511) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(-2.4070807) q[2];
sx q[2];
rz(-2.8399158) q[2];
sx q[2];
rz(0.5010571) q[2];
rz(-1.3737804) q[3];
sx q[3];
rz(-1.1548629) q[3];
sx q[3];
rz(-2.9668273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
