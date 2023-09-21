OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(3.5482121) q[0];
sx q[0];
rz(9.1756048) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(4.1133147) q[1];
sx q[1];
rz(9.9749554) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8367856) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(-3.1138793) q[0];
rz(-2.7031541) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(-3.0868798) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6018511) q[1];
sx q[1];
rz(-1.9508233) q[1];
sx q[1];
rz(-2.3439581) q[1];
rz(-pi) q[2];
rz(1.4518277) q[3];
sx q[3];
rz(-2.527378) q[3];
sx q[3];
rz(0.14106942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249811) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(-2.7532817) q[0];
x q[1];
rz(1.5741882) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-2.0397182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0257033) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(2.2504836) q[1];
rz(-pi) q[2];
rz(0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2597044) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(2.3957516) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5999818) q[0];
sx q[0];
rz(-2.9020502) q[0];
sx q[0];
rz(-1.9639067) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8638641) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(0.03253983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5513788) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(-2.6022634) q[1];
rz(-0.45332076) q[3];
sx q[3];
rz(-1.7954149) q[3];
sx q[3];
rz(2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-2.8059373) q[0];
sx q[0];
rz(1.2980952) q[0];
rz(-pi) q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8371493) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(-1.3346346) q[1];
rz(-2.4920739) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(2.3424145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(-2.1544429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6325729) q[0];
sx q[0];
rz(-2.1633254) q[0];
sx q[0];
rz(-0.49125262) q[0];
rz(-1.5802025) q[2];
sx q[2];
rz(-2.214589) q[2];
sx q[2];
rz(-2.8807092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.8835888) q[1];
sx q[1];
rz(0.66004628) q[1];
rz(-pi) q[2];
rz(-1.6721252) q[3];
sx q[3];
rz(-0.39468995) q[3];
sx q[3];
rz(0.24995382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-2.3497537) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(-1.2794367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6364481) q[0];
sx q[0];
rz(-0.2821869) q[0];
sx q[0];
rz(-1.6237153) q[0];
x q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(-2.9786125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2020859) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(-0.034394666) q[1];
rz(-pi) q[2];
rz(0.1412973) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768893) q[0];
sx q[0];
rz(-1.9386374) q[0];
sx q[0];
rz(-2.6951615) q[0];
rz(-pi) q[1];
rz(0.41530208) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(-1.0559527) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6047302) q[1];
sx q[1];
rz(-0.68568789) q[1];
sx q[1];
rz(-1.0889978) q[1];
rz(-pi) q[2];
rz(1.3725029) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(-1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.7238808) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6060651) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(-2.8019964) q[0];
rz(-pi) q[1];
rz(2.0754201) q[2];
sx q[2];
rz(-0.80427158) q[2];
sx q[2];
rz(1.7240766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8933967) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(-3.0688973) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33903665) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(2.8807004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(2.6146467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(-2.7816539) q[0];
rz(-pi) q[1];
rz(0.91295816) q[2];
sx q[2];
rz(-1.3557938) q[2];
sx q[2];
rz(1.022033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37104169) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(0.71055926) q[1];
rz(-2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4496778) q[0];
sx q[0];
rz(-1.5458108) q[0];
sx q[0];
rz(-1.4033532) q[0];
rz(-pi) q[1];
rz(-1.0656409) q[2];
sx q[2];
rz(-2.897122) q[2];
sx q[2];
rz(2.5214362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39721397) q[1];
sx q[1];
rz(-0.75512868) q[1];
sx q[1];
rz(0.20882864) q[1];
rz(-pi) q[2];
rz(-2.9721857) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(1.3225079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5671134) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-1.3654937) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(2.4102224) q[3];
sx q[3];
rz(-2.2435355) q[3];
sx q[3];
rz(0.13767195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];