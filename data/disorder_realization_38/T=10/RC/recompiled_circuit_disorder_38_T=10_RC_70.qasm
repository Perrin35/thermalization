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
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2599517) q[0];
sx q[0];
rz(-1.5437484) q[0];
sx q[0];
rz(-1.7904439) q[0];
x q[1];
rz(-0.89994853) q[2];
sx q[2];
rz(-2.4953825) q[2];
sx q[2];
rz(-0.83713573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.9508233) q[1];
sx q[1];
rz(-0.7976346) q[1];
rz(-1.4518277) q[3];
sx q[3];
rz(-0.61421466) q[3];
sx q[3];
rz(0.14106942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-0.63981167) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249811) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(2.7532817) q[0];
x q[1];
rz(-0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-0.46797215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1158893) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(-2.2504836) q[1];
rz(-1.8566425) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(2.3957516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5416109) q[0];
sx q[0];
rz(-0.2395425) q[0];
sx q[0];
rz(1.177686) q[0];
rz(-pi) q[1];
rz(2.6150871) q[2];
sx q[2];
rz(-1.7130934) q[2];
sx q[2];
rz(-1.7775747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5513788) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(0.53932921) q[1];
x q[2];
rz(-0.45332076) q[3];
sx q[3];
rz(-1.7954149) q[3];
sx q[3];
rz(2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-0.050343242) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(1.5215727) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.2980952) q[0];
rz(0.64812135) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(1.8304706) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30444333) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(-1.8069581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9076505) q[3];
sx q[3];
rz(-0.66286874) q[3];
sx q[3];
rz(-2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-0.50393528) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6325729) q[0];
sx q[0];
rz(-0.97826725) q[0];
sx q[0];
rz(-2.65034) q[0];
x q[1];
rz(0.012533112) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(-2.8963793) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0070237006) q[1];
sx q[1];
rz(-2.1937074) q[1];
sx q[1];
rz(1.182215) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9636642) q[3];
sx q[3];
rz(-1.6097027) q[3];
sx q[3];
rz(1.9143357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8905028) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.8621559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148221) q[0];
sx q[0];
rz(-1.5855256) q[0];
sx q[0];
rz(1.8526088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7369667) q[2];
sx q[2];
rz(-2.4936516) q[2];
sx q[2];
rz(-2.4857156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64843236) q[1];
sx q[1];
rz(-1.5409768) q[1];
sx q[1];
rz(-1.0489419) q[1];
x q[2];
rz(3.0002954) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(2.3395204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-0.77077579) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-0.055667002) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-0.0035704426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(0.44643114) q[0];
rz(-1.1793169) q[2];
sx q[2];
rz(-2.2841095) q[2];
sx q[2];
rz(-1.6187402) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7230941) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(2.198092) q[1];
rz(-1.3725029) q[3];
sx q[3];
rz(-1.4962422) q[3];
sx q[3];
rz(-1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5355276) q[0];
sx q[0];
rz(-1.3074271) q[0];
sx q[0];
rz(0.33959629) q[0];
x q[1];
rz(-2.0754201) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(-1.417516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39290909) q[1];
sx q[1];
rz(-2.6137685) q[1];
sx q[1];
rz(1.6960359) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4123165) q[3];
sx q[3];
rz(-1.1500119) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(2.881799) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-0.52694595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5848815) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(-0.35993872) q[0];
rz(-1.9138463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(0.81817852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.305294) q[1];
sx q[1];
rz(0.32151476) q[1];
x q[2];
rz(1.280904) q[3];
sx q[3];
rz(-1.999951) q[3];
sx q[3];
rz(0.030062519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(0.33716831) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11689582) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(0.025339729) q[0];
rz(1.0656409) q[2];
sx q[2];
rz(-0.2444707) q[2];
sx q[2];
rz(-0.62015647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4611778) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(1.7635036) q[1];
rz(-pi) q[2];
rz(1.4233227) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(-2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.7452314) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(0.097638771) q[2];
sx q[2];
rz(-1.7751481) q[2];
sx q[2];
rz(-2.912599) q[2];
rz(0.75136649) q[3];
sx q[3];
rz(-1.02117) q[3];
sx q[3];
rz(1.1985967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];