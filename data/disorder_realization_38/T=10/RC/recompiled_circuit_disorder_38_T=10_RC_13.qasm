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
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(-0.5501774) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.881641) q[0];
sx q[0];
rz(-1.5437484) q[0];
sx q[0];
rz(1.3511488) q[0];
rz(-0.43843856) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(3.0868798) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5397415) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(2.3439581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1816741) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(-2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0523895) q[0];
sx q[0];
rz(-0.69141483) q[0];
sx q[0];
rz(-1.0538488) q[0];
rz(1.5674044) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-1.1018745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43803793) q[1];
sx q[1];
rz(-2.2433271) q[1];
sx q[1];
rz(-2.9707675) q[1];
rz(-0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(2.3957516) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999818) q[0];
sx q[0];
rz(-0.2395425) q[0];
sx q[0];
rz(1.177686) q[0];
x q[1];
rz(-2.8638641) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(-0.03253983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97400372) q[1];
sx q[1];
rz(-2.0473192) q[1];
sx q[1];
rz(-1.2567026) q[1];
rz(-pi) q[2];
rz(2.6607473) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(-1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8400152) q[1];
sx q[1];
rz(-1.8043648) q[1];
sx q[1];
rz(0.15109269) q[1];
rz(-pi) q[2];
rz(2.4920739) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(-2.3424145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(-0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-0.082745634) q[0];
rz(-0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894831) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(2.2228918) q[0];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(2.5226749) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(0.48536761) q[1];
rz(-pi) q[2];
rz(-1.4694674) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(-2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(0.791839) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.2794367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1267705) q[0];
sx q[0];
rz(-1.5855256) q[0];
sx q[0];
rz(-1.8526088) q[0];
rz(-0.40462599) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(-0.65587703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2710323) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.51103) q[1];
x q[2];
rz(-3.0002954) q[3];
sx q[3];
rz(-1.391489) q[3];
sx q[3];
rz(2.3395204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(0.77077579) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-0.055667002) q[0];
rz(2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(2.4128782) q[0];
rz(-2.3891719) q[2];
sx q[2];
rz(-1.2781029) q[2];
sx q[2];
rz(0.21586403) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6047302) q[1];
sx q[1];
rz(-0.68568789) q[1];
sx q[1];
rz(2.0525949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2083863) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(-0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(1.8493269) q[0];
x q[1];
rz(-2.6762814) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(-1.0516143) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24819599) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(0.072695331) q[1];
rz(2.802556) q[3];
sx q[3];
rz(-0.44796523) q[3];
sx q[3];
rz(0.26089222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53370595) q[0];
sx q[0];
rz(-2.7237646) q[0];
sx q[0];
rz(2.582344) q[0];
x q[1];
rz(-1.9138463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(2.3234141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.867013) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(0.32151476) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6960877) q[3];
sx q[3];
rz(-1.3078491) q[3];
sx q[3];
rz(-1.477369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(2.7808166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6919149) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(-1.4033532) q[0];
x q[1];
rz(3.0214494) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(-3.0394768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39721397) q[1];
sx q[1];
rz(-2.386464) q[1];
sx q[1];
rz(0.20882864) q[1];
rz(1.4233227) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(-1.7815331) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(-1.1311244) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(-0.87313575) q[3];
sx q[3];
rz(-2.1921039) q[3];
sx q[3];
rz(2.3154142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];