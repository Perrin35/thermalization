OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(-0.5501774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30480706) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(0.027713393) q[0];
rz(-0.43843856) q[2];
sx q[2];
rz(-1.0796094) q[2];
sx q[2];
rz(0.05471281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6018511) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(-2.3439581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1816741) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.5126022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249811) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(-2.7532817) q[0];
rz(1.5674044) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-1.1018745) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.70817972) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(-1.3604926) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2643331) q[3];
sx q[3];
rz(-0.3650529) q[3];
sx q[3];
rz(0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(-2.3382323) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.196516) q[0];
sx q[0];
rz(-1.7917504) q[0];
sx q[0];
rz(-3.0483079) q[0];
x q[1];
rz(-0.27772851) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(-0.03253983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5902139) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(-2.6022634) q[1];
rz(-pi) q[2];
rz(0.48084534) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(-2.0744827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21928366) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(0.10198378) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(0.27329683) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057987) q[0];
sx q[0];
rz(-1.6596284) q[0];
sx q[0];
rz(1.2466794) q[0];
x q[1];
rz(2.4934713) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(-1.8304706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3015775) q[1];
sx q[1];
rz(-1.8043648) q[1];
sx q[1];
rz(-2.9905) q[1];
x q[2];
rz(-2.9076505) q[3];
sx q[3];
rz(-0.66286874) q[3];
sx q[3];
rz(2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-0.50393528) q[2];
rz(-3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824317) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(3.058847) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8866681) q[0];
sx q[0];
rz(-2.3912171) q[0];
sx q[0];
rz(0.95959856) q[0];
rz(-pi) q[1];
rz(1.5802025) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(-2.8807092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.134569) q[1];
sx q[1];
rz(-2.1937074) q[1];
sx q[1];
rz(-1.182215) q[1];
rz(-3.0994814) q[3];
sx q[3];
rz(-1.9633506) q[3];
sx q[3];
rz(-0.3596572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(0.791839) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.8621559) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5813542) q[0];
sx q[0];
rz(-1.8525774) q[0];
sx q[0];
rz(3.1262585) q[0];
rz(-1.2811786) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(2.9786125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.6305627) q[1];
x q[2];
rz(-2.2313314) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(-1.4753301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(3.1380222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768893) q[0];
sx q[0];
rz(-1.9386374) q[0];
sx q[0];
rz(0.44643114) q[0];
rz(2.7262906) q[2];
sx q[2];
rz(-2.3447782) q[2];
sx q[2];
rz(-1.0559527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0842972) q[1];
sx q[1];
rz(-2.1665386) q[1];
sx q[1];
rz(0.36235313) q[1];
rz(-pi) q[2];
rz(1.9332063) q[3];
sx q[3];
rz(-2.9299195) q[3];
sx q[3];
rz(3.1162804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(2.4628941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40020254) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(-2.461117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83301021) q[2];
sx q[2];
rz(-1.2150803) q[2];
sx q[2];
rz(0.21258159) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.855367) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(-1.0463868) q[1];
x q[2];
rz(-1.7292761) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-0.93723047) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(2.881799) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(-2.7816539) q[0];
rz(0.91295816) q[2];
sx q[2];
rz(-1.3557938) q[2];
sx q[2];
rz(-2.1195597) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.305294) q[1];
sx q[1];
rz(0.32151476) q[1];
x q[2];
rz(2.6960877) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(1.6642237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-0.3607761) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11689582) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(-0.025339729) q[0];
rz(3.0214494) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(0.10211589) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1211179) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(-2.3974182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9721857) q[3];
sx q[3];
rz(-1.425404) q[3];
sx q[3];
rz(1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5671134) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(-1.3600596) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-0.097638771) q[2];
sx q[2];
rz(-1.3664445) q[2];
sx q[2];
rz(0.22899361) q[2];
rz(-2.4102224) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
