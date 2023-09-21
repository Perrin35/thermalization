OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(-2.7116099) q[1];
sx q[1];
rz(-2.4584682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8389991) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(-1.6291314) q[0];
rz(-pi) q[1];
rz(0.97857742) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(0.65939553) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4343623) q[1];
sx q[1];
rz(-2.011236) q[1];
sx q[1];
rz(2.4643154) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1300415) q[3];
sx q[3];
rz(-0.64239255) q[3];
sx q[3];
rz(-1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10509051) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1128088) q[2];
sx q[2];
rz(-0.53456842) q[2];
sx q[2];
rz(-2.2528258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1190471) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(-0.70336282) q[1];
x q[2];
rz(-2.7705455) q[3];
sx q[3];
rz(-2.3608532) q[3];
sx q[3];
rz(-2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96317545) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8974509) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.1600526) q[0];
rz(-pi) q[1];
rz(-0.66652253) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(-0.18596622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.264818) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-3.1289711) q[1];
x q[2];
rz(0.3018467) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(-0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(-0.96757403) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(-0.46491369) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(1.0850614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5269055) q[0];
sx q[0];
rz(-2.5589716) q[0];
sx q[0];
rz(-2.7132062) q[0];
rz(-pi) q[1];
rz(0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(2.0476066) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81740582) q[1];
sx q[1];
rz(-1.4293912) q[1];
sx q[1];
rz(-1.4694587) q[1];
x q[2];
rz(-2.5303909) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(-1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-2.1599105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22589332) q[0];
sx q[0];
rz(-2.3908273) q[0];
sx q[0];
rz(2.4322926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1443411) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(-0.85948247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9202068) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(1.7358001) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97335191) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(-0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.107347) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(-0.29883595) q[0];
rz(-1.500962) q[2];
sx q[2];
rz(-2.2572821) q[2];
sx q[2];
rz(-1.3697461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11152553) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(1.0982606) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1092471) q[3];
sx q[3];
rz(-0.49406067) q[3];
sx q[3];
rz(0.1915313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2281987) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(2.4051037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.526028) q[0];
sx q[0];
rz(-0.051858735) q[0];
sx q[0];
rz(1.2349013) q[0];
rz(1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(2.2155857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0349717) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(0.76645318) q[1];
rz(-pi) q[2];
rz(2.001708) q[3];
sx q[3];
rz(-2.0322554) q[3];
sx q[3];
rz(-1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066276) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(-pi) q[1];
rz(0.83306649) q[2];
sx q[2];
rz(-1.7852011) q[2];
sx q[2];
rz(0.094878541) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39107716) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(-0.68117546) q[1];
rz(-0.86254085) q[3];
sx q[3];
rz(-1.7077912) q[3];
sx q[3];
rz(0.18602895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.4484891) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1595759) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(-0.52635877) q[0];
x q[1];
rz(0.18025132) q[2];
sx q[2];
rz(-2.2580574) q[2];
sx q[2];
rz(2.7625411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19787994) q[1];
sx q[1];
rz(-0.76849741) q[1];
sx q[1];
rz(-2.6430623) q[1];
rz(-2.4217442) q[3];
sx q[3];
rz(-2.3206629) q[3];
sx q[3];
rz(-1.5035226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(-0.71406281) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6274174) q[0];
sx q[0];
rz(-1.6192993) q[0];
sx q[0];
rz(1.4072627) q[0];
x q[1];
rz(-2.5942957) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(1.0123569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51917045) q[1];
sx q[1];
rz(-1.6455489) q[1];
sx q[1];
rz(-2.1291817) q[1];
rz(-pi) q[2];
x q[2];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7174299) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(2.4049315) q[2];
sx q[2];
rz(-1.3623289) q[2];
sx q[2];
rz(1.3883432) q[2];
rz(-0.35076326) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
