OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1448016) q[0];
sx q[0];
rz(0.15455833) q[0];
sx q[0];
rz(6.9757087) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614852) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(1.1397584) q[0];
x q[1];
rz(-0.55180438) q[2];
sx q[2];
rz(-1.8066415) q[2];
sx q[2];
rz(0.15197309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4298809) q[1];
sx q[1];
rz(-1.8036588) q[1];
sx q[1];
rz(-1.4038605) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65269835) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(2.936426) q[2];
rz(-2.3702879) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(2.0390959) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40760621) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(-1.0247963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.227238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7176712) q[0];
sx q[0];
rz(-2.9917891) q[0];
sx q[0];
rz(2.1013837) q[0];
rz(-pi) q[1];
rz(0.51867698) q[2];
sx q[2];
rz(-1.1653324) q[2];
sx q[2];
rz(-0.60207089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2482359) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(2.8398819) q[1];
x q[2];
rz(-1.4144054) q[3];
sx q[3];
rz(-1.3723433) q[3];
sx q[3];
rz(0.95201492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(0.36519095) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(-0.99575106) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(2.8083037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5211398) q[0];
sx q[0];
rz(-2.1846909) q[0];
sx q[0];
rz(-0.99434538) q[0];
x q[1];
rz(-0.097924175) q[2];
sx q[2];
rz(-1.7469179) q[2];
sx q[2];
rz(-0.91913659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7283199) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(2.1876213) q[1];
rz(-2.1214478) q[3];
sx q[3];
rz(-1.2216976) q[3];
sx q[3];
rz(0.89494866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(-1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(-0.68471318) q[0];
rz(-1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.9365786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689777) q[0];
sx q[0];
rz(-2.2616771) q[0];
sx q[0];
rz(-0.0432424) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8393458) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(0.72999398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4971784) q[1];
sx q[1];
rz(-0.89343151) q[1];
sx q[1];
rz(0.35269423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2756696) q[3];
sx q[3];
rz(-1.5161627) q[3];
sx q[3];
rz(1.6161402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(3.0084685) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(-0.55508074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984891) q[0];
sx q[0];
rz(-1.886133) q[0];
sx q[0];
rz(-0.01339162) q[0];
x q[1];
rz(-1.7246036) q[2];
sx q[2];
rz(-0.4193192) q[2];
sx q[2];
rz(-2.9749982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2808387) q[1];
sx q[1];
rz(-2.6959531) q[1];
sx q[1];
rz(1.526236) q[1];
rz(-0.080901905) q[3];
sx q[3];
rz(-1.6228075) q[3];
sx q[3];
rz(-2.4018283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(-2.1991918) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(0.55148235) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72551661) q[0];
sx q[0];
rz(-2.2726739) q[0];
sx q[0];
rz(-0.26728018) q[0];
rz(-pi) q[1];
rz(-0.13055735) q[2];
sx q[2];
rz(-1.6160384) q[2];
sx q[2];
rz(-0.10938489) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3547937) q[1];
sx q[1];
rz(-1.2878294) q[1];
sx q[1];
rz(-1.1250886) q[1];
rz(-0.88605373) q[3];
sx q[3];
rz(-1.3430809) q[3];
sx q[3];
rz(-2.7453604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55398983) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(0.51876846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9414026) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(0.83394136) q[0];
rz(2.9154645) q[2];
sx q[2];
rz(-0.78352189) q[2];
sx q[2];
rz(-2.4353611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7873951) q[1];
sx q[1];
rz(-2.0320315) q[1];
sx q[1];
rz(-2.4128777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026168907) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(-2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(2.5781412) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.7425591) q[0];
rz(2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-0.74434892) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8101013) q[0];
sx q[0];
rz(-1.6252675) q[0];
sx q[0];
rz(-0.14793747) q[0];
rz(-pi) q[1];
rz(1.3049576) q[2];
sx q[2];
rz(-0.53005866) q[2];
sx q[2];
rz(-2.838138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8213615) q[1];
sx q[1];
rz(-1.3470955) q[1];
sx q[1];
rz(-0.51364586) q[1];
rz(0.45458557) q[3];
sx q[3];
rz(-1.5496407) q[3];
sx q[3];
rz(1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4259592) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(2.5320833) q[2];
rz(2.4842747) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(2.8801584) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(-1.9001182) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(2.3666568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938213) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(2.2872778) q[0];
x q[1];
rz(-1.2504134) q[2];
sx q[2];
rz(-0.46386007) q[2];
sx q[2];
rz(1.8360209) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1800268) q[1];
sx q[1];
rz(-2.0225836) q[1];
sx q[1];
rz(1.7801442) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3516515) q[3];
sx q[3];
rz(-1.7100167) q[3];
sx q[3];
rz(-2.5185891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(-2.9837218) q[2];
rz(-1.212451) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(2.4841323) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(2.0956031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5770618) q[0];
sx q[0];
rz(-0.37527592) q[0];
sx q[0];
rz(-0.29348404) q[0];
x q[1];
rz(-2.5506053) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(-1.9901333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54088456) q[1];
sx q[1];
rz(-0.67544671) q[1];
sx q[1];
rz(0.9687959) q[1];
rz(0.06185992) q[3];
sx q[3];
rz(-0.40611551) q[3];
sx q[3];
rz(-1.5861685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(-1.0894758) q[2];
rz(-1.5661092) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2789223) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.6090341) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(0.63411843) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(3.0388721) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
