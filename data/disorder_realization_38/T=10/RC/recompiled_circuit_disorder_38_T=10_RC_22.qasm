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
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(-2.1698706) q[1];
sx q[1];
rz(0.5501774) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2599517) q[0];
sx q[0];
rz(-1.5437484) q[0];
sx q[0];
rz(1.3511488) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0371738) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(-1.2984315) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8093811) q[1];
sx q[1];
rz(-2.2978133) q[1];
sx q[1];
rz(-1.0512645) q[1];
x q[2];
rz(-1.689765) q[3];
sx q[3];
rz(-0.61421466) q[3];
sx q[3];
rz(3.0005232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(-2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0892031) q[0];
sx q[0];
rz(-2.4501778) q[0];
sx q[0];
rz(-1.0538488) q[0];
rz(-pi) q[1];
rz(0.28378758) q[2];
sx q[2];
rz(-1.5675401) q[2];
sx q[2];
rz(-0.46797215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4334129) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(-1.3604926) q[1];
rz(-pi) q[2];
rz(-1.2849502) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(-2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2597044) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(-0.23513901) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.495372) q[0];
sx q[0];
rz(-1.6618068) q[0];
sx q[0];
rz(-1.3489086) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8638641) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(3.1090528) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1675889) q[1];
sx q[1];
rz(-2.0473192) q[1];
sx q[1];
rz(1.2567026) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6607473) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(-1.2980952) q[0];
rz(-2.4457473) q[2];
sx q[2];
rz(-2.3629284) q[2];
sx q[2];
rz(2.3455182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8400152) q[1];
sx q[1];
rz(-1.3372278) q[1];
sx q[1];
rz(2.9905) q[1];
rz(1.3917543) q[3];
sx q[3];
rz(-0.92902196) q[3];
sx q[3];
rz(-0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(-0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(-0.98714978) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894831) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(0.91870086) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5802025) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(-0.26088342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0070237006) q[1];
sx q[1];
rz(-0.94788523) q[1];
sx q[1];
rz(-1.9593777) q[1];
x q[2];
rz(0.042111245) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(0.3596572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(0.063407272) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148221) q[0];
sx q[0];
rz(-1.5855256) q[0];
sx q[0];
rz(-1.2889839) q[0];
rz(-1.2811786) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(-2.9786125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8705604) q[1];
sx q[1];
rz(-2.6189657) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-0.91026129) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(1.4753301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(0.77077579) q[2];
rz(-1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-0.055667002) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768893) q[0];
sx q[0];
rz(-1.9386374) q[0];
sx q[0];
rz(2.6951615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(1.0559527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6047302) q[1];
sx q[1];
rz(-0.68568789) q[1];
sx q[1];
rz(-1.0889978) q[1];
x q[2];
rz(1.9332063) q[3];
sx q[3];
rz(-2.9299195) q[3];
sx q[3];
rz(3.1162804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253117) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(-2.4628941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020254) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(2.461117) q[0];
rz(0.46531123) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(2.0899783) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7486836) q[1];
sx q[1];
rz(-0.52782413) q[1];
sx q[1];
rz(1.4455568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.802556) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(-0.26089222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-0.52694595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(-0.35993872) q[0];
x q[1];
rz(-0.2692659) q[2];
sx q[2];
rz(-0.93062799) q[2];
sx q[2];
rz(0.38538853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2745797) q[1];
sx q[1];
rz(-1.305294) q[1];
sx q[1];
rz(-0.32151476) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(0.3607761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11689582) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(0.025339729) q[0];
rz(-pi) q[1];
rz(-2.0759517) q[2];
sx q[2];
rz(-0.2444707) q[2];
sx q[2];
rz(2.5214362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4611778) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(-1.3780891) q[1];
rz(-pi) q[2];
rz(2.4264614) q[3];
sx q[3];
rz(-2.9188041) q[3];
sx q[3];
rz(0.95105329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0439904) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(-1.7815331) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.7760989) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(-2.3902262) q[3];
sx q[3];
rz(-1.02117) q[3];
sx q[3];
rz(1.1985967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
