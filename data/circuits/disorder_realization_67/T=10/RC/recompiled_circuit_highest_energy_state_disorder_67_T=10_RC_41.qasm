OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.500904) q[0];
sx q[0];
rz(-0.7839497) q[0];
sx q[0];
rz(3.0408903) q[0];
rz(-4.5667629) q[1];
sx q[1];
rz(6.9520091) q[1];
sx q[1];
rz(7.8520757) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9933321) q[0];
sx q[0];
rz(-1.0251704) q[0];
sx q[0];
rz(1.7094897) q[0];
x q[1];
rz(-2.8694081) q[2];
sx q[2];
rz(-1.4708752) q[2];
sx q[2];
rz(1.3306432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4922759) q[1];
sx q[1];
rz(-1.7360949) q[1];
sx q[1];
rz(0.16053546) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4675191) q[3];
sx q[3];
rz(-2.010502) q[3];
sx q[3];
rz(1.7252418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
sx q[2];
rz(-1.5300306) q[2];
rz(-1.642646) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(-2.4422395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4636369) q[0];
sx q[0];
rz(-1.7706484) q[0];
sx q[0];
rz(-0.26588765) q[0];
rz(-1.5222585) q[1];
sx q[1];
rz(-1.8306754) q[1];
sx q[1];
rz(-1.486385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7911969) q[0];
sx q[0];
rz(-0.7753709) q[0];
sx q[0];
rz(0.010678963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56358524) q[2];
sx q[2];
rz(-1.775309) q[2];
sx q[2];
rz(-2.9737008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7922272) q[1];
sx q[1];
rz(-1.1046032) q[1];
sx q[1];
rz(-1.8833877) q[1];
x q[2];
rz(-0.19274917) q[3];
sx q[3];
rz(-1.1067353) q[3];
sx q[3];
rz(2.131586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3050401) q[2];
sx q[2];
rz(-1.0779447) q[2];
sx q[2];
rz(0.067122785) q[2];
rz(0.51131311) q[3];
sx q[3];
rz(-1.7917683) q[3];
sx q[3];
rz(-0.93222031) q[3];
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
rz(2.9028645) q[0];
sx q[0];
rz(-1.4816062) q[0];
sx q[0];
rz(2.9929152) q[0];
rz(2.6909434) q[1];
sx q[1];
rz(-1.5573749) q[1];
sx q[1];
rz(0.91948909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7203569) q[0];
sx q[0];
rz(-1.8381869) q[0];
sx q[0];
rz(-1.8084099) q[0];
rz(1.3344263) q[2];
sx q[2];
rz(-0.70118517) q[2];
sx q[2];
rz(-1.2126306) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0362404) q[1];
sx q[1];
rz(-2.4215464) q[1];
sx q[1];
rz(2.4612263) q[1];
x q[2];
rz(2.6387003) q[3];
sx q[3];
rz(-0.38030312) q[3];
sx q[3];
rz(-1.2740434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.273169) q[2];
sx q[2];
rz(-1.0852852) q[2];
sx q[2];
rz(-2.2853509) q[2];
rz(2.26561) q[3];
sx q[3];
rz(-2.7464505) q[3];
sx q[3];
rz(-1.1538848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3530537) q[0];
sx q[0];
rz(-1.5748698) q[0];
sx q[0];
rz(-2.5153611) q[0];
rz(-1.7010472) q[1];
sx q[1];
rz(-0.82653058) q[1];
sx q[1];
rz(-0.70762077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514446) q[0];
sx q[0];
rz(-2.2934545) q[0];
sx q[0];
rz(0.11369205) q[0];
rz(-0.78379102) q[2];
sx q[2];
rz(-1.3221127) q[2];
sx q[2];
rz(0.75088203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7654082) q[1];
sx q[1];
rz(-1.4896684) q[1];
sx q[1];
rz(2.8721832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.510101) q[3];
sx q[3];
rz(-0.48572054) q[3];
sx q[3];
rz(-0.39050366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6669199) q[2];
sx q[2];
rz(-1.2771353) q[2];
sx q[2];
rz(0.39371583) q[2];
rz(-2.3422824) q[3];
sx q[3];
rz(-2.7768713) q[3];
sx q[3];
rz(-1.5289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9202775) q[0];
sx q[0];
rz(-2.3096313) q[0];
sx q[0];
rz(-3.063391) q[0];
rz(-0.68866628) q[1];
sx q[1];
rz(-2.2024901) q[1];
sx q[1];
rz(0.92855531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9445223) q[0];
sx q[0];
rz(-1.8816948) q[0];
sx q[0];
rz(0.59493712) q[0];
rz(-2.2735394) q[2];
sx q[2];
rz(-2.6871024) q[2];
sx q[2];
rz(-2.8502685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9282058) q[1];
sx q[1];
rz(-1.3524071) q[1];
sx q[1];
rz(2.7534927) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6070091) q[3];
sx q[3];
rz(-2.4494736) q[3];
sx q[3];
rz(0.87833709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3277305) q[2];
sx q[2];
rz(-2.6321415) q[2];
sx q[2];
rz(0.65208411) q[2];
rz(-0.70513519) q[3];
sx q[3];
rz(-1.4950246) q[3];
sx q[3];
rz(-2.7130073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9014277) q[0];
sx q[0];
rz(-0.21900284) q[0];
sx q[0];
rz(1.3448311) q[0];
rz(0.52974686) q[1];
sx q[1];
rz(-1.8579282) q[1];
sx q[1];
rz(1.8240428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12152762) q[0];
sx q[0];
rz(-2.1449267) q[0];
sx q[0];
rz(-3.0453389) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7129301) q[2];
sx q[2];
rz(-0.33944079) q[2];
sx q[2];
rz(1.0975099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3287813) q[1];
sx q[1];
rz(-1.9095632) q[1];
sx q[1];
rz(-3.0121606) q[1];
rz(-pi) q[2];
rz(-2.0798912) q[3];
sx q[3];
rz(-1.7359043) q[3];
sx q[3];
rz(-2.3632339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.162447) q[2];
sx q[2];
rz(-1.2787168) q[2];
sx q[2];
rz(-2.0549959) q[2];
rz(0.094001683) q[3];
sx q[3];
rz(-2.0407929) q[3];
sx q[3];
rz(0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(3.060044) q[0];
rz(-1.6185919) q[1];
sx q[1];
rz(-2.1269507) q[1];
sx q[1];
rz(0.60752216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.994068) q[0];
sx q[0];
rz(-2.3660064) q[0];
sx q[0];
rz(-1.8344384) q[0];
rz(-1.188368) q[2];
sx q[2];
rz(-2.0086346) q[2];
sx q[2];
rz(3.0944648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0237777) q[1];
sx q[1];
rz(-1.5842321) q[1];
sx q[1];
rz(-1.4680844) q[1];
rz(-pi) q[2];
rz(2.7681302) q[3];
sx q[3];
rz(-1.6420396) q[3];
sx q[3];
rz(0.072865818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0532694) q[2];
sx q[2];
rz(-2.93556) q[2];
sx q[2];
rz(-0.44331178) q[2];
rz(-0.045470227) q[3];
sx q[3];
rz(-1.1683522) q[3];
sx q[3];
rz(-0.62062353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745673) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(-2.5734651) q[0];
rz(1.8727632) q[1];
sx q[1];
rz(-1.9867691) q[1];
sx q[1];
rz(0.61187977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9107008) q[0];
sx q[0];
rz(-2.7534427) q[0];
sx q[0];
rz(2.9064473) q[0];
rz(-pi) q[1];
x q[1];
rz(2.955825) q[2];
sx q[2];
rz(-1.298873) q[2];
sx q[2];
rz(2.853554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.001732262) q[1];
sx q[1];
rz(-2.0563158) q[1];
sx q[1];
rz(2.8667169) q[1];
rz(-2.2886226) q[3];
sx q[3];
rz(-0.26085869) q[3];
sx q[3];
rz(-2.0776231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.890471) q[2];
sx q[2];
rz(-1.2476363) q[2];
sx q[2];
rz(-2.9337511) q[2];
rz(-0.017092997) q[3];
sx q[3];
rz(-2.1199675) q[3];
sx q[3];
rz(-2.0489073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51181483) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(1.7275607) q[0];
rz(0.066970197) q[1];
sx q[1];
rz(-1.8194865) q[1];
sx q[1];
rz(-3.0019143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0020376677) q[0];
sx q[0];
rz(-1.25601) q[0];
sx q[0];
rz(-2.806753) q[0];
rz(-pi) q[1];
rz(2.5508444) q[2];
sx q[2];
rz(-2.4409238) q[2];
sx q[2];
rz(-0.76393647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5433257) q[1];
sx q[1];
rz(-1.7964592) q[1];
sx q[1];
rz(1.3545827) q[1];
rz(-pi) q[2];
rz(-1.2788676) q[3];
sx q[3];
rz(-0.17660429) q[3];
sx q[3];
rz(0.3738598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.70860538) q[2];
sx q[2];
rz(-1.8726417) q[2];
sx q[2];
rz(2.404786) q[2];
rz(-1.6144729) q[3];
sx q[3];
rz(-1.0090642) q[3];
sx q[3];
rz(1.2705151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0438185) q[0];
sx q[0];
rz(-0.96776861) q[0];
sx q[0];
rz(-1.9737825) q[0];
rz(-2.4763926) q[1];
sx q[1];
rz(-1.2780814) q[1];
sx q[1];
rz(2.1591689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2806399) q[0];
sx q[0];
rz(-2.6053944) q[0];
sx q[0];
rz(1.9704738) q[0];
rz(1.8389687) q[2];
sx q[2];
rz(-0.94899584) q[2];
sx q[2];
rz(1.5121429) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9537045) q[1];
sx q[1];
rz(-2.4361103) q[1];
sx q[1];
rz(3.0445552) q[1];
rz(-pi) q[2];
rz(1.1725964) q[3];
sx q[3];
rz(-1.5086344) q[3];
sx q[3];
rz(-2.8650995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3490225) q[2];
sx q[2];
rz(-2.4945365) q[2];
sx q[2];
rz(1.6046074) q[2];
rz(-2.1215306) q[3];
sx q[3];
rz(-1.5217109) q[3];
sx q[3];
rz(2.9807239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0237324) q[0];
sx q[0];
rz(-1.5845789) q[0];
sx q[0];
rz(-1.3316863) q[0];
rz(2.3125519) q[1];
sx q[1];
rz(-1.6571028) q[1];
sx q[1];
rz(2.170457) q[1];
rz(2.7132158) q[2];
sx q[2];
rz(-0.72441341) q[2];
sx q[2];
rz(-2.7073467) q[2];
rz(0.42468023) q[3];
sx q[3];
rz(-1.1241354) q[3];
sx q[3];
rz(-0.507171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
