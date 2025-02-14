OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(-2.9232803) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(2.6638439) q[1];
sx q[1];
rz(12.01241) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57620063) q[0];
sx q[0];
rz(-1.0668653) q[0];
sx q[0];
rz(2.6070057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5114097) q[2];
sx q[2];
rz(-2.0989387) q[2];
sx q[2];
rz(0.33231653) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3240149) q[1];
sx q[1];
rz(-2.2715873) q[1];
sx q[1];
rz(0.27591095) q[1];
x q[2];
rz(1.9945108) q[3];
sx q[3];
rz(-1.0726352) q[3];
sx q[3];
rz(0.26396449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37609425) q[2];
sx q[2];
rz(-2.8480397) q[2];
sx q[2];
rz(-1.2676839) q[2];
rz(2.3119161) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(-0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0881398) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.7470737) q[0];
rz(-2.246619) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5825533) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5916067) q[0];
sx q[0];
rz(-0.9261407) q[0];
sx q[0];
rz(0.87199535) q[0];
rz(-1.9910013) q[2];
sx q[2];
rz(-1.8981427) q[2];
sx q[2];
rz(-0.98132747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5176508) q[1];
sx q[1];
rz(-1.2576767) q[1];
sx q[1];
rz(0.70065686) q[1];
x q[2];
rz(-2.8361735) q[3];
sx q[3];
rz(-1.9807743) q[3];
sx q[3];
rz(2.3596606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6643657) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(0.030755432) q[2];
rz(2.6886046) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(2.9636813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396616) q[0];
sx q[0];
rz(-0.10098305) q[0];
sx q[0];
rz(-2.3133551) q[0];
rz(-0.047686934) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(1.9140859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2313034) q[0];
sx q[0];
rz(-2.4132015) q[0];
sx q[0];
rz(0.69407082) q[0];
rz(-pi) q[1];
rz(-1.9768049) q[2];
sx q[2];
rz(-1.3719402) q[2];
sx q[2];
rz(0.59954294) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3674254) q[1];
sx q[1];
rz(-1.8313421) q[1];
sx q[1];
rz(-2.0266286) q[1];
rz(-0.84215409) q[3];
sx q[3];
rz(-1.6098566) q[3];
sx q[3];
rz(1.0674764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.092827169) q[2];
sx q[2];
rz(-0.91492492) q[2];
sx q[2];
rz(-2.0406593) q[2];
rz(-0.01384211) q[3];
sx q[3];
rz(-1.775454) q[3];
sx q[3];
rz(0.83465105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58532995) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(-0.334326) q[0];
rz(0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(-0.11071959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3363655) q[0];
sx q[0];
rz(-2.3799161) q[0];
sx q[0];
rz(0.021999981) q[0];
rz(-pi) q[1];
rz(2.6312073) q[2];
sx q[2];
rz(-2.9523627) q[2];
sx q[2];
rz(0.58923474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19921215) q[1];
sx q[1];
rz(-0.95006493) q[1];
sx q[1];
rz(2.278028) q[1];
rz(0.7752876) q[3];
sx q[3];
rz(-1.2167769) q[3];
sx q[3];
rz(0.86590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4577786) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(-1.4078183) q[2];
rz(1.9862566) q[3];
sx q[3];
rz(-1.1563533) q[3];
sx q[3];
rz(-2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487576) q[0];
sx q[0];
rz(-0.91788569) q[0];
sx q[0];
rz(-0.99739972) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(1.4195199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62904692) q[0];
sx q[0];
rz(-1.7897494) q[0];
sx q[0];
rz(-0.16694582) q[0];
rz(-pi) q[1];
rz(2.3676374) q[2];
sx q[2];
rz(-1.4593235) q[2];
sx q[2];
rz(-0.74711266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0599654) q[1];
sx q[1];
rz(-0.94725376) q[1];
sx q[1];
rz(0.39638955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2577293) q[3];
sx q[3];
rz(-0.12774865) q[3];
sx q[3];
rz(1.7614438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2898966) q[2];
sx q[2];
rz(-3.0471314) q[2];
sx q[2];
rz(-2.8186901) q[2];
rz(1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-2.7285301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.7776529) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(-0.80663484) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-2.0239425) q[1];
sx q[1];
rz(1.1899828) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34790137) q[0];
sx q[0];
rz(-2.7712203) q[0];
sx q[0];
rz(-1.9166975) q[0];
rz(-pi) q[1];
rz(-0.34752589) q[2];
sx q[2];
rz(-0.56171562) q[2];
sx q[2];
rz(2.7507741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2014043) q[1];
sx q[1];
rz(-1.40213) q[1];
sx q[1];
rz(-3.0894214) q[1];
x q[2];
rz(1.4011995) q[3];
sx q[3];
rz(-0.84753321) q[3];
sx q[3];
rz(-2.4270428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7334062) q[2];
sx q[2];
rz(-1.6383645) q[2];
sx q[2];
rz(-0.1850941) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(-0.68814284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9516893) q[0];
sx q[0];
rz(-2.1898495) q[0];
sx q[0];
rz(0.21251799) q[0];
rz(0.7849794) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(0.77883887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76701421) q[0];
sx q[0];
rz(-0.85142577) q[0];
sx q[0];
rz(-0.79921754) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1981702) q[2];
sx q[2];
rz(-1.0154795) q[2];
sx q[2];
rz(0.057387847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0743588) q[1];
sx q[1];
rz(-1.7869084) q[1];
sx q[1];
rz(-1.8317779) q[1];
rz(-2.6565348) q[3];
sx q[3];
rz(-0.9204671) q[3];
sx q[3];
rz(-1.2984087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51398858) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(0.40840515) q[2];
rz(-0.70185316) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.7931165) q[0];
sx q[0];
rz(-1.9261253) q[0];
sx q[0];
rz(0.066019639) q[0];
rz(1.5090212) q[1];
sx q[1];
rz(-1.0450109) q[1];
sx q[1];
rz(0.95796934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1692266) q[0];
sx q[0];
rz(-2.027085) q[0];
sx q[0];
rz(1.8877939) q[0];
rz(0.81664576) q[2];
sx q[2];
rz(-1.9672251) q[2];
sx q[2];
rz(2.3289837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66855318) q[1];
sx q[1];
rz(-0.91382342) q[1];
sx q[1];
rz(1.9690352) q[1];
rz(-pi) q[2];
rz(2.1470966) q[3];
sx q[3];
rz(-0.18261431) q[3];
sx q[3];
rz(-1.8908569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8336739) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(-1.4109122) q[2];
rz(0.69502568) q[3];
sx q[3];
rz(-1.6618988) q[3];
sx q[3];
rz(-3.1237349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0787635) q[0];
sx q[0];
rz(-1.6595027) q[0];
sx q[0];
rz(-0.98989809) q[0];
rz(2.6784189) q[1];
sx q[1];
rz(-1.6894692) q[1];
sx q[1];
rz(-1.90082) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.038886) q[0];
sx q[0];
rz(-2.8409344) q[0];
sx q[0];
rz(1.7666398) q[0];
x q[1];
rz(-1.1616917) q[2];
sx q[2];
rz(-0.4767524) q[2];
sx q[2];
rz(-2.7331405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5338143) q[1];
sx q[1];
rz(-1.8318212) q[1];
sx q[1];
rz(0.29284524) q[1];
rz(-pi) q[2];
rz(2.1911591) q[3];
sx q[3];
rz(-1.9470805) q[3];
sx q[3];
rz(0.84445124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3339633) q[2];
sx q[2];
rz(-1.3796076) q[2];
sx q[2];
rz(-2.4228952) q[2];
rz(-1.9888318) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54866791) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(0.80192178) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(-0.71802872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336003) q[0];
sx q[0];
rz(-1.2498858) q[0];
sx q[0];
rz(-0.72707392) q[0];
x q[1];
rz(2.3921591) q[2];
sx q[2];
rz(-2.4735056) q[2];
sx q[2];
rz(1.4632478) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1433318) q[1];
sx q[1];
rz(-1.623053) q[1];
sx q[1];
rz(1.9027756) q[1];
rz(-0.84661412) q[3];
sx q[3];
rz(-2.1606956) q[3];
sx q[3];
rz(0.55800948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7903018) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(2.6289319) q[2];
rz(-0.74448186) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(0.7640394) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9164593) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(2.4254639) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(0.11083435) q[2];
sx q[2];
rz(-1.7946984) q[2];
sx q[2];
rz(2.6738965) q[2];
rz(-1.2342831) q[3];
sx q[3];
rz(-1.6675622) q[3];
sx q[3];
rz(0.29940816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
