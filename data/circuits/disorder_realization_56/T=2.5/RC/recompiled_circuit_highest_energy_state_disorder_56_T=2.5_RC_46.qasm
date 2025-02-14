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
rz(2.9958148) q[0];
sx q[0];
rz(-1.79359) q[0];
sx q[0];
rz(-0.3682799) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(-2.2987125) q[1];
sx q[1];
rz(0.24922961) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6347353) q[0];
sx q[0];
rz(-1.0930499) q[0];
sx q[0];
rz(1.696582) q[0];
rz(-2.1051718) q[2];
sx q[2];
rz(-1.8412565) q[2];
sx q[2];
rz(1.1391183) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8650444) q[1];
sx q[1];
rz(-2.2575592) q[1];
sx q[1];
rz(-1.264132) q[1];
rz(-pi) q[2];
rz(0.35855583) q[3];
sx q[3];
rz(-1.6772146) q[3];
sx q[3];
rz(0.88332435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2117846) q[2];
sx q[2];
rz(-0.70694184) q[2];
sx q[2];
rz(0.24016538) q[2];
rz(-2.5943878) q[3];
sx q[3];
rz(-1.5144843) q[3];
sx q[3];
rz(2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7996247) q[0];
sx q[0];
rz(-0.37536055) q[0];
sx q[0];
rz(2.4976835) q[0];
rz(-1.0470692) q[1];
sx q[1];
rz(-1.964485) q[1];
sx q[1];
rz(0.26328304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851259) q[0];
sx q[0];
rz(-1.5432285) q[0];
sx q[0];
rz(-3.1113202) q[0];
rz(-1.5059169) q[2];
sx q[2];
rz(-1.5047964) q[2];
sx q[2];
rz(-2.9213816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4810679) q[1];
sx q[1];
rz(-2.1612289) q[1];
sx q[1];
rz(0.51312889) q[1];
x q[2];
rz(-2.1016779) q[3];
sx q[3];
rz(-0.71456075) q[3];
sx q[3];
rz(-1.2604063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(2.4158884) q[2];
rz(2.8343685) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(1.0544302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211408) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(-2.3725574) q[0];
rz(-0.10737315) q[1];
sx q[1];
rz(-1.702405) q[1];
sx q[1];
rz(-0.69951397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.526992) q[0];
sx q[0];
rz(-1.8934947) q[0];
sx q[0];
rz(-1.931577) q[0];
rz(-1.5753463) q[2];
sx q[2];
rz(-2.1855436) q[2];
sx q[2];
rz(-0.67753032) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2638638) q[1];
sx q[1];
rz(-2.6843964) q[1];
sx q[1];
rz(1.3844107) q[1];
x q[2];
rz(2.3885257) q[3];
sx q[3];
rz(-2.1478466) q[3];
sx q[3];
rz(0.53423877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37381441) q[2];
sx q[2];
rz(-0.57743293) q[2];
sx q[2];
rz(-0.90816298) q[2];
rz(2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813756) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(-2.6326219) q[0];
rz(3.0834037) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(-2.0740654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876388) q[0];
sx q[0];
rz(-1.2042856) q[0];
sx q[0];
rz(-2.3838759) q[0];
rz(1.7420962) q[2];
sx q[2];
rz(-0.41738415) q[2];
sx q[2];
rz(-2.4462552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4719676) q[1];
sx q[1];
rz(-1.4635007) q[1];
sx q[1];
rz(1.4221684) q[1];
rz(-1.8640609) q[3];
sx q[3];
rz(-1.6365657) q[3];
sx q[3];
rz(-1.6848967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1099781) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(0.99008375) q[2];
rz(1.3958942) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(-1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083967) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(2.0268832) q[1];
sx q[1];
rz(-1.474294) q[1];
sx q[1];
rz(1.4310736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1170332) q[0];
sx q[0];
rz(-0.22732553) q[0];
sx q[0];
rz(2.2973799) q[0];
rz(-pi) q[1];
rz(-2.0189907) q[2];
sx q[2];
rz(-2.5914874) q[2];
sx q[2];
rz(-3.0193626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.614735) q[1];
sx q[1];
rz(-2.4522374) q[1];
sx q[1];
rz(-0.3853285) q[1];
rz(-pi) q[2];
rz(1.4214755) q[3];
sx q[3];
rz(-1.3589371) q[3];
sx q[3];
rz(-0.52313357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7913671) q[2];
sx q[2];
rz(-0.75608772) q[2];
sx q[2];
rz(1.4678601) q[2];
rz(1.0427467) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(0.79152542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4333711) q[0];
sx q[0];
rz(-1.8480166) q[0];
sx q[0];
rz(0.21155393) q[0];
rz(-0.64282203) q[1];
sx q[1];
rz(-2.0170409) q[1];
sx q[1];
rz(-2.0932253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12473561) q[0];
sx q[0];
rz(-2.2679459) q[0];
sx q[0];
rz(-1.7378896) q[0];
x q[1];
rz(2.8564212) q[2];
sx q[2];
rz(-2.9551571) q[2];
sx q[2];
rz(0.16229072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.590946) q[1];
sx q[1];
rz(-0.63884402) q[1];
sx q[1];
rz(-0.039012564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5023175) q[3];
sx q[3];
rz(-2.5279495) q[3];
sx q[3];
rz(-0.63921038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6362777) q[2];
sx q[2];
rz(-1.2039528) q[2];
sx q[2];
rz(2.8875276) q[2];
rz(0.17357477) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(1.180163) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59488615) q[0];
sx q[0];
rz(-0.43989023) q[0];
sx q[0];
rz(-0.79676262) q[0];
rz(-2.2548389) q[1];
sx q[1];
rz(-1.3462857) q[1];
sx q[1];
rz(-2.5409882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89413801) q[0];
sx q[0];
rz(-2.8464212) q[0];
sx q[0];
rz(2.9516793) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7037292) q[2];
sx q[2];
rz(-2.3355451) q[2];
sx q[2];
rz(1.9423167) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.826304) q[1];
sx q[1];
rz(-0.70543324) q[1];
sx q[1];
rz(-1.3194196) q[1];
rz(-2.10403) q[3];
sx q[3];
rz(-1.8775465) q[3];
sx q[3];
rz(0.83707929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.1116921) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(2.0491484) q[2];
rz(1.7732636) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261037) q[0];
sx q[0];
rz(-1.1142718) q[0];
sx q[0];
rz(2.6901167) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(2.480004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115727) q[0];
sx q[0];
rz(-1.9826188) q[0];
sx q[0];
rz(-2.9020082) q[0];
x q[1];
rz(1.6986474) q[2];
sx q[2];
rz(-0.50179728) q[2];
sx q[2];
rz(1.0867829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3319131) q[1];
sx q[1];
rz(-0.47985489) q[1];
sx q[1];
rz(-1.859568) q[1];
rz(0.39644661) q[3];
sx q[3];
rz(-1.3260734) q[3];
sx q[3];
rz(2.1978767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0852069) q[2];
sx q[2];
rz(-1.8737917) q[2];
sx q[2];
rz(0.80643225) q[2];
rz(1.8993529) q[3];
sx q[3];
rz(-2.9759088) q[3];
sx q[3];
rz(3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-0.30271444) q[0];
sx q[0];
rz(-2.0063722) q[0];
sx q[0];
rz(0.43933991) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-0.83963436) q[1];
sx q[1];
rz(2.1913948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882975) q[0];
sx q[0];
rz(-1.9425009) q[0];
sx q[0];
rz(-0.88754334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15301159) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(0.52545122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62225759) q[1];
sx q[1];
rz(-0.37168113) q[1];
sx q[1];
rz(1.8628013) q[1];
rz(2.1860649) q[3];
sx q[3];
rz(-0.63378171) q[3];
sx q[3];
rz(-1.2891226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83303678) q[2];
sx q[2];
rz(-0.50806442) q[2];
sx q[2];
rz(1.9527831) q[2];
rz(-3.0360119) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(-1.0134491) q[3];
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
rz(-0.078876615) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(-2.4421332) q[0];
rz(-1.7309011) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(-2.9404822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5546306) q[0];
sx q[0];
rz(-2.2788958) q[0];
sx q[0];
rz(0.59815852) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2314448) q[2];
sx q[2];
rz(-2.4322699) q[2];
sx q[2];
rz(-0.34914474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0929151) q[1];
sx q[1];
rz(-0.72775048) q[1];
sx q[1];
rz(-1.2341094) q[1];
rz(-pi) q[2];
rz(0.46259677) q[3];
sx q[3];
rz(-2.9411081) q[3];
sx q[3];
rz(-1.5104664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.19645277) q[2];
sx q[2];
rz(-1.2597193) q[2];
sx q[2];
rz(0.067616612) q[2];
rz(-1.128528) q[3];
sx q[3];
rz(-2.0566548) q[3];
sx q[3];
rz(1.5170001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45831281) q[0];
sx q[0];
rz(-1.0673609) q[0];
sx q[0];
rz(1.3491032) q[0];
rz(1.6364527) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(-1.3624874) q[2];
sx q[2];
rz(-0.97803331) q[2];
sx q[2];
rz(-2.5292653) q[2];
rz(1.296386) q[3];
sx q[3];
rz(-1.8897703) q[3];
sx q[3];
rz(-0.29714656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
