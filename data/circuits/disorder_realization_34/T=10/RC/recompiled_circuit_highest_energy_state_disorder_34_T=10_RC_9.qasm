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
rz(1.0518987) q[0];
sx q[0];
rz(2.9419152) q[0];
sx q[0];
rz(10.022104) q[0];
rz(-0.53961331) q[1];
sx q[1];
rz(-0.55066723) q[1];
sx q[1];
rz(0.82672969) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4524627) q[0];
sx q[0];
rz(-2.0553382) q[0];
sx q[0];
rz(-1.9758805) q[0];
rz(-0.8491926) q[2];
sx q[2];
rz(-1.7623644) q[2];
sx q[2];
rz(1.220773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.93365) q[1];
sx q[1];
rz(-2.4550555) q[1];
sx q[1];
rz(2.0523493) q[1];
rz(-pi) q[2];
rz(0.22181819) q[3];
sx q[3];
rz(-2.6008587) q[3];
sx q[3];
rz(2.3774035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4248767) q[2];
sx q[2];
rz(-2.0108607) q[2];
sx q[2];
rz(2.7190599) q[2];
rz(-1.9033963) q[3];
sx q[3];
rz(-2.7178552) q[3];
sx q[3];
rz(-2.197926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0331405) q[0];
sx q[0];
rz(-0.55522454) q[0];
sx q[0];
rz(2.2174477) q[0];
rz(2.4657501) q[1];
sx q[1];
rz(-1.5166413) q[1];
sx q[1];
rz(-2.1088375) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076457523) q[0];
sx q[0];
rz(-1.9357114) q[0];
sx q[0];
rz(2.5315842) q[0];
rz(-0.087914596) q[2];
sx q[2];
rz(-2.9682142) q[2];
sx q[2];
rz(-1.8176469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38235462) q[1];
sx q[1];
rz(-1.865631) q[1];
sx q[1];
rz(1.0947541) q[1];
rz(-pi) q[2];
rz(2.1274125) q[3];
sx q[3];
rz(-0.92837205) q[3];
sx q[3];
rz(-3.0730221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5702901) q[2];
sx q[2];
rz(-1.4455659) q[2];
sx q[2];
rz(-0.22682133) q[2];
rz(-1.4443385) q[3];
sx q[3];
rz(-3.0310013) q[3];
sx q[3];
rz(0.9816696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.084376) q[0];
sx q[0];
rz(-2.9380874) q[0];
sx q[0];
rz(0.89619005) q[0];
rz(2.8016688) q[1];
sx q[1];
rz(-1.8819921) q[1];
sx q[1];
rz(2.6673754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8086169) q[0];
sx q[0];
rz(-1.364437) q[0];
sx q[0];
rz(-1.1550857) q[0];
x q[1];
rz(-1.3856085) q[2];
sx q[2];
rz(-1.5366925) q[2];
sx q[2];
rz(-2.2540384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5061298) q[1];
sx q[1];
rz(-0.78054192) q[1];
sx q[1];
rz(0.64029633) q[1];
rz(2.3044448) q[3];
sx q[3];
rz(-1.7230534) q[3];
sx q[3];
rz(2.4796132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.077896) q[2];
sx q[2];
rz(-1.6113969) q[2];
sx q[2];
rz(-1.873675) q[2];
rz(2.7481713) q[3];
sx q[3];
rz(-2.7712517) q[3];
sx q[3];
rz(-0.89746499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1657408) q[0];
sx q[0];
rz(-0.535088) q[0];
sx q[0];
rz(-2.72056) q[0];
rz(-0.16972217) q[1];
sx q[1];
rz(-1.38849) q[1];
sx q[1];
rz(1.1441182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4108858) q[0];
sx q[0];
rz(-0.050844103) q[0];
sx q[0];
rz(0.4836785) q[0];
rz(-1.3058) q[2];
sx q[2];
rz(-0.41171968) q[2];
sx q[2];
rz(0.078814332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1269234) q[1];
sx q[1];
rz(-1.3433392) q[1];
sx q[1];
rz(2.3421351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3973976) q[3];
sx q[3];
rz(-2.5928516) q[3];
sx q[3];
rz(-2.5732267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0890395) q[2];
sx q[2];
rz(-1.7567987) q[2];
sx q[2];
rz(2.4844737) q[2];
rz(-1.1013912) q[3];
sx q[3];
rz(-1.907932) q[3];
sx q[3];
rz(-1.5677946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9927486) q[0];
sx q[0];
rz(-1.0555457) q[0];
sx q[0];
rz(1.2626935) q[0];
rz(2.8008723) q[1];
sx q[1];
rz(-1.6993399) q[1];
sx q[1];
rz(2.5943894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1520098) q[0];
sx q[0];
rz(-0.18704601) q[0];
sx q[0];
rz(1.90045) q[0];
rz(-1.5261602) q[2];
sx q[2];
rz(-1.6812162) q[2];
sx q[2];
rz(-2.8912889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62740993) q[1];
sx q[1];
rz(-1.0189432) q[1];
sx q[1];
rz(2.7937006) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30468632) q[3];
sx q[3];
rz(-1.2126367) q[3];
sx q[3];
rz(0.55613925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1269232) q[2];
sx q[2];
rz(-1.1875443) q[2];
sx q[2];
rz(-2.8200601) q[2];
rz(2.1899636) q[3];
sx q[3];
rz(-1.8961597) q[3];
sx q[3];
rz(1.2661701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6283918) q[0];
sx q[0];
rz(-0.43190792) q[0];
sx q[0];
rz(0.54650724) q[0];
rz(2.6642117) q[1];
sx q[1];
rz(-0.42196208) q[1];
sx q[1];
rz(-1.2786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0437677) q[0];
sx q[0];
rz(-1.861933) q[0];
sx q[0];
rz(1.2173675) q[0];
rz(3.0307604) q[2];
sx q[2];
rz(-2.5907907) q[2];
sx q[2];
rz(0.11167311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1977935) q[1];
sx q[1];
rz(-1.1778801) q[1];
sx q[1];
rz(2.7019634) q[1];
rz(2.8987721) q[3];
sx q[3];
rz(-0.87174635) q[3];
sx q[3];
rz(0.2350014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4930341) q[2];
sx q[2];
rz(-2.9849122) q[2];
sx q[2];
rz(2.1873761) q[2];
rz(2.3844125) q[3];
sx q[3];
rz(-1.4933519) q[3];
sx q[3];
rz(-0.89603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2467982) q[0];
sx q[0];
rz(-3.09258) q[0];
sx q[0];
rz(-3.0931296) q[0];
rz(0.63352829) q[1];
sx q[1];
rz(-2.3506479) q[1];
sx q[1];
rz(-1.6652426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8970743) q[0];
sx q[0];
rz(-0.71673043) q[0];
sx q[0];
rz(1.1900363) q[0];
rz(-pi) q[1];
rz(-2.9505312) q[2];
sx q[2];
rz(-1.3891561) q[2];
sx q[2];
rz(-1.5921648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39346972) q[1];
sx q[1];
rz(-0.49897614) q[1];
sx q[1];
rz(-0.62403719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7853297) q[3];
sx q[3];
rz(-1.8275785) q[3];
sx q[3];
rz(0.40865001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1518636) q[2];
sx q[2];
rz(-2.180474) q[2];
sx q[2];
rz(0.1667008) q[2];
rz(1.9549595) q[3];
sx q[3];
rz(-1.229076) q[3];
sx q[3];
rz(2.1738539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1879021) q[0];
sx q[0];
rz(-0.1583651) q[0];
sx q[0];
rz(0.70575356) q[0];
rz(-1.3450274) q[1];
sx q[1];
rz(-1.8555565) q[1];
sx q[1];
rz(0.3261303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46485422) q[0];
sx q[0];
rz(-0.87552545) q[0];
sx q[0];
rz(-2.763163) q[0];
rz(1.5495315) q[2];
sx q[2];
rz(-2.6274963) q[2];
sx q[2];
rz(-1.7618881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5183954) q[1];
sx q[1];
rz(-2.9120486) q[1];
sx q[1];
rz(-2.3238682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9097861) q[3];
sx q[3];
rz(-0.44054809) q[3];
sx q[3];
rz(-2.6772628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2976133) q[2];
sx q[2];
rz(-1.16301) q[2];
sx q[2];
rz(-2.6824717) q[2];
rz(1.8504359) q[3];
sx q[3];
rz(-1.922311) q[3];
sx q[3];
rz(0.084376924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5463663) q[0];
sx q[0];
rz(-1.3331174) q[0];
sx q[0];
rz(3.1072733) q[0];
rz(-1.8345087) q[1];
sx q[1];
rz(-2.391075) q[1];
sx q[1];
rz(-1.7274571) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.14489) q[0];
sx q[0];
rz(-1.7976947) q[0];
sx q[0];
rz(-2.4064896) q[0];
x q[1];
rz(-2.6488536) q[2];
sx q[2];
rz(-0.73046267) q[2];
sx q[2];
rz(-1.7809389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.203277) q[1];
sx q[1];
rz(-1.1585981) q[1];
sx q[1];
rz(2.5286872) q[1];
rz(-pi) q[2];
rz(1.5758547) q[3];
sx q[3];
rz(-1.3225319) q[3];
sx q[3];
rz(0.92456532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2109334) q[2];
sx q[2];
rz(-2.2948269) q[2];
sx q[2];
rz(-2.6625114) q[2];
rz(0.72862285) q[3];
sx q[3];
rz(-1.9727547) q[3];
sx q[3];
rz(2.5463879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625192) q[0];
sx q[0];
rz(-0.74051028) q[0];
sx q[0];
rz(0.24653521) q[0];
rz(2.0138228) q[1];
sx q[1];
rz(-1.133254) q[1];
sx q[1];
rz(-0.22334982) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52106954) q[0];
sx q[0];
rz(-0.4979254) q[0];
sx q[0];
rz(-3.0715406) q[0];
rz(-pi) q[1];
rz(1.3904273) q[2];
sx q[2];
rz(-1.0374238) q[2];
sx q[2];
rz(0.37796256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1829233) q[1];
sx q[1];
rz(-2.8374817) q[1];
sx q[1];
rz(1.3728549) q[1];
rz(-1.8605883) q[3];
sx q[3];
rz(-1.1613356) q[3];
sx q[3];
rz(2.0807165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.25536209) q[2];
sx q[2];
rz(-0.88408154) q[2];
sx q[2];
rz(-0.13038005) q[2];
rz(-2.5731795) q[3];
sx q[3];
rz(-1.3610898) q[3];
sx q[3];
rz(-0.4140678) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2100621) q[0];
sx q[0];
rz(-2.2097586) q[0];
sx q[0];
rz(-2.7037822) q[0];
rz(1.9991649) q[1];
sx q[1];
rz(-1.8985959) q[1];
sx q[1];
rz(1.7172145) q[1];
rz(0.78442153) q[2];
sx q[2];
rz(-0.57195819) q[2];
sx q[2];
rz(2.3637003) q[2];
rz(-2.4655399) q[3];
sx q[3];
rz(-2.65391) q[3];
sx q[3];
rz(0.20269453) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
