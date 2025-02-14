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
rz(2.6019793) q[1];
sx q[1];
rz(-2.5909254) q[1];
sx q[1];
rz(-0.82672969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4524627) q[0];
sx q[0];
rz(-1.0862545) q[0];
sx q[0];
rz(1.1657122) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8887822) q[2];
sx q[2];
rz(-0.86517715) q[2];
sx q[2];
rz(-2.9575612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3394471) q[1];
sx q[1];
rz(-0.97426322) q[1];
sx q[1];
rz(-0.36277186) q[1];
rz(-pi) q[2];
rz(-0.22181819) q[3];
sx q[3];
rz(-2.6008587) q[3];
sx q[3];
rz(-2.3774035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4248767) q[2];
sx q[2];
rz(-1.1307319) q[2];
sx q[2];
rz(0.42253271) q[2];
rz(1.9033963) q[3];
sx q[3];
rz(-0.4237375) q[3];
sx q[3];
rz(0.94366664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0331405) q[0];
sx q[0];
rz(-0.55522454) q[0];
sx q[0];
rz(0.92414498) q[0];
rz(0.67584258) q[1];
sx q[1];
rz(-1.5166413) q[1];
sx q[1];
rz(2.1088375) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1193405) q[0];
sx q[0];
rz(-2.4428833) q[0];
sx q[0];
rz(-0.58813372) q[0];
rz(-pi) q[1];
rz(-2.9688705) q[2];
sx q[2];
rz(-1.5859436) q[2];
sx q[2];
rz(-2.8081389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38235462) q[1];
sx q[1];
rz(-1.865631) q[1];
sx q[1];
rz(-2.0468385) q[1];
x q[2];
rz(2.5268484) q[3];
sx q[3];
rz(-2.3182676) q[3];
sx q[3];
rz(-2.4058482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5702901) q[2];
sx q[2];
rz(-1.6960267) q[2];
sx q[2];
rz(2.9147713) q[2];
rz(1.4443385) q[3];
sx q[3];
rz(-3.0310013) q[3];
sx q[3];
rz(-0.9816696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572166) q[0];
sx q[0];
rz(-0.20350525) q[0];
sx q[0];
rz(0.89619005) q[0];
rz(0.33992386) q[1];
sx q[1];
rz(-1.2596005) q[1];
sx q[1];
rz(2.6673754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6723876) q[0];
sx q[0];
rz(-0.46142277) q[0];
sx q[0];
rz(-1.0925596) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.034696636) q[2];
sx q[2];
rz(-1.7558753) q[2];
sx q[2];
rz(-2.4519631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82638559) q[1];
sx q[1];
rz(-0.97123324) q[1];
sx q[1];
rz(2.1050598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9378652) q[3];
sx q[3];
rz(-2.2940638) q[3];
sx q[3];
rz(2.096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0636966) q[2];
sx q[2];
rz(-1.6113969) q[2];
sx q[2];
rz(1.873675) q[2];
rz(-2.7481713) q[3];
sx q[3];
rz(-2.7712517) q[3];
sx q[3];
rz(0.89746499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1657408) q[0];
sx q[0];
rz(-0.535088) q[0];
sx q[0];
rz(-2.72056) q[0];
rz(-2.9718705) q[1];
sx q[1];
rz(-1.7531027) q[1];
sx q[1];
rz(-1.9974744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73070684) q[0];
sx q[0];
rz(-3.0907486) q[0];
sx q[0];
rz(2.6579141) q[0];
rz(1.969643) q[2];
sx q[2];
rz(-1.4657925) q[2];
sx q[2];
rz(1.4058553) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3575625) q[1];
sx q[1];
rz(-0.79748193) q[1];
sx q[1];
rz(-1.2502115) q[1];
x q[2];
rz(-2.7189488) q[3];
sx q[3];
rz(-1.9319252) q[3];
sx q[3];
rz(-0.33651957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.052553136) q[2];
sx q[2];
rz(-1.384794) q[2];
sx q[2];
rz(-2.4844737) q[2];
rz(-1.1013912) q[3];
sx q[3];
rz(-1.907932) q[3];
sx q[3];
rz(-1.5677946) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9927486) q[0];
sx q[0];
rz(-1.0555457) q[0];
sx q[0];
rz(1.2626935) q[0];
rz(0.34072033) q[1];
sx q[1];
rz(-1.4422528) q[1];
sx q[1];
rz(-0.54720324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094484821) q[0];
sx q[0];
rz(-1.63103) q[0];
sx q[0];
rz(-1.747986) q[0];
rz(1.5261602) q[2];
sx q[2];
rz(-1.4603764) q[2];
sx q[2];
rz(0.25030372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3860847) q[1];
sx q[1];
rz(-1.2762462) q[1];
sx q[1];
rz(-0.99099116) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9447127) q[3];
sx q[3];
rz(-1.2860048) q[3];
sx q[3];
rz(-2.2367331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0146694) q[2];
sx q[2];
rz(-1.1875443) q[2];
sx q[2];
rz(2.8200601) q[2];
rz(-2.1899636) q[3];
sx q[3];
rz(-1.8961597) q[3];
sx q[3];
rz(1.8754225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132009) q[0];
sx q[0];
rz(-2.7096847) q[0];
sx q[0];
rz(0.54650724) q[0];
rz(2.6642117) q[1];
sx q[1];
rz(-0.42196208) q[1];
sx q[1];
rz(1.8629927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.097825) q[0];
sx q[0];
rz(-1.861933) q[0];
sx q[0];
rz(-1.9242251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.638627) q[2];
sx q[2];
rz(-2.1178341) q[2];
sx q[2];
rz(-0.018195823) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9467612) q[1];
sx q[1];
rz(-1.166718) q[1];
sx q[1];
rz(-2.0003009) q[1];
rz(2.2845479) q[3];
sx q[3];
rz(-1.3856944) q[3];
sx q[3];
rz(-1.1777267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64855856) q[2];
sx q[2];
rz(-0.15668046) q[2];
sx q[2];
rz(-2.1873761) q[2];
rz(-0.75718015) q[3];
sx q[3];
rz(-1.4933519) q[3];
sx q[3];
rz(2.2455588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2467982) q[0];
sx q[0];
rz(-3.09258) q[0];
sx q[0];
rz(-3.0931296) q[0];
rz(2.5080644) q[1];
sx q[1];
rz(-2.3506479) q[1];
sx q[1];
rz(1.6652426) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7325981) q[0];
sx q[0];
rz(-2.2267003) q[0];
sx q[0];
rz(0.31314416) q[0];
x q[1];
rz(1.3858651) q[2];
sx q[2];
rz(-1.7586767) q[2];
sx q[2];
rz(0.056294346) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4004835) q[1];
sx q[1];
rz(-1.2874075) q[1];
sx q[1];
rz(0.41639911) q[1];
rz(-2.879056) q[3];
sx q[3];
rz(-1.3634014) q[3];
sx q[3];
rz(1.2174264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9897291) q[2];
sx q[2];
rz(-2.180474) q[2];
sx q[2];
rz(-0.1667008) q[2];
rz(-1.9549595) q[3];
sx q[3];
rz(-1.9125166) q[3];
sx q[3];
rz(-0.96773875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1879021) q[0];
sx q[0];
rz(-0.1583651) q[0];
sx q[0];
rz(-0.70575356) q[0];
rz(1.7965652) q[1];
sx q[1];
rz(-1.2860362) q[1];
sx q[1];
rz(2.8154624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46485422) q[0];
sx q[0];
rz(-2.2660672) q[0];
sx q[0];
rz(-0.37842964) q[0];
rz(-pi) q[1];
rz(-2.0847959) q[2];
sx q[2];
rz(-1.56034) q[2];
sx q[2];
rz(0.20960854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20777648) q[1];
sx q[1];
rz(-1.414555) q[1];
sx q[1];
rz(1.7396512) q[1];
x q[2];
rz(-2.9860849) q[3];
sx q[3];
rz(-1.9846791) q[3];
sx q[3];
rz(-0.83603102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8439794) q[2];
sx q[2];
rz(-1.16301) q[2];
sx q[2];
rz(0.45912099) q[2];
rz(-1.8504359) q[3];
sx q[3];
rz(-1.922311) q[3];
sx q[3];
rz(-0.084376924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5952263) q[0];
sx q[0];
rz(-1.3331174) q[0];
sx q[0];
rz(3.1072733) q[0];
rz(1.307084) q[1];
sx q[1];
rz(-0.75051761) q[1];
sx q[1];
rz(1.7274571) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33024055) q[0];
sx q[0];
rz(-0.76302397) q[0];
sx q[0];
rz(0.33154063) q[0];
rz(-pi) q[1];
rz(-1.1700045) q[2];
sx q[2];
rz(-2.199186) q[2];
sx q[2];
rz(-2.4054995) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88534158) q[1];
sx q[1];
rz(-2.4180545) q[1];
sx q[1];
rz(-2.4916562) q[1];
rz(1.565738) q[3];
sx q[3];
rz(-1.8190608) q[3];
sx q[3];
rz(-2.2170273) q[3];
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
rz(-0.59520477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0625192) q[0];
sx q[0];
rz(-0.74051028) q[0];
sx q[0];
rz(-0.24653521) q[0];
rz(1.1277699) q[1];
sx q[1];
rz(-2.0083387) q[1];
sx q[1];
rz(2.9182428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6205231) q[0];
sx q[0];
rz(-2.6436673) q[0];
sx q[0];
rz(0.070052043) q[0];
x q[1];
rz(-1.7511654) q[2];
sx q[2];
rz(-1.0374238) q[2];
sx q[2];
rz(-2.7636301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97572647) q[1];
sx q[1];
rz(-1.8687848) q[1];
sx q[1];
rz(-0.061640114) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8605883) q[3];
sx q[3];
rz(-1.9802571) q[3];
sx q[3];
rz(-1.0608761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8862306) q[2];
sx q[2];
rz(-0.88408154) q[2];
sx q[2];
rz(-3.0112126) q[2];
rz(-2.5731795) q[3];
sx q[3];
rz(-1.3610898) q[3];
sx q[3];
rz(-0.4140678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100621) q[0];
sx q[0];
rz(-0.93183403) q[0];
sx q[0];
rz(0.43781042) q[0];
rz(-1.9991649) q[1];
sx q[1];
rz(-1.2429968) q[1];
sx q[1];
rz(-1.4243781) q[1];
rz(0.42752888) q[2];
sx q[2];
rz(-1.1784381) q[2];
sx q[2];
rz(0.094712275) q[2];
rz(-1.891248) q[3];
sx q[3];
rz(-1.9449825) q[3];
sx q[3];
rz(-2.2016761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
