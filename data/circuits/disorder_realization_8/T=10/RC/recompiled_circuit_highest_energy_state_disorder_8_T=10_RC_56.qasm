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
rz(5.4665718) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(10.583034) q[0];
rz(0.090016063) q[1];
sx q[1];
rz(-2.6675192) q[1];
sx q[1];
rz(-0.83516821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4157566) q[0];
sx q[0];
rz(-0.8716363) q[0];
sx q[0];
rz(1.1582202) q[0];
rz(-2.5005241) q[2];
sx q[2];
rz(-1.3136697) q[2];
sx q[2];
rz(-0.4349622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.07552927) q[1];
sx q[1];
rz(-1.8154241) q[1];
sx q[1];
rz(0.3800769) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6938758) q[3];
sx q[3];
rz(-0.78842794) q[3];
sx q[3];
rz(-0.83886787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-2.2455402) q[2];
sx q[2];
rz(-0.33207616) q[2];
rz(-2.8927228) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2700972) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(-0.53263295) q[0];
rz(0.10781413) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(-2.8210988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924393) q[0];
sx q[0];
rz(-2.4863613) q[0];
sx q[0];
rz(-1.4661319) q[0];
rz(2.2079289) q[2];
sx q[2];
rz(-0.28544989) q[2];
sx q[2];
rz(0.97412005) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1793666) q[1];
sx q[1];
rz(-1.4186064) q[1];
sx q[1];
rz(2.9801912) q[1];
x q[2];
rz(-0.56902253) q[3];
sx q[3];
rz(-2.7561765) q[3];
sx q[3];
rz(-0.33447124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33114854) q[2];
sx q[2];
rz(-2.836477) q[2];
sx q[2];
rz(-1.6395052) q[2];
rz(2.8034927) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(0.52687183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51760393) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(-2.7841618) q[0];
rz(-2.7492211) q[1];
sx q[1];
rz(-2.3452499) q[1];
sx q[1];
rz(1.9270814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93335184) q[0];
sx q[0];
rz(-0.14259556) q[0];
sx q[0];
rz(1.4487793) q[0];
rz(0.54698555) q[2];
sx q[2];
rz(-1.437709) q[2];
sx q[2];
rz(-2.2625429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9704772) q[1];
sx q[1];
rz(-1.5737281) q[1];
sx q[1];
rz(1.9582002) q[1];
rz(-pi) q[2];
rz(0.9830832) q[3];
sx q[3];
rz(-1.6020892) q[3];
sx q[3];
rz(0.33558095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3452722) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(2.7447682) q[2];
rz(-1.2567629) q[3];
sx q[3];
rz(-1.4182914) q[3];
sx q[3];
rz(2.5313012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20763718) q[0];
sx q[0];
rz(-1.3960681) q[0];
sx q[0];
rz(1.2285832) q[0];
rz(-1.2127016) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(-0.62087762) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29781326) q[0];
sx q[0];
rz(-0.50619805) q[0];
sx q[0];
rz(1.9202581) q[0];
rz(-pi) q[1];
rz(0.75350301) q[2];
sx q[2];
rz(-2.916159) q[2];
sx q[2];
rz(-1.6627251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3132857) q[1];
sx q[1];
rz(-1.8290214) q[1];
sx q[1];
rz(-0.26744803) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6895164) q[3];
sx q[3];
rz(-1.9510036) q[3];
sx q[3];
rz(-0.79352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2584201) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(2.9849198) q[2];
rz(2.2251718) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85585344) q[0];
sx q[0];
rz(-1.9802977) q[0];
sx q[0];
rz(2.7440985) q[0];
rz(-3.1187348) q[1];
sx q[1];
rz(-0.50067478) q[1];
sx q[1];
rz(2.4668677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0521654) q[0];
sx q[0];
rz(-2.8244655) q[0];
sx q[0];
rz(-2.2990312) q[0];
rz(1.4762474) q[2];
sx q[2];
rz(-1.4473905) q[2];
sx q[2];
rz(-2.9837554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0583084) q[1];
sx q[1];
rz(-1.7232019) q[1];
sx q[1];
rz(0.17249523) q[1];
rz(-pi) q[2];
rz(1.8800903) q[3];
sx q[3];
rz(-1.3095244) q[3];
sx q[3];
rz(2.3051534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(-2.442339) q[2];
rz(-1.3462542) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(0.7114555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(-2.7432192) q[0];
rz(0.97459546) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(-1.7030565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4476623) q[0];
sx q[0];
rz(-1.2691783) q[0];
sx q[0];
rz(-1.680879) q[0];
rz(-pi) q[1];
rz(-1.6104524) q[2];
sx q[2];
rz(-0.25601124) q[2];
sx q[2];
rz(-2.2289952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7185091) q[1];
sx q[1];
rz(-0.72186493) q[1];
sx q[1];
rz(-1.2150498) q[1];
rz(-pi) q[2];
rz(-1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(-1.069266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4633816) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(1.8360651) q[2];
rz(1.8259004) q[3];
sx q[3];
rz(-1.3633599) q[3];
sx q[3];
rz(-0.8684043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618133) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(1.4965936) q[0];
rz(-2.1553701) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(-0.49759069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80579306) q[0];
sx q[0];
rz(-1.8247274) q[0];
sx q[0];
rz(-3.0428314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6025101) q[2];
sx q[2];
rz(-1.0175127) q[2];
sx q[2];
rz(-2.6607571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0329602) q[1];
sx q[1];
rz(-1.7268306) q[1];
sx q[1];
rz(-1.1724657) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4769745) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(2.97992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8387973) q[2];
sx q[2];
rz(-1.5668198) q[2];
sx q[2];
rz(-0.29067579) q[2];
rz(0.21026462) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(1.0342342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656089) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(1.959311) q[0];
rz(-0.15444175) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(2.2363037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2700896) q[0];
sx q[0];
rz(-0.2729899) q[0];
sx q[0];
rz(1.6156625) q[0];
x q[1];
rz(-2.6851875) q[2];
sx q[2];
rz(-1.3564566) q[2];
sx q[2];
rz(2.6021007) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0539848) q[1];
sx q[1];
rz(-1.4248669) q[1];
sx q[1];
rz(-2.3838359) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0175627) q[3];
sx q[3];
rz(-0.39504566) q[3];
sx q[3];
rz(1.0329122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0514544) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(2.1913989) q[2];
rz(3.0692302) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(-2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41910928) q[0];
sx q[0];
rz(-2.1868732) q[0];
sx q[0];
rz(2.0461653) q[0];
rz(1.0911881) q[1];
sx q[1];
rz(-2.2406816) q[1];
sx q[1];
rz(0.99517623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93928775) q[0];
sx q[0];
rz(-1.5065333) q[0];
sx q[0];
rz(2.4758384) q[0];
rz(-pi) q[1];
rz(-1.9068933) q[2];
sx q[2];
rz(-1.8310412) q[2];
sx q[2];
rz(-1.1059831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51982626) q[1];
sx q[1];
rz(-1.4859746) q[1];
sx q[1];
rz(1.6165401) q[1];
rz(0.18539683) q[3];
sx q[3];
rz(-1.3747921) q[3];
sx q[3];
rz(-0.40962266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5114078) q[2];
sx q[2];
rz(-2.6738622) q[2];
sx q[2];
rz(0.67997813) q[2];
rz(-1.4558815) q[3];
sx q[3];
rz(-0.89434353) q[3];
sx q[3];
rz(-1.1792012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7162914) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(0.61538482) q[0];
rz(-1.3353434) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(-1.7937484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772661) q[0];
sx q[0];
rz(-1.6166256) q[0];
sx q[0];
rz(2.7861676) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87839076) q[2];
sx q[2];
rz(-2.2579402) q[2];
sx q[2];
rz(-2.2855482) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19619689) q[1];
sx q[1];
rz(-0.41204231) q[1];
sx q[1];
rz(0.41007385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12935454) q[3];
sx q[3];
rz(-0.760303) q[3];
sx q[3];
rz(-0.018825104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24796692) q[2];
sx q[2];
rz(-1.3004356) q[2];
sx q[2];
rz(-0.51353961) q[2];
rz(2.9945471) q[3];
sx q[3];
rz(-0.5439609) q[3];
sx q[3];
rz(-1.2646382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87179398) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(-1.3006032) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(-1.945449) q[2];
sx q[2];
rz(-1.7032663) q[2];
sx q[2];
rz(1.8864529) q[2];
rz(2.5138598) q[3];
sx q[3];
rz(-1.9398324) q[3];
sx q[3];
rz(2.0331665) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
