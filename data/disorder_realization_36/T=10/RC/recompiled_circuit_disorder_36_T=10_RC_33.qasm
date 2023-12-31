OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(0.51529348) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(2.5640092) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-2.8773017) q[0];
rz(-2.0787813) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(-2.9302772) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88749332) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(-1.9894132) q[1];
rz(-pi) q[2];
rz(-0.76831423) q[3];
sx q[3];
rz(-0.81211219) q[3];
sx q[3];
rz(2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(-0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0613522) q[0];
sx q[0];
rz(-0.70525673) q[0];
sx q[0];
rz(0.7028701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3697853) q[2];
sx q[2];
rz(-2.9559921) q[2];
sx q[2];
rz(1.1112569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3413275) q[1];
sx q[1];
rz(-2.6086573) q[1];
sx q[1];
rz(-0.37377263) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94988471) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(1.8267531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(1.0864331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(3.1397318) q[0];
x q[1];
rz(-2.2181182) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(2.0331241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74124386) q[1];
sx q[1];
rz(-1.2819918) q[1];
sx q[1];
rz(1.5653533) q[1];
rz(-pi) q[2];
rz(1.9757189) q[3];
sx q[3];
rz(-1.7305264) q[3];
sx q[3];
rz(-0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(0.80254054) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.6960467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(-0.66288373) q[0];
rz(2.9252058) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(-0.83425922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2913937) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.5047969) q[1];
rz(2.8533832) q[3];
sx q[3];
rz(-2.0734348) q[3];
sx q[3];
rz(2.1484745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4399453) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-2.2873926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094164205) q[0];
sx q[0];
rz(-1.5625956) q[0];
sx q[0];
rz(-2.2846983) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1335667) q[2];
sx q[2];
rz(-2.3587583) q[2];
sx q[2];
rz(2.7340739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90480587) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(-0.76101117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0075931) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(-1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(0.577315) q[2];
rz(-0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-3.0153826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835411) q[0];
sx q[0];
rz(-1.7833033) q[0];
sx q[0];
rz(1.538518) q[0];
rz(-pi) q[1];
rz(2.7563165) q[2];
sx q[2];
rz(-2.3284973) q[2];
sx q[2];
rz(0.77020459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8773552) q[1];
sx q[1];
rz(-0.74062956) q[1];
sx q[1];
rz(-1.8514368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.24228) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(1.1920284) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647117) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-0.60140951) q[0];
x q[1];
rz(0.26142188) q[2];
sx q[2];
rz(-1.520642) q[2];
sx q[2];
rz(0.23471552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(-0.97538235) q[1];
rz(1.7667174) q[3];
sx q[3];
rz(-0.94983263) q[3];
sx q[3];
rz(1.2697112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(0.34379488) q[2];
rz(2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5430785) q[0];
sx q[0];
rz(-3.0656272) q[0];
sx q[0];
rz(-0.11673467) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7224738) q[2];
sx q[2];
rz(-1.1693839) q[2];
sx q[2];
rz(-2.7733208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1069113) q[1];
sx q[1];
rz(-2.2195663) q[1];
sx q[1];
rz(-1.1120863) q[1];
x q[2];
rz(3.001508) q[3];
sx q[3];
rz(-2.1685371) q[3];
sx q[3];
rz(2.7261581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(0.16601673) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.5244012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7302007) q[0];
sx q[0];
rz(-3.054266) q[0];
sx q[0];
rz(1.9090396) q[0];
x q[1];
rz(2.6644601) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(2.0810623) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2251687) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(-1.9041063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7306261) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(-0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2828335) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(1.1742924) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9099405) q[2];
sx q[2];
rz(-1.7288496) q[2];
sx q[2];
rz(-0.16976419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(2.5245689) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57755034) q[3];
sx q[3];
rz(-0.69637075) q[3];
sx q[3];
rz(1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-2.3989427) q[2];
sx q[2];
rz(-2.7682318) q[2];
sx q[2];
rz(-2.9329185) q[2];
rz(-0.75541227) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
