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
rz(-1.5440829) q[0];
sx q[0];
rz(4.5151526) q[0];
sx q[0];
rz(8.406352) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(0.63176027) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951183) q[0];
sx q[0];
rz(-0.24226878) q[0];
sx q[0];
rz(1.0231859) q[0];
rz(-pi) q[1];
rz(0.024876923) q[2];
sx q[2];
rz(-1.8734249) q[2];
sx q[2];
rz(2.6278815) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5093928) q[1];
sx q[1];
rz(-2.191847) q[1];
sx q[1];
rz(-2.7953447) q[1];
rz(-3.0491563) q[3];
sx q[3];
rz(-1.8957924) q[3];
sx q[3];
rz(2.6324038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.53039256) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(-1.4872888) q[2];
rz(-0.90855956) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1210043) q[0];
sx q[0];
rz(-1.4326743) q[0];
sx q[0];
rz(-2.9793136) q[0];
rz(-0.24457112) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(2.8499106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8957386) q[0];
sx q[0];
rz(-1.3783921) q[0];
sx q[0];
rz(-2.8537575) q[0];
x q[1];
rz(2.2333916) q[2];
sx q[2];
rz(-2.78763) q[2];
sx q[2];
rz(-2.5733054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6704935) q[1];
sx q[1];
rz(-2.3355977) q[1];
sx q[1];
rz(-0.32260311) q[1];
x q[2];
rz(1.1532182) q[3];
sx q[3];
rz(-1.1912701) q[3];
sx q[3];
rz(0.97009995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67819277) q[2];
sx q[2];
rz(-0.86898154) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(-1.894527) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(-1.6943078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70188824) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(-0.18371789) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(2.9768129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8863109) q[0];
sx q[0];
rz(-1.5802529) q[0];
sx q[0];
rz(2.1566118) q[0];
x q[1];
rz(-2.7585331) q[2];
sx q[2];
rz(-2.428741) q[2];
sx q[2];
rz(1.0733611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2307869) q[1];
sx q[1];
rz(-1.9525098) q[1];
sx q[1];
rz(-0.64818212) q[1];
rz(-pi) q[2];
rz(-1.0298877) q[3];
sx q[3];
rz(-0.60324429) q[3];
sx q[3];
rz(-1.3689976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(1.1064233) q[2];
rz(0.70118457) q[3];
sx q[3];
rz(-1.8132352) q[3];
sx q[3];
rz(2.6760694) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(-2.1970774) q[0];
rz(1.5185897) q[1];
sx q[1];
rz(-2.1606162) q[1];
sx q[1];
rz(-1.5400344) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0991411) q[0];
sx q[0];
rz(-0.97223982) q[0];
sx q[0];
rz(0.70229806) q[0];
rz(2.3643199) q[2];
sx q[2];
rz(-2.9748355) q[2];
sx q[2];
rz(0.47698944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67805144) q[1];
sx q[1];
rz(-1.9423331) q[1];
sx q[1];
rz(2.3080565) q[1];
rz(-pi) q[2];
rz(-3.0936315) q[3];
sx q[3];
rz(-2.3905919) q[3];
sx q[3];
rz(-2.6329071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1314142) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(0.16769257) q[2];
rz(0.0532648) q[3];
sx q[3];
rz(-1.1054509) q[3];
sx q[3];
rz(1.3767327) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(-2.8422624) q[0];
rz(-0.37295595) q[1];
sx q[1];
rz(-1.8920218) q[1];
sx q[1];
rz(-1.6711055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0051992) q[0];
sx q[0];
rz(-1.0233425) q[0];
sx q[0];
rz(-1.0639079) q[0];
x q[1];
rz(-1.4506571) q[2];
sx q[2];
rz(-1.1098301) q[2];
sx q[2];
rz(1.9751939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59488397) q[1];
sx q[1];
rz(-2.1404033) q[1];
sx q[1];
rz(2.1807266) q[1];
x q[2];
rz(-1.3277131) q[3];
sx q[3];
rz(-2.5057) q[3];
sx q[3];
rz(-2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9386998) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(-1.1104442) q[2];
rz(0.86197305) q[3];
sx q[3];
rz(-1.6317261) q[3];
sx q[3];
rz(-1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2170169) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(0.4050912) q[0];
rz(-2.8054667) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(2.1902671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66621214) q[0];
sx q[0];
rz(-0.48298353) q[0];
sx q[0];
rz(-0.57421143) q[0];
x q[1];
rz(-1.7826005) q[2];
sx q[2];
rz(-2.2957605) q[2];
sx q[2];
rz(2.447809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8321633) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(1.840074) q[1];
x q[2];
rz(-0.27928593) q[3];
sx q[3];
rz(-1.3853867) q[3];
sx q[3];
rz(-0.11073555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0508017) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(2.9108099) q[2];
rz(-0.18203059) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35293216) q[0];
sx q[0];
rz(-3.1020628) q[0];
sx q[0];
rz(-0.93210644) q[0];
rz(-2.508029) q[1];
sx q[1];
rz(-1.1195868) q[1];
sx q[1];
rz(-1.3379785) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18620488) q[0];
sx q[0];
rz(-3.0078631) q[0];
sx q[0];
rz(0.9132847) q[0];
x q[1];
rz(2.8678738) q[2];
sx q[2];
rz(-1.459957) q[2];
sx q[2];
rz(2.1901166) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4949957) q[1];
sx q[1];
rz(-2.3892168) q[1];
sx q[1];
rz(0.97480358) q[1];
rz(-0.058480992) q[3];
sx q[3];
rz(-2.6425458) q[3];
sx q[3];
rz(2.8935695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6771217) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(-1.202549) q[2];
rz(2.7770212) q[3];
sx q[3];
rz(-0.69260827) q[3];
sx q[3];
rz(-0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637591) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(-0.17663503) q[0];
rz(0.44003507) q[1];
sx q[1];
rz(-1.9562079) q[1];
sx q[1];
rz(0.96493351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595183) q[0];
sx q[0];
rz(-1.4384369) q[0];
sx q[0];
rz(2.9455723) q[0];
rz(-pi) q[1];
rz(-2.1148483) q[2];
sx q[2];
rz(-1.0389345) q[2];
sx q[2];
rz(-0.085467664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1159462) q[1];
sx q[1];
rz(-1.0315511) q[1];
sx q[1];
rz(0.42614062) q[1];
rz(-pi) q[2];
x q[2];
rz(0.075841622) q[3];
sx q[3];
rz(-1.3170529) q[3];
sx q[3];
rz(-0.50204078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21260103) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(-3.0250004) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(-2.5998083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56124878) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(-0.62193459) q[0];
rz(1.7033345) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(-0.66351801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7166876) q[0];
sx q[0];
rz(-0.7640673) q[0];
sx q[0];
rz(0.31830799) q[0];
rz(-pi) q[1];
rz(-1.2873257) q[2];
sx q[2];
rz(-1.6890235) q[2];
sx q[2];
rz(-0.17519874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89412824) q[1];
sx q[1];
rz(-0.93527764) q[1];
sx q[1];
rz(2.7429917) q[1];
rz(-pi) q[2];
rz(-2.8264753) q[3];
sx q[3];
rz(-0.2956008) q[3];
sx q[3];
rz(-0.35741266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62349391) q[2];
sx q[2];
rz(-1.7524717) q[2];
sx q[2];
rz(-0.36250472) q[2];
rz(0.47719657) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(-1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057673205) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(2.2398563) q[0];
rz(2.4665191) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(-2.9973082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5675303) q[0];
sx q[0];
rz(-1.6533378) q[0];
sx q[0];
rz(-1.922419) q[0];
x q[1];
rz(-1.6727757) q[2];
sx q[2];
rz(-1.7893409) q[2];
sx q[2];
rz(-0.25221014) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.044202494) q[1];
sx q[1];
rz(-2.4501738) q[1];
sx q[1];
rz(2.7425062) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7561109) q[3];
sx q[3];
rz(-0.47244888) q[3];
sx q[3];
rz(2.2466898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6361864) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(-2.9409161) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.5526445) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(0.48925346) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(-0.092838661) q[2];
sx q[2];
rz(-1.9839109) q[2];
sx q[2];
rz(-2.4674923) q[2];
rz(-3.0808385) q[3];
sx q[3];
rz(-0.44296064) q[3];
sx q[3];
rz(0.62173494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
