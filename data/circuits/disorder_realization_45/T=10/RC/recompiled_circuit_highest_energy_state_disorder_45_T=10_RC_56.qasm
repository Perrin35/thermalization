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
rz(-1.2019914) q[0];
sx q[0];
rz(-2.658598) q[0];
sx q[0];
rz(1.510409) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18099526) q[0];
sx q[0];
rz(-2.1707524) q[0];
sx q[0];
rz(2.9636895) q[0];
rz(-pi) q[1];
rz(0.59487409) q[2];
sx q[2];
rz(-1.383179) q[2];
sx q[2];
rz(1.7930195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12355676) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(-0.058079795) q[1];
rz(-pi) q[2];
rz(-2.9128051) q[3];
sx q[3];
rz(-1.2902033) q[3];
sx q[3];
rz(-0.56786637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55945021) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-0.5994125) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(0.74554602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6744989) q[0];
sx q[0];
rz(-1.281453) q[0];
sx q[0];
rz(1.4379005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7600766) q[2];
sx q[2];
rz(-1.3323297) q[2];
sx q[2];
rz(0.60205215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9245431) q[1];
sx q[1];
rz(-2.1991208) q[1];
sx q[1];
rz(-1.373686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.638444) q[3];
sx q[3];
rz(-2.1981578) q[3];
sx q[3];
rz(-2.7706551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45592371) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.7342742) q[2];
rz(-2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(-2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5582964) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(2.1287647) q[0];
rz(2.5335675) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(1.2695405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69508775) q[0];
sx q[0];
rz(-1.1881314) q[0];
sx q[0];
rz(1.0951359) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082398675) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(-1.5562039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0996272) q[1];
sx q[1];
rz(-1.3746975) q[1];
sx q[1];
rz(0.72814299) q[1];
rz(-pi) q[2];
rz(2.5699411) q[3];
sx q[3];
rz(-2.5749675) q[3];
sx q[3];
rz(-1.0039312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6262007) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(-2.9279809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13852791) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277058) q[0];
sx q[0];
rz(-1.9235652) q[0];
sx q[0];
rz(1.1350495) q[0];
x q[1];
rz(2.1221913) q[2];
sx q[2];
rz(-1.2605091) q[2];
sx q[2];
rz(-2.8648368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5508073) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(0.6145668) q[1];
rz(-pi) q[2];
rz(2.6163231) q[3];
sx q[3];
rz(-2.161918) q[3];
sx q[3];
rz(-1.0159462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5248096) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(-1.3516124) q[2];
rz(2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(-1.9701689) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(3.0426262) q[0];
rz(-2.4348266) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(0.4695355) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564559) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(1.3993457) q[0];
x q[1];
rz(2.7623873) q[2];
sx q[2];
rz(-2.0742886) q[2];
sx q[2];
rz(-0.45631726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81988702) q[1];
sx q[1];
rz(-2.4155136) q[1];
sx q[1];
rz(-0.5989845) q[1];
rz(-pi) q[2];
rz(2.9088777) q[3];
sx q[3];
rz(-2.5732627) q[3];
sx q[3];
rz(-1.9459023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(2.7169054) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(-0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.11923085) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(-1.7599531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19162543) q[0];
sx q[0];
rz(-0.81474308) q[0];
sx q[0];
rz(-0.11884584) q[0];
x q[1];
rz(-2.6082615) q[2];
sx q[2];
rz(-1.0600277) q[2];
sx q[2];
rz(0.086670808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7834251) q[1];
sx q[1];
rz(-1.3853711) q[1];
sx q[1];
rz(-2.3802451) q[1];
rz(-pi) q[2];
rz(-1.6448037) q[3];
sx q[3];
rz(-1.2701708) q[3];
sx q[3];
rz(-3.0174676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(2.6324658) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(3.0522108) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.4559454) q[1];
sx q[1];
rz(-2.1481029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8315961) q[0];
sx q[0];
rz(-0.78959268) q[0];
sx q[0];
rz(1.5168651) q[0];
rz(-pi) q[1];
rz(-0.69665945) q[2];
sx q[2];
rz(-2.2242745) q[2];
sx q[2];
rz(-2.4652664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3741039) q[1];
sx q[1];
rz(-1.933177) q[1];
sx q[1];
rz(2.9844173) q[1];
rz(0.66331373) q[3];
sx q[3];
rz(-2.3938826) q[3];
sx q[3];
rz(1.6689036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(-0.40204027) q[2];
rz(-2.0461931) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(-0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47356975) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(1.8079669) q[0];
rz(2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(-0.91748253) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79865341) q[0];
sx q[0];
rz(-0.22428939) q[0];
sx q[0];
rz(2.0464315) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94709227) q[2];
sx q[2];
rz(-1.6674041) q[2];
sx q[2];
rz(1.7283224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0304347) q[1];
sx q[1];
rz(-1.7923755) q[1];
sx q[1];
rz(1.9244003) q[1];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.9894674) q[3];
sx q[3];
rz(0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(1.8514006) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(-0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(0.27715096) q[0];
rz(1.9174891) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(-1.3714429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087589892) q[0];
sx q[0];
rz(-1.3773019) q[0];
sx q[0];
rz(-0.67275472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4898289) q[2];
sx q[2];
rz(-2.4874176) q[2];
sx q[2];
rz(0.68133611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0321694) q[1];
sx q[1];
rz(-0.9309097) q[1];
sx q[1];
rz(2.8830322) q[1];
rz(-pi) q[2];
rz(2.9670466) q[3];
sx q[3];
rz(-2.2851599) q[3];
sx q[3];
rz(-2.3256231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93854967) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(0.31442434) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.2615874) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(-0.57149291) q[0];
rz(2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(-2.4933955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2637973) q[0];
sx q[0];
rz(-1.5938984) q[0];
sx q[0];
rz(-0.013851555) q[0];
rz(2.4785751) q[2];
sx q[2];
rz(-0.94919357) q[2];
sx q[2];
rz(-0.66689516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98013377) q[1];
sx q[1];
rz(-2.4422283) q[1];
sx q[1];
rz(-3.1070263) q[1];
rz(-0.51697124) q[3];
sx q[3];
rz(-1.3992157) q[3];
sx q[3];
rz(-1.7652546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8615243) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(-1.466922) q[2];
rz(-2.9649949) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27547729) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(2.7571309) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(2.4118829) q[2];
sx q[2];
rz(-2.6860102) q[2];
sx q[2];
rz(-0.44853733) q[2];
rz(1.0064784) q[3];
sx q[3];
rz(-2.0340393) q[3];
sx q[3];
rz(-1.5599193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
