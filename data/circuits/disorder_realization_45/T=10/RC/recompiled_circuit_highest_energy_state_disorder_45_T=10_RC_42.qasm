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
rz(3.6245873) q[0];
sx q[0];
rz(10.935187) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529634) q[0];
sx q[0];
rz(-1.7173816) q[0];
sx q[0];
rz(0.96340553) q[0];
x q[1];
rz(1.3454622) q[2];
sx q[2];
rz(-0.98774922) q[2];
sx q[2];
rz(-3.0449113) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3986721) q[1];
sx q[1];
rz(-1.5389331) q[1];
sx q[1];
rz(-2.5614061) q[1];
rz(-pi) q[2];
rz(-1.2830986) q[3];
sx q[3];
rz(-1.7904864) q[3];
sx q[3];
rz(0.9385329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1175179) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(-0.81895858) q[2];
rz(3.135318) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55945021) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-2.5421802) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760147) q[0];
sx q[0];
rz(-1.6981372) q[0];
sx q[0];
rz(-0.29177427) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2471871) q[2];
sx q[2];
rz(-2.3043046) q[2];
sx q[2];
rz(-2.3937283) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47059575) q[1];
sx q[1];
rz(-1.4116916) q[1];
sx q[1];
rz(-2.5039423) q[1];
x q[2];
rz(2.2621731) q[3];
sx q[3];
rz(-1.169765) q[3];
sx q[3];
rz(-1.629231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(-1.4073184) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(-1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5832962) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(1.012828) q[0];
rz(2.5335675) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(1.8720522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2484754) q[0];
sx q[0];
rz(-2.5404544) q[0];
sx q[0];
rz(-0.84971835) q[0];
rz(-1.467134) q[2];
sx q[2];
rz(-2.4680063) q[2];
sx q[2];
rz(1.4530593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3568282) q[1];
sx q[1];
rz(-0.85961378) q[1];
sx q[1];
rz(-1.8309092) q[1];
rz(2.6502887) q[3];
sx q[3];
rz(-1.2761371) q[3];
sx q[3];
rz(-0.06959411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(3.1332704) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(-1.3767161) q[0];
rz(1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6277058) q[0];
sx q[0];
rz(-1.2180274) q[0];
sx q[0];
rz(-2.0065432) q[0];
rz(-pi) q[1];
rz(-1.0194014) q[2];
sx q[2];
rz(-1.2605091) q[2];
sx q[2];
rz(-2.8648368) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59078539) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(2.5270259) q[1];
rz(-2.6163231) q[3];
sx q[3];
rz(-0.97967463) q[3];
sx q[3];
rz(-1.0159462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6167831) q[2];
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
rz(-pi) q[1];
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
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(0.098966448) q[0];
rz(-0.70676604) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(2.6720572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75824518) q[0];
sx q[0];
rz(-1.7086141) q[0];
sx q[0];
rz(3.1175749) q[0];
rz(-pi) q[1];
rz(-2.1060313) q[2];
sx q[2];
rz(-1.9010086) q[2];
sx q[2];
rz(-0.92452985) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5818122) q[1];
sx q[1];
rz(-0.99039927) q[1];
sx q[1];
rz(2.0349166) q[1];
x q[2];
rz(-0.55600663) q[3];
sx q[3];
rz(-1.4463498) q[3];
sx q[3];
rz(0.57226244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2609451) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(1.3060695) q[2];
rz(-2.7169054) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(-1.3223883) q[0];
rz(2.0139096) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.3816396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36395006) q[0];
sx q[0];
rz(-2.3780883) q[0];
sx q[0];
rz(-1.6958773) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99389771) q[2];
sx q[2];
rz(-1.111278) q[2];
sx q[2];
rz(-1.3764868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35816757) q[1];
sx q[1];
rz(-1.3853711) q[1];
sx q[1];
rz(-0.76134759) q[1];
rz(-0.30140169) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(1.4247198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7606925) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(-2.6324658) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(3.0522108) q[0];
rz(2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(0.99348974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334796) q[0];
sx q[0];
rz(-0.78266875) q[0];
sx q[0];
rz(0.054305768) q[0];
rz(-pi) q[1];
rz(2.3553576) q[2];
sx q[2];
rz(-2.1053542) q[2];
sx q[2];
rz(1.7766118) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.74718432) q[1];
sx q[1];
rz(-1.7176873) q[1];
sx q[1];
rz(-1.2042852) q[1];
x q[2];
rz(-2.4782789) q[3];
sx q[3];
rz(-0.74771008) q[3];
sx q[3];
rz(1.472689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(2.0461931) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680229) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(-1.3336257) q[0];
rz(-2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(-2.2241101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30675754) q[0];
sx q[0];
rz(-1.6728171) q[0];
sx q[0];
rz(-1.3707042) q[0];
rz(-pi) q[1];
rz(-1.4063666) q[2];
sx q[2];
rz(-0.63015579) q[2];
sx q[2];
rz(3.117331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6010103) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(0.2356727) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(2.6961117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(1.290192) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057587) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(-1.2241036) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(-1.3714429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4215721) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(-2.8369342) q[0];
rz(-pi) q[1];
rz(2.5941237) q[2];
sx q[2];
rz(-1.9488504) q[2];
sx q[2];
rz(-2.7965656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.692457) q[1];
sx q[1];
rz(-0.68329158) q[1];
sx q[1];
rz(1.9016674) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84884642) q[3];
sx q[3];
rz(-1.4392142) q[3];
sx q[3];
rz(-2.2717486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.203043) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(-0.68823632) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(1.8800053) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(0.64819711) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.307319) q[0];
sx q[0];
rz(-1.5569485) q[0];
sx q[0];
rz(1.547692) q[0];
x q[1];
rz(2.3085015) q[2];
sx q[2];
rz(-1.0467741) q[2];
sx q[2];
rz(1.3307216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1163016) q[1];
sx q[1];
rz(-2.2696583) q[1];
sx q[1];
rz(1.5998597) q[1];
rz(-pi) q[2];
rz(-0.33721029) q[3];
sx q[3];
rz(-2.5993532) q[3];
sx q[3];
rz(0.48619871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28006831) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(2.9649949) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(1.2944029) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(2.7914417) q[2];
sx q[2];
rz(-1.8684917) q[2];
sx q[2];
rz(-2.6960052) q[2];
rz(-2.3220358) q[3];
sx q[3];
rz(-2.4278276) q[3];
sx q[3];
rz(0.62558382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
