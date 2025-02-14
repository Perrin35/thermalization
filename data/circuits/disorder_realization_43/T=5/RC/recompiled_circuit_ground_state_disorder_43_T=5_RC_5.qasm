OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(2.6743968) q[0];
sx q[0];
rz(13.462486) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2786699) q[0];
sx q[0];
rz(-0.76919829) q[0];
sx q[0];
rz(1.0095566) q[0];
rz(2.5619123) q[2];
sx q[2];
rz(-2.2319921) q[2];
sx q[2];
rz(1.5722317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7973229) q[1];
sx q[1];
rz(-1.8472278) q[1];
sx q[1];
rz(1.396342) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39405502) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(1.8436028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2200615) q[2];
sx q[2];
rz(-1.0342197) q[2];
sx q[2];
rz(-2.6424778) q[2];
rz(2.3257997) q[3];
sx q[3];
rz(-0.59320265) q[3];
sx q[3];
rz(0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4895184) q[0];
sx q[0];
rz(-0.3638142) q[0];
sx q[0];
rz(2.9079085) q[0];
rz(3.0417327) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(-1.608009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31029345) q[0];
sx q[0];
rz(-1.6673894) q[0];
sx q[0];
rz(1.4486758) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38012114) q[2];
sx q[2];
rz(-0.46762662) q[2];
sx q[2];
rz(-1.1884226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93615426) q[1];
sx q[1];
rz(-0.9673276) q[1];
sx q[1];
rz(-0.11118576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30939524) q[3];
sx q[3];
rz(-1.1077048) q[3];
sx q[3];
rz(2.9675296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0443772) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(0.054917939) q[2];
rz(-0.6461668) q[3];
sx q[3];
rz(-1.61444) q[3];
sx q[3];
rz(-3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.9570479) q[0];
sx q[0];
rz(-2.8962729) q[0];
sx q[0];
rz(0.74485892) q[0];
rz(1.8648196) q[1];
sx q[1];
rz(-2.1121912) q[1];
sx q[1];
rz(2.2778817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833763) q[0];
sx q[0];
rz(-1.5124707) q[0];
sx q[0];
rz(-0.44022473) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3529007) q[2];
sx q[2];
rz(-1.9448819) q[2];
sx q[2];
rz(-1.7497334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5472467) q[1];
sx q[1];
rz(-0.9556831) q[1];
sx q[1];
rz(-2.1662461) q[1];
x q[2];
rz(2.0710122) q[3];
sx q[3];
rz(-2.1693834) q[3];
sx q[3];
rz(2.9974496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6700217) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(-2.7460597) q[2];
rz(1.0223355) q[3];
sx q[3];
rz(-0.46415713) q[3];
sx q[3];
rz(0.49183229) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96524298) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(-0.055971948) q[0];
rz(-2.5296027) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(-2.2956119) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2200945) q[0];
sx q[0];
rz(-1.3125696) q[0];
sx q[0];
rz(1.0357712) q[0];
rz(-pi) q[1];
rz(2.1899784) q[2];
sx q[2];
rz(-1.2834833) q[2];
sx q[2];
rz(3.0653022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8236134) q[1];
sx q[1];
rz(-2.7732964) q[1];
sx q[1];
rz(2.5846534) q[1];
x q[2];
rz(0.37644551) q[3];
sx q[3];
rz(-1.1247702) q[3];
sx q[3];
rz(-0.50807398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.446898) q[2];
sx q[2];
rz(2.4212627) q[2];
rz(-0.35308853) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(-2.8928355) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.4530797) q[0];
sx q[0];
rz(-0.98689669) q[0];
rz(1.997021) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(-0.95485895) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3564044) q[0];
sx q[0];
rz(-2.4739728) q[0];
sx q[0];
rz(1.5285399) q[0];
x q[1];
rz(2.8138732) q[2];
sx q[2];
rz(-1.0340889) q[2];
sx q[2];
rz(-2.2959054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2209619) q[1];
sx q[1];
rz(-2.536533) q[1];
sx q[1];
rz(0.69843881) q[1];
rz(-pi) q[2];
rz(-2.0935861) q[3];
sx q[3];
rz(-2.8510333) q[3];
sx q[3];
rz(-2.235242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5130875) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(2.0950192) q[2];
rz(0.70932499) q[3];
sx q[3];
rz(-1.1449292) q[3];
sx q[3];
rz(-2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.3777622) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(-0.27107006) q[0];
rz(0.079004869) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(-1.7105506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56667152) q[0];
sx q[0];
rz(-1.780176) q[0];
sx q[0];
rz(-1.4230595) q[0];
x q[1];
rz(2.1151401) q[2];
sx q[2];
rz(-1.0751131) q[2];
sx q[2];
rz(3.0718671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9985314) q[1];
sx q[1];
rz(-2.240671) q[1];
sx q[1];
rz(-1.8865442) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9810199) q[3];
sx q[3];
rz(-1.1571006) q[3];
sx q[3];
rz(-0.90908066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.048910353) q[2];
sx q[2];
rz(-1.1834669) q[2];
sx q[2];
rz(-0.023822039) q[2];
rz(-1.228099) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(0.14321271) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.989885) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(0.09224961) q[0];
rz(0.30091885) q[1];
sx q[1];
rz(-1.2291127) q[1];
sx q[1];
rz(2.0613964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0475104) q[0];
sx q[0];
rz(-1.2258397) q[0];
sx q[0];
rz(-0.10070078) q[0];
x q[1];
rz(2.2713234) q[2];
sx q[2];
rz(-0.53127161) q[2];
sx q[2];
rz(-2.9722555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2094592) q[1];
sx q[1];
rz(-1.5725296) q[1];
sx q[1];
rz(-2.9936547) q[1];
rz(-pi) q[2];
rz(2.4467434) q[3];
sx q[3];
rz(-1.3245305) q[3];
sx q[3];
rz(-1.6978557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(-0.62823137) q[2];
rz(0.96588165) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-0.32111827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1107165) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(-1.734717) q[0];
rz(3.0137317) q[1];
sx q[1];
rz(-1.4297994) q[1];
sx q[1];
rz(-2.0157287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073407452) q[0];
sx q[0];
rz(-2.3940384) q[0];
sx q[0];
rz(-2.8697467) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01861219) q[2];
sx q[2];
rz(-1.1191223) q[2];
sx q[2];
rz(1.9960595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8496823) q[1];
sx q[1];
rz(-2.9314802) q[1];
sx q[1];
rz(0.30847524) q[1];
x q[2];
rz(0.58784318) q[3];
sx q[3];
rz(-0.70826642) q[3];
sx q[3];
rz(-1.4204587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0109978) q[2];
sx q[2];
rz(-2.0086292) q[2];
sx q[2];
rz(-0.093718378) q[2];
rz(-1.7081918) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(1.3933498) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.691064) q[0];
sx q[0];
rz(-1.0858303) q[0];
sx q[0];
rz(1.0040671) q[0];
rz(-1.1035236) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(-1.4080661) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133014) q[0];
sx q[0];
rz(-2.6264418) q[0];
sx q[0];
rz(-3.1293082) q[0];
x q[1];
rz(1.682071) q[2];
sx q[2];
rz(-1.2603052) q[2];
sx q[2];
rz(1.7686012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4346659) q[1];
sx q[1];
rz(-1.7499171) q[1];
sx q[1];
rz(2.0113025) q[1];
rz(-pi) q[2];
rz(1.5212112) q[3];
sx q[3];
rz(-1.0523497) q[3];
sx q[3];
rz(0.29274669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39190009) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(0.54769546) q[2];
rz(2.8294166) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(-1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373229) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(2.3626589) q[0];
rz(-2.7670822) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(0.83190727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0078985) q[0];
sx q[0];
rz(-1.8323032) q[0];
sx q[0];
rz(2.4228839) q[0];
rz(-pi) q[1];
rz(-2.4291762) q[2];
sx q[2];
rz(-2.2210026) q[2];
sx q[2];
rz(-0.67931108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9495957) q[1];
sx q[1];
rz(-1.0284541) q[1];
sx q[1];
rz(-2.3344912) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5396773) q[3];
sx q[3];
rz(-2.3672769) q[3];
sx q[3];
rz(1.7208791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(-3.1331151) q[2];
rz(-0.40555412) q[3];
sx q[3];
rz(-1.4529934) q[3];
sx q[3];
rz(1.6976374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(2.0472732) q[1];
sx q[1];
rz(-0.37847395) q[1];
sx q[1];
rz(2.7056221) q[1];
rz(-3.0158184) q[2];
sx q[2];
rz(-1.3604506) q[2];
sx q[2];
rz(2.4764555) q[2];
rz(3.0892298) q[3];
sx q[3];
rz(-3.047245) q[3];
sx q[3];
rz(0.17518763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
