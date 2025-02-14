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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(-2.3469143) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2993752) q[0];
sx q[0];
rz(-1.2769852) q[0];
sx q[0];
rz(-0.037328193) q[0];
rz(1.8914521) q[2];
sx q[2];
rz(-0.76092623) q[2];
sx q[2];
rz(-2.2087532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2268035) q[1];
sx q[1];
rz(-2.7046596) q[1];
sx q[1];
rz(0.41586693) q[1];
rz(-pi) q[2];
rz(-2.6618979) q[3];
sx q[3];
rz(-0.51691662) q[3];
sx q[3];
rz(1.2933047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(-2.0726974) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(-1.6525035) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25817961) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(-2.6302443) q[0];
rz(-1.5072352) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(1.3828329) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6561474) q[0];
sx q[0];
rz(-1.1880179) q[0];
sx q[0];
rz(-2.9188927) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9963919) q[2];
sx q[2];
rz(-2.4239102) q[2];
sx q[2];
rz(-1.2956144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4310337) q[1];
sx q[1];
rz(-1.9703428) q[1];
sx q[1];
rz(0.33627681) q[1];
rz(-pi) q[2];
rz(2.6182943) q[3];
sx q[3];
rz(-1.4804258) q[3];
sx q[3];
rz(3.0449611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6766659) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(0.10291544) q[2];
rz(1.0627221) q[3];
sx q[3];
rz(-1.9461742) q[3];
sx q[3];
rz(3.0265871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677143) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(-3.0318731) q[0];
rz(-0.32490718) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(-0.5330162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4895297) q[0];
sx q[0];
rz(-2.8127413) q[0];
sx q[0];
rz(2.7016599) q[0];
rz(-0.43230482) q[2];
sx q[2];
rz(-0.057261618) q[2];
sx q[2];
rz(0.5106411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.24025) q[1];
sx q[1];
rz(-2.2247215) q[1];
sx q[1];
rz(-2.4010977) q[1];
rz(-pi) q[2];
rz(-1.1148861) q[3];
sx q[3];
rz(-0.81509903) q[3];
sx q[3];
rz(0.86642735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.816232) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(1.6858961) q[2];
rz(-3.0560737) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(2.1948309) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(-2.4192659) q[0];
rz(-1.5707312) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(3.0991203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52666506) q[0];
sx q[0];
rz(-0.4836429) q[0];
sx q[0];
rz(-2.7381104) q[0];
rz(-pi) q[1];
rz(3.1257079) q[2];
sx q[2];
rz(-1.4804891) q[2];
sx q[2];
rz(1.7866194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.153923) q[1];
sx q[1];
rz(-0.81036416) q[1];
sx q[1];
rz(2.8778879) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3298395) q[3];
sx q[3];
rz(-0.78380871) q[3];
sx q[3];
rz(-0.93894053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(-0.70146504) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0303845) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(0.17330387) q[0];
rz(2.2068842) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-0.25397837) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7807863) q[0];
sx q[0];
rz(-0.82366952) q[0];
sx q[0];
rz(-0.15374462) q[0];
rz(-pi) q[1];
rz(1.9333657) q[2];
sx q[2];
rz(-2.1340279) q[2];
sx q[2];
rz(-2.8993487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7813999) q[1];
sx q[1];
rz(-1.5500229) q[1];
sx q[1];
rz(-0.69333496) q[1];
rz(-1.1831102) q[3];
sx q[3];
rz(-1.5106987) q[3];
sx q[3];
rz(-2.5340486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0737334) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(-1.4402639) q[2];
rz(-0.67617792) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(0.079782709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9192231) q[0];
sx q[0];
rz(-1.8310522) q[0];
sx q[0];
rz(-2.4134912) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(-2.5808835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9084204) q[0];
sx q[0];
rz(-1.5058555) q[0];
sx q[0];
rz(-0.12663314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9496983) q[2];
sx q[2];
rz(-0.55804306) q[2];
sx q[2];
rz(3.0241338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.079472629) q[1];
sx q[1];
rz(-1.9626856) q[1];
sx q[1];
rz(-0.23418871) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0440935) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(1.9815552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5460983) q[2];
sx q[2];
rz(-1.6807154) q[2];
sx q[2];
rz(2.0687436) q[2];
rz(0.32235518) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(0.391092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3748465) q[0];
sx q[0];
rz(-1.5006737) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(0.32632581) q[1];
sx q[1];
rz(-0.34714454) q[1];
sx q[1];
rz(-1.6654642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107619) q[0];
sx q[0];
rz(-2.9617887) q[0];
sx q[0];
rz(-1.3489978) q[0];
x q[1];
rz(-2.9020643) q[2];
sx q[2];
rz(-1.0205185) q[2];
sx q[2];
rz(0.36877647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0261111) q[1];
sx q[1];
rz(-3.012595) q[1];
sx q[1];
rz(-1.1327101) q[1];
rz(-pi) q[2];
rz(1.7404775) q[3];
sx q[3];
rz(-1.7998003) q[3];
sx q[3];
rz(-1.9226216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8866715) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(2.2754748) q[2];
rz(-0.73650375) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(0.88732639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0199468) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.8628927) q[0];
rz(-2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-2.6967646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747975) q[0];
sx q[0];
rz(-0.54682362) q[0];
sx q[0];
rz(1.8575791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1110531) q[2];
sx q[2];
rz(-0.90219775) q[2];
sx q[2];
rz(2.8313314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0792744) q[1];
sx q[1];
rz(-1.6898815) q[1];
sx q[1];
rz(-0.80715413) q[1];
rz(-1.3108389) q[3];
sx q[3];
rz(-1.9366067) q[3];
sx q[3];
rz(-1.8577675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7030846) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-2.6355696) q[2];
rz(2.1972726) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29189062) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(-3.131409) q[0];
rz(2.5306375) q[1];
sx q[1];
rz(-0.98038951) q[1];
sx q[1];
rz(-3.0593061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5389494) q[0];
sx q[0];
rz(-1.2017439) q[0];
sx q[0];
rz(-2.6239388) q[0];
x q[1];
rz(-2.250483) q[2];
sx q[2];
rz(-1.3902877) q[2];
sx q[2];
rz(-2.0497104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.391606) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(0.27952607) q[1];
x q[2];
rz(1.3547586) q[3];
sx q[3];
rz(-0.91751614) q[3];
sx q[3];
rz(-0.21524425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-2.4324379) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(-0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3191147) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(2.5122232) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.1782882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8346658) q[0];
sx q[0];
rz(-2.3884058) q[0];
sx q[0];
rz(-2.1801394) q[0];
x q[1];
rz(2.7897116) q[2];
sx q[2];
rz(-2.0469672) q[2];
sx q[2];
rz(-1.9989524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3523184) q[1];
sx q[1];
rz(-1.6802009) q[1];
sx q[1];
rz(-2.3477702) q[1];
x q[2];
rz(2.516516) q[3];
sx q[3];
rz(-0.82603329) q[3];
sx q[3];
rz(0.99398289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(-0.81864041) q[2];
rz(1.5107752) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(-0.40217933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(1.2661487) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(-0.25898375) q[2];
sx q[2];
rz(-0.91194802) q[2];
sx q[2];
rz(2.9695445) q[2];
rz(1.6615909) q[3];
sx q[3];
rz(-1.3868854) q[3];
sx q[3];
rz(-1.8793061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
