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
rz(1.4752969) q[0];
sx q[0];
rz(-1.2694321) q[0];
sx q[0];
rz(-2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2494932) q[0];
sx q[0];
rz(-1.6803015) q[0];
sx q[0];
rz(3.0654869) q[0];
x q[1];
rz(2.640921) q[2];
sx q[2];
rz(-1.8953875) q[2];
sx q[2];
rz(-2.2364834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1266844) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(1.3247299) q[1];
rz(0.69656482) q[3];
sx q[3];
rz(-1.0961354) q[3];
sx q[3];
rz(-1.153095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(1.7993571) q[2];
rz(1.7736769) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.20464483) q[0];
sx q[0];
rz(-2.1277675) q[0];
sx q[0];
rz(-2.192705) q[0];
rz(-0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(-0.92996517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052005589) q[0];
sx q[0];
rz(-1.6723335) q[0];
sx q[0];
rz(2.9278281) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5783159) q[2];
sx q[2];
rz(-2.0759656) q[2];
sx q[2];
rz(-1.9746321) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2163712) q[1];
sx q[1];
rz(-2.9489046) q[1];
sx q[1];
rz(-0.57740258) q[1];
rz(-pi) q[2];
rz(-2.5245776) q[3];
sx q[3];
rz(-2.2301144) q[3];
sx q[3];
rz(0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13964222) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(0.742221) q[2];
rz(-2.6774075) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6716229) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.4416913) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(-0.86732078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7881667) q[0];
sx q[0];
rz(-1.5713619) q[0];
sx q[0];
rz(-4.0325469e-05) q[0];
rz(-pi) q[1];
rz(-0.27133743) q[2];
sx q[2];
rz(-1.0808766) q[2];
sx q[2];
rz(-2.4957531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2677541) q[1];
sx q[1];
rz(-1.8513362) q[1];
sx q[1];
rz(2.8024017) q[1];
rz(-pi) q[2];
rz(2.5977449) q[3];
sx q[3];
rz(-2.2529054) q[3];
sx q[3];
rz(0.92095514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9637588) q[2];
sx q[2];
rz(-1.1744262) q[2];
sx q[2];
rz(0.7589232) q[2];
rz(-0.80398503) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(-0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3114965) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(2.6575644) q[0];
rz(1.6054035) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(1.4253634) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40544505) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(0.59595705) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31230782) q[2];
sx q[2];
rz(-1.3538401) q[2];
sx q[2];
rz(1.4258476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63923525) q[1];
sx q[1];
rz(-2.4662152) q[1];
sx q[1];
rz(1.2137854) q[1];
rz(-pi) q[2];
rz(0.45141545) q[3];
sx q[3];
rz(-1.9115931) q[3];
sx q[3];
rz(1.6926852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(-0.29328406) q[2];
rz(0.048132345) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(2.8369246) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83679477) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(-2.0378713) q[0];
rz(-0.91148218) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(0.10890659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2809362) q[0];
sx q[0];
rz(-2.0514279) q[0];
sx q[0];
rz(-0.74451609) q[0];
rz(-pi) q[1];
rz(0.53498603) q[2];
sx q[2];
rz(-1.4145383) q[2];
sx q[2];
rz(-1.6339169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8629294) q[1];
sx q[1];
rz(-0.68611523) q[1];
sx q[1];
rz(1.3115694) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57753528) q[3];
sx q[3];
rz(-1.3157237) q[3];
sx q[3];
rz(2.5437742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4917422) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447164) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(-1.1599524) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(2.8181308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2787298) q[0];
sx q[0];
rz(-0.6375618) q[0];
sx q[0];
rz(2.3414073) q[0];
rz(-pi) q[1];
rz(0.49655621) q[2];
sx q[2];
rz(-1.3471896) q[2];
sx q[2];
rz(1.3204492) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66993078) q[1];
sx q[1];
rz(-2.1403229) q[1];
sx q[1];
rz(-1.8402151) q[1];
x q[2];
rz(-2.5497421) q[3];
sx q[3];
rz(-0.53880063) q[3];
sx q[3];
rz(-1.1065229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9385927) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(-0.8482376) q[2];
rz(-1.2674468) q[3];
sx q[3];
rz(-1.4013314) q[3];
sx q[3];
rz(-2.447824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(0.87345901) q[0];
rz(-1.2365485) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(-0.083560856) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1655052) q[0];
sx q[0];
rz(-0.3437316) q[0];
sx q[0];
rz(2.7718839) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9026743) q[2];
sx q[2];
rz(-1.0537801) q[2];
sx q[2];
rz(-0.97415249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2703184) q[1];
sx q[1];
rz(-1.9650998) q[1];
sx q[1];
rz(-0.81812268) q[1];
x q[2];
rz(1.8734345) q[3];
sx q[3];
rz(-1.2250966) q[3];
sx q[3];
rz(-1.8283591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7053335) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(1.7748888) q[2];
rz(2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90683872) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(-0.93836623) q[0];
rz(-2.3387108) q[1];
sx q[1];
rz(-0.96790853) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9296449) q[0];
sx q[0];
rz(-1.4370455) q[0];
sx q[0];
rz(-1.1991713) q[0];
x q[1];
rz(-1.9782045) q[2];
sx q[2];
rz(-0.40336868) q[2];
sx q[2];
rz(2.5591116) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8269315) q[1];
sx q[1];
rz(-0.95035997) q[1];
sx q[1];
rz(3.0112096) q[1];
x q[2];
rz(-1.9441767) q[3];
sx q[3];
rz(-0.42694091) q[3];
sx q[3];
rz(-0.77746848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20251033) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(-2.0738156) q[2];
rz(-2.0104525) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(0.40618968) q[0];
rz(1.4093026) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(1.0850151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34124464) q[0];
sx q[0];
rz(-2.4946176) q[0];
sx q[0];
rz(-2.8498309) q[0];
rz(-0.61030881) q[2];
sx q[2];
rz(-1.5745192) q[2];
sx q[2];
rz(-3.0368254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0250562) q[1];
sx q[1];
rz(-2.3485314) q[1];
sx q[1];
rz(-2.039394) q[1];
rz(-pi) q[2];
rz(2.5625251) q[3];
sx q[3];
rz(-2.0430312) q[3];
sx q[3];
rz(-3.0621665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5857508) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(1.3961004) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(-0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68468204) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(1.8344301) q[0];
rz(-2.7556509) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(2.1070259) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.912303) q[0];
sx q[0];
rz(-1.6737537) q[0];
sx q[0];
rz(1.7454171) q[0];
rz(-pi) q[1];
rz(-2.6808591) q[2];
sx q[2];
rz(-1.3855532) q[2];
sx q[2];
rz(-0.68250193) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9683796) q[1];
sx q[1];
rz(-2.3189841) q[1];
sx q[1];
rz(-0.93507018) q[1];
x q[2];
rz(-3.0183082) q[3];
sx q[3];
rz(-0.81466952) q[3];
sx q[3];
rz(2.2565763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6941541) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(-0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601396) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(3.0837334) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(-2.9686676) q[2];
sx q[2];
rz(-1.6676219) q[2];
sx q[2];
rz(0.22990083) q[2];
rz(0.41856159) q[3];
sx q[3];
rz(-2.3971315) q[3];
sx q[3];
rz(-2.286398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
