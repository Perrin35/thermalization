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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286228) q[0];
sx q[0];
rz(-1.6464454) q[0];
sx q[0];
rz(-1.4609758) q[0];
rz(-pi) q[1];
rz(-0.50067164) q[2];
sx q[2];
rz(-1.8953875) q[2];
sx q[2];
rz(0.90510923) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6554101) q[1];
sx q[1];
rz(-1.3453801) q[1];
sx q[1];
rz(0.42014007) q[1];
rz(-2.1610307) q[3];
sx q[3];
rz(-0.96348982) q[3];
sx q[3];
rz(0.052562873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0483094) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(1.3422356) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(-1.5526519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(2.192705) q[0];
rz(0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(0.92996517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0895871) q[0];
sx q[0];
rz(-1.4692592) q[0];
sx q[0];
rz(0.2137645) q[0];
x q[1];
rz(1.5783159) q[2];
sx q[2];
rz(-1.0656271) q[2];
sx q[2];
rz(-1.1669605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9252214) q[1];
sx q[1];
rz(-0.19268806) q[1];
sx q[1];
rz(-2.5641901) q[1];
rz(-pi) q[2];
rz(-0.61701507) q[3];
sx q[3];
rz(-0.91147826) q[3];
sx q[3];
rz(0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13964222) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(2.3993717) q[2];
rz(0.46418515) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(1.5531042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4699698) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(1.6999014) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(-2.2742719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21737032) q[0];
sx q[0];
rz(-1.5708367) q[0];
sx q[0];
rz(1.5702308) q[0];
x q[1];
rz(0.27133743) q[2];
sx q[2];
rz(-1.0808766) q[2];
sx q[2];
rz(2.4957531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2056624) q[1];
sx q[1];
rz(-1.2453658) q[1];
sx q[1];
rz(1.8673351) q[1];
rz(-0.54384772) q[3];
sx q[3];
rz(-2.2529054) q[3];
sx q[3];
rz(-2.2206375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9637588) q[2];
sx q[2];
rz(-1.1744262) q[2];
sx q[2];
rz(-0.7589232) q[2];
rz(-0.80398503) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3114965) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-2.6575644) q[0];
rz(-1.6054035) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(1.4253634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4771627) q[0];
sx q[0];
rz(-1.2331404) q[0];
sx q[0];
rz(-2.5974524) q[0];
rz(-2.5192833) q[2];
sx q[2];
rz(-2.7633585) q[2];
sx q[2];
rz(2.6985199) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63923525) q[1];
sx q[1];
rz(-2.4662152) q[1];
sx q[1];
rz(1.2137854) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4590153) q[3];
sx q[3];
rz(-2.5831476) q[3];
sx q[3];
rz(-0.48170939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75590762) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(-2.8483086) q[2];
rz(0.048132345) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(2.8369246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3047979) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(1.1037214) q[0];
rz(0.91148218) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(-0.10890659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24507228) q[0];
sx q[0];
rz(-0.86045107) q[0];
sx q[0];
rz(2.4857387) q[0];
x q[1];
rz(-2.6066066) q[2];
sx q[2];
rz(-1.4145383) q[2];
sx q[2];
rz(-1.6339169) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8629294) q[1];
sx q[1];
rz(-0.68611523) q[1];
sx q[1];
rz(1.3115694) q[1];
rz(-pi) q[2];
rz(-0.44556983) q[3];
sx q[3];
rz(-0.62544367) q[3];
sx q[3];
rz(-1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4917422) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(2.7598925) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3968762) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(1.1599524) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(-2.8181308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8628629) q[0];
sx q[0];
rz(-0.6375618) q[0];
sx q[0];
rz(2.3414073) q[0];
x q[1];
rz(1.3176962) q[2];
sx q[2];
rz(-2.0539114) q[2];
sx q[2];
rz(0.3699257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0486794) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(-2.5552555) q[1];
rz(-0.46053912) q[3];
sx q[3];
rz(-1.8611188) q[3];
sx q[3];
rz(-0.059062171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20299992) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(0.8482376) q[2];
rz(-1.2674468) q[3];
sx q[3];
rz(-1.4013314) q[3];
sx q[3];
rz(0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11442014) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(-0.87345901) q[0];
rz(1.2365485) q[1];
sx q[1];
rz(-0.97573391) q[1];
sx q[1];
rz(3.0580318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75514166) q[0];
sx q[0];
rz(-1.6928732) q[0];
sx q[0];
rz(-0.32213078) q[0];
rz(-pi) q[1];
rz(0.23891831) q[2];
sx q[2];
rz(-2.0878125) q[2];
sx q[2];
rz(2.1674402) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8712743) q[1];
sx q[1];
rz(-1.9650998) q[1];
sx q[1];
rz(-2.32347) q[1];
rz(-0.6912937) q[3];
sx q[3];
rz(-0.45540998) q[3];
sx q[3];
rz(-2.5728855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43625912) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(-0.27240917) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.2347539) q[0];
sx q[0];
rz(-2.315157) q[0];
sx q[0];
rz(-2.2032264) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-0.96790853) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7308189) q[0];
sx q[0];
rz(-1.2026498) q[0];
sx q[0];
rz(-2.9981722) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16751473) q[2];
sx q[2];
rz(-1.9394839) q[2];
sx q[2];
rz(2.9978254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31466111) q[1];
sx q[1];
rz(-2.1912327) q[1];
sx q[1];
rz(0.13038306) q[1];
rz(-0.16444178) q[3];
sx q[3];
rz(-1.9666082) q[3];
sx q[3];
rz(-0.37104097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9390823) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(-1.067777) q[2];
rz(-2.0104525) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-0.080032674) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.985567) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(-0.40618968) q[0];
rz(-1.73229) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-2.0565775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34124464) q[0];
sx q[0];
rz(-0.64697504) q[0];
sx q[0];
rz(0.2917618) q[0];
x q[1];
rz(1.5753393) q[2];
sx q[2];
rz(-2.1811003) q[2];
sx q[2];
rz(-1.4686327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0250562) q[1];
sx q[1];
rz(-2.3485314) q[1];
sx q[1];
rz(-1.1021986) q[1];
rz(-1.0228588) q[3];
sx q[3];
rz(-2.0798488) q[3];
sx q[3];
rz(-1.7804543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5857508) q[2];
sx q[2];
rz(-0.61183524) q[2];
sx q[2];
rz(0.93070585) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68468204) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(1.8344301) q[0];
rz(-0.3859418) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(1.0345667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86901122) q[0];
sx q[0];
rz(-2.9391461) q[0];
sx q[0];
rz(-2.1073209) q[0];
x q[1];
rz(-2.7427086) q[2];
sx q[2];
rz(-2.6475057) q[2];
sx q[2];
rz(1.8979771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99920995) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(0.56908619) q[1];
rz(-2.3307269) q[3];
sx q[3];
rz(-1.6603784) q[3];
sx q[3];
rz(2.5406264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6941541) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(-0.36699692) q[2];
rz(1.1191818) q[3];
sx q[3];
rz(-0.39803353) q[3];
sx q[3];
rz(-1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(3.0837334) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(-0.51390263) q[2];
sx q[2];
rz(-0.19795098) q[2];
sx q[2];
rz(1.2951938) q[2];
rz(2.7230311) q[3];
sx q[3];
rz(-0.7444612) q[3];
sx q[3];
rz(0.85519467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
