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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2494932) q[0];
sx q[0];
rz(-1.6803015) q[0];
sx q[0];
rz(3.0654869) q[0];
rz(-pi) q[1];
rz(2.5301928) q[2];
sx q[2];
rz(-2.5525103) q[2];
sx q[2];
rz(3.0036199) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54851531) q[1];
sx q[1];
rz(-2.6679949) q[1];
sx q[1];
rz(-0.51215902) q[1];
rz(-pi) q[2];
rz(-0.67528649) q[3];
sx q[3];
rz(-0.8199586) q[3];
sx q[3];
rz(2.2238071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0483094) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(1.7993571) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(1.5526519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20464483) q[0];
sx q[0];
rz(-2.1277675) q[0];
sx q[0];
rz(2.192705) q[0];
rz(3.0401547) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(-0.92996517) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052005589) q[0];
sx q[0];
rz(-1.6723335) q[0];
sx q[0];
rz(-2.9278281) q[0];
rz(-0.50518121) q[2];
sx q[2];
rz(-1.5773767) q[2];
sx q[2];
rz(2.7413961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3392433) q[1];
sx q[1];
rz(-1.4096469) q[1];
sx q[1];
rz(-1.4646962) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61701507) q[3];
sx q[3];
rz(-2.2301144) q[3];
sx q[3];
rz(-0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(0.742221) q[2];
rz(-0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5531042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6716229) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(-2.9761159) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(0.86732078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21737032) q[0];
sx q[0];
rz(-1.570756) q[0];
sx q[0];
rz(-1.5713619) q[0];
rz(2.0763511) q[2];
sx q[2];
rz(-1.3320413) q[2];
sx q[2];
rz(2.3468034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8738385) q[1];
sx q[1];
rz(-1.8513362) q[1];
sx q[1];
rz(0.33919097) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1380566) q[3];
sx q[3];
rz(-2.2972882) q[3];
sx q[3];
rz(-2.9851825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1778339) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(-2.3826694) q[2];
rz(2.3376076) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(2.4132531) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(0.48402825) q[0];
rz(-1.5361891) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(1.7162292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1041906) q[0];
sx q[0];
rz(-1.0604825) q[0];
sx q[0];
rz(1.9602106) q[0];
rz(-pi) q[1];
rz(-1.7984125) q[2];
sx q[2];
rz(-1.8755442) q[2];
sx q[2];
rz(-2.9272542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5023574) q[1];
sx q[1];
rz(-0.67537748) q[1];
sx q[1];
rz(-1.9278072) q[1];
x q[2];
rz(2.4590153) q[3];
sx q[3];
rz(-2.5831476) q[3];
sx q[3];
rz(0.48170939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(-2.8483086) q[2];
rz(3.0934603) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3047979) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(2.0378713) q[0];
rz(2.2301105) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(-0.10890659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24507228) q[0];
sx q[0];
rz(-2.2811416) q[0];
sx q[0];
rz(-0.65585391) q[0];
rz(-pi) q[1];
rz(1.3896732) q[2];
sx q[2];
rz(-1.0430153) q[2];
sx q[2];
rz(2.9865047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6088691) q[1];
sx q[1];
rz(-2.2298498) q[1];
sx q[1];
rz(-2.9347035) q[1];
rz(-pi) q[2];
rz(2.6960228) q[3];
sx q[3];
rz(-2.516149) q[3];
sx q[3];
rz(1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.3968762) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(1.1599524) q[0];
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
rz(-1.8628629) q[0];
sx q[0];
rz(-2.5040309) q[0];
sx q[0];
rz(-0.80018534) q[0];
rz(-0.44536369) q[2];
sx q[2];
rz(-2.6008743) q[2];
sx q[2];
rz(-3.0034844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9449294) q[1];
sx q[1];
rz(-2.5179836) q[1];
sx q[1];
rz(0.39396472) q[1];
x q[2];
rz(2.6810535) q[3];
sx q[3];
rz(-1.8611188) q[3];
sx q[3];
rz(-0.059062171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9385927) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(-2.2933551) q[2];
rz(1.8741459) q[3];
sx q[3];
rz(-1.4013314) q[3];
sx q[3];
rz(-2.447824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(-0.87345901) q[0];
rz(1.2365485) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75514166) q[0];
sx q[0];
rz(-1.6928732) q[0];
sx q[0];
rz(0.32213078) q[0];
rz(-pi) q[1];
rz(1.0413076) q[2];
sx q[2];
rz(-1.7779967) q[2];
sx q[2];
rz(0.71646128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4518731) q[1];
sx q[1];
rz(-0.83137935) q[1];
sx q[1];
rz(-2.1175794) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2681581) q[3];
sx q[3];
rz(-1.916496) q[3];
sx q[3];
rz(1.3132335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43625912) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(1.3667038) q[2];
rz(2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(-0.75850707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.90683872) q[0];
sx q[0];
rz(-2.315157) q[0];
sx q[0];
rz(0.93836623) q[0];
rz(2.3387108) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(0.82069194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4107738) q[0];
sx q[0];
rz(-1.2026498) q[0];
sx q[0];
rz(0.14342043) q[0];
rz(-pi) q[1];
rz(1.1633881) q[2];
sx q[2];
rz(-2.738224) q[2];
sx q[2];
rz(0.58248108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8269315) q[1];
sx q[1];
rz(-2.1912327) q[1];
sx q[1];
rz(0.13038306) q[1];
rz(-pi) q[2];
rz(1.971463) q[3];
sx q[3];
rz(-1.4191711) q[3];
sx q[3];
rz(-1.135863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9390823) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(-2.0738156) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-3.06156) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(-0.40618968) q[0];
rz(1.4093026) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(2.0565775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4403518) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(1.7847654) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61030881) q[2];
sx q[2];
rz(-1.5670735) q[2];
sx q[2];
rz(3.0368254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20445261) q[1];
sx q[1];
rz(-1.2431743) q[1];
sx q[1];
rz(-2.3069068) q[1];
x q[2];
rz(-0.75094838) q[3];
sx q[3];
rz(-0.72970684) q[3];
sx q[3];
rz(0.88353222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(-2.2108868) q[2];
rz(-1.3961004) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(-0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(-1.3071625) q[0];
rz(-2.7556509) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-2.1070259) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32337727) q[0];
sx q[0];
rz(-1.3971097) q[0];
sx q[0];
rz(-0.10453577) q[0];
rz(-pi) q[1];
rz(2.6808591) q[2];
sx q[2];
rz(-1.3855532) q[2];
sx q[2];
rz(-2.4590907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17321302) q[1];
sx q[1];
rz(-0.82260859) q[1];
sx q[1];
rz(-2.2065225) q[1];
rz(-pi) q[2];
rz(-1.7004556) q[3];
sx q[3];
rz(-0.76414062) q[3];
sx q[3];
rz(2.0779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(-0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781453) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-0.057859261) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(-1.6690785) q[2];
sx q[2];
rz(-1.7429033) q[2];
sx q[2];
rz(1.8175816) q[2];
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
