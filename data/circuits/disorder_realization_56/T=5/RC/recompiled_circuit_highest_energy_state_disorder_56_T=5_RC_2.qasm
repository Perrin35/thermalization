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
rz(-1.6662958) q[0];
sx q[0];
rz(-1.8721606) q[0];
sx q[0];
rz(2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(-2.3840005) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8920994) q[0];
sx q[0];
rz(-1.6803015) q[0];
sx q[0];
rz(-3.0654869) q[0];
x q[1];
rz(1.2045317) q[2];
sx q[2];
rz(-2.0431113) q[2];
sx q[2];
rz(0.83845316) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4861826) q[1];
sx q[1];
rz(-1.7962126) q[1];
sx q[1];
rz(2.7214526) q[1];
rz(-pi) q[2];
rz(-2.4663062) q[3];
sx q[3];
rz(-2.3216341) q[3];
sx q[3];
rz(-0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0483094) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(1.3422356) q[2];
rz(1.3679158) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(-1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.20464483) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(0.94888765) q[0];
rz(0.10143796) q[1];
sx q[1];
rz(-2.0815492) q[1];
sx q[1];
rz(2.2116275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0819433) q[0];
sx q[0];
rz(-0.23632061) q[0];
sx q[0];
rz(0.44775072) q[0];
rz(0.013596046) q[2];
sx q[2];
rz(-0.50522035) q[2];
sx q[2];
rz(-1.9590953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21446642) q[1];
sx q[1];
rz(-7*pi/15) q[1];
sx q[1];
rz(-0.16204496) q[1];
rz(-pi) q[2];
rz(-2.5245776) q[3];
sx q[3];
rz(-0.91147826) q[3];
sx q[3];
rz(-0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0019504) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(0.742221) q[2];
rz(2.6774075) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5531042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716229) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(-1.4416913) q[0];
rz(-2.9761159) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(-2.2742719) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246068) q[0];
sx q[0];
rz(-3.1410257) q[0];
sx q[0];
rz(1.4996155) q[0];
rz(-pi) q[1];
rz(0.27133743) q[2];
sx q[2];
rz(-1.0808766) q[2];
sx q[2];
rz(2.4957531) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9359303) q[1];
sx q[1];
rz(-1.8962269) q[1];
sx q[1];
rz(-1.2742576) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3300784) q[3];
sx q[3];
rz(-1.9841188) q[3];
sx q[3];
rz(1.0140401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1778339) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(-2.3826694) q[2];
rz(0.80398503) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(2.4132531) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8300962) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-2.6575644) q[0];
rz(1.5361891) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(-1.7162292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7361476) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(2.5456356) q[0];
x q[1];
rz(-2.8292848) q[2];
sx q[2];
rz(-1.3538401) q[2];
sx q[2];
rz(1.4258476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0565383) q[1];
sx q[1];
rz(-0.94496545) q[1];
sx q[1];
rz(2.8686348) q[1];
x q[2];
rz(-1.9462162) q[3];
sx q[3];
rz(-1.1470801) q[3];
sx q[3];
rz(2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(2.8483086) q[2];
rz(-0.048132345) q[3];
sx q[3];
rz(-1.8354514) q[3];
sx q[3];
rz(-0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3047979) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(2.0378713) q[0];
rz(0.91148218) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(0.10890659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965204) q[0];
sx q[0];
rz(-0.86045107) q[0];
sx q[0];
rz(2.4857387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7519195) q[2];
sx q[2];
rz(-2.0985773) q[2];
sx q[2];
rz(0.15508791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.089757925) q[1];
sx q[1];
rz(-1.4076774) q[1];
sx q[1];
rz(0.9012797) q[1];
x q[2];
rz(-2.6960228) q[3];
sx q[3];
rz(-0.62544367) q[3];
sx q[3];
rz(1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4917422) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(-2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3968762) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(-1.9816403) q[0];
rz(-2.5634735) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(-0.32346183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8628629) q[0];
sx q[0];
rz(-0.6375618) q[0];
sx q[0];
rz(0.80018534) q[0];
rz(-pi) q[1];
rz(0.44536369) q[2];
sx q[2];
rz(-2.6008743) q[2];
sx q[2];
rz(-0.13810829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0486794) q[1];
sx q[1];
rz(-1.7968751) q[1];
sx q[1];
rz(-0.58633713) q[1];
rz(-pi) q[2];
rz(0.46053912) q[3];
sx q[3];
rz(-1.2804739) q[3];
sx q[3];
rz(3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9385927) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(0.8482376) q[2];
rz(-1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11442014) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(0.87345901) q[0];
rz(-1.2365485) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(3.0580318) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.386451) q[0];
sx q[0];
rz(-1.4487195) q[0];
sx q[0];
rz(2.8194619) q[0];
x q[1];
rz(1.9651863) q[2];
sx q[2];
rz(-0.56496921) q[2];
sx q[2];
rz(2.625287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68971953) q[1];
sx q[1];
rz(-0.83137935) q[1];
sx q[1];
rz(-2.1175794) q[1];
x q[2];
rz(-0.6912937) q[3];
sx q[3];
rz(-2.6861827) q[3];
sx q[3];
rz(-0.56870715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43625912) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(-1.3667038) q[2];
rz(0.27240917) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(-0.75850707) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90683872) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(2.2032264) q[0];
rz(-2.3387108) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(2.3209007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4107738) q[0];
sx q[0];
rz(-1.9389429) q[0];
sx q[0];
rz(-0.14342043) q[0];
rz(1.9442476) q[2];
sx q[2];
rz(-1.7269584) q[2];
sx q[2];
rz(1.3661623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0487602) q[1];
sx q[1];
rz(-0.63221778) q[1];
sx q[1];
rz(1.7507751) q[1];
rz(2.9771509) q[3];
sx q[3];
rz(-1.1749845) q[3];
sx q[3];
rz(-2.7705517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9390823) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(1.067777) q[2];
rz(2.0104525) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.0850151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4403518) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(-1.7847654) q[0];
rz(-pi) q[1];
rz(-2.5312838) q[2];
sx q[2];
rz(-1.5745192) q[2];
sx q[2];
rz(-0.10476724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20445261) q[1];
sx q[1];
rz(-1.2431743) q[1];
sx q[1];
rz(0.83468584) q[1];
rz(1.0228588) q[3];
sx q[3];
rz(-1.0617439) q[3];
sx q[3];
rz(1.3611384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(0.93070585) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(-0.11955587) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4569106) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(-1.3071625) q[0];
rz(-2.7556509) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-2.1070259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912303) q[0];
sx q[0];
rz(-1.6737537) q[0];
sx q[0];
rz(1.3961755) q[0];
rz(-pi) q[1];
rz(-1.3645646) q[2];
sx q[2];
rz(-2.0230556) q[2];
sx q[2];
rz(-2.3444676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99920995) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(0.56908619) q[1];
rz(-pi) q[2];
rz(-1.7004556) q[3];
sx q[3];
rz(-2.377452) q[3];
sx q[3];
rz(-2.0779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(2.7745957) q[2];
rz(1.1191818) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.057859261) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(2.62769) q[2];
sx q[2];
rz(-0.19795098) q[2];
sx q[2];
rz(1.2951938) q[2];
rz(-0.41856159) q[3];
sx q[3];
rz(-0.7444612) q[3];
sx q[3];
rz(0.85519467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
