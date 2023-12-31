OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.035064) q[0];
sx q[0];
rz(-2.0523235) q[0];
sx q[0];
rz(2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0954477) q[0];
sx q[0];
rz(-2.5382222) q[0];
sx q[0];
rz(1.7741816) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3688989) q[2];
sx q[2];
rz(-1.2374094) q[2];
sx q[2];
rz(-0.29793973) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6160994) q[1];
sx q[1];
rz(-1.6615189) q[1];
sx q[1];
rz(2.4331122) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7339488) q[3];
sx q[3];
rz(-1.4789494) q[3];
sx q[3];
rz(2.7013005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.5637406) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6710885) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(-1.6862306) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73194506) q[2];
sx q[2];
rz(-2.7436896) q[2];
sx q[2];
rz(-1.0674455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2240552) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(-0.097150306) q[1];
x q[2];
rz(1.7277667) q[3];
sx q[3];
rz(-1.4179215) q[3];
sx q[3];
rz(-0.65519858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-2.0593026) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(-2.6813566) q[0];
rz(-pi) q[1];
rz(0.41855721) q[2];
sx q[2];
rz(-1.3658128) q[2];
sx q[2];
rz(0.95450729) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5971165) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(2.4878923) q[1];
rz(1.5738437) q[3];
sx q[3];
rz(-1.8584195) q[3];
sx q[3];
rz(-2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-0.90467492) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(2.604565) q[0];
rz(2.5771192) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(-2.0828473) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0996453) q[1];
sx q[1];
rz(-1.6640267) q[1];
sx q[1];
rz(-2.9779742) q[1];
rz(-pi) q[2];
rz(-2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.4286208) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(0.3516745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85580101) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-2.48824) q[0];
rz(-pi) q[1];
rz(2.2064477) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(0.93174975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56375757) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(-1.5549591) q[1];
x q[2];
rz(-1.0293803) q[3];
sx q[3];
rz(-0.57313985) q[3];
sx q[3];
rz(2.074632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11739843) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-1.0995964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0352286) q[0];
sx q[0];
rz(-2.0237676) q[0];
sx q[0];
rz(0.41082541) q[0];
rz(-pi) q[1];
rz(1.0397644) q[2];
sx q[2];
rz(-1.4073939) q[2];
sx q[2];
rz(0.49887564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6096566) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(-0.68105662) q[1];
rz(0.78824708) q[3];
sx q[3];
rz(-1.8133834) q[3];
sx q[3];
rz(-1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-0.15730102) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29042127) q[0];
sx q[0];
rz(-2.3010845) q[0];
sx q[0];
rz(1.69676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.332875) q[2];
sx q[2];
rz(-0.750713) q[2];
sx q[2];
rz(0.44516341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.09307043) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(-2.1919769) q[1];
rz(-2.5141684) q[3];
sx q[3];
rz(-1.809123) q[3];
sx q[3];
rz(1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(-3.1307401) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5787443) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(-1.5638652) q[0];
rz(-pi) q[1];
rz(1.3515527) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(-3.0041681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7640904) q[1];
sx q[1];
rz(-1.9891771) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(0.71473748) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(1.294567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-1.0639169) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.3508505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0161184) q[0];
sx q[0];
rz(-3.0258412) q[0];
sx q[0];
rz(0.8158169) q[0];
rz(-pi) q[1];
rz(1.6532134) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(1.8694307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1119712) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(3.1236468) q[1];
rz(-pi) q[2];
rz(2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(2.9719877) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039625) q[0];
sx q[0];
rz(-0.61015218) q[0];
sx q[0];
rz(2.2274341) q[0];
rz(-1.595396) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(1.8116902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1093724) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(2.1857775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(-2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-2.8725273) q[2];
sx q[2];
rz(-1.8620327) q[2];
sx q[2];
rz(-3.0897683) q[2];
rz(1.7789755) q[3];
sx q[3];
rz(-1.6856185) q[3];
sx q[3];
rz(0.57193397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
