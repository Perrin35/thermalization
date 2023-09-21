OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015895122) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(0.35981052) q[0];
x q[1];
rz(-1.4719226) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82422772) q[1];
sx q[1];
rz(-2.1050276) q[1];
sx q[1];
rz(-2.6298916) q[1];
rz(0.37681864) q[3];
sx q[3];
rz(-1.0816649) q[3];
sx q[3];
rz(-0.34987846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.6072134) q[2];
rz(-2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(-0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.672357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(-0.83067466) q[0];
rz(0.13164481) q[2];
sx q[2];
rz(-0.83877124) q[2];
sx q[2];
rz(-1.6694348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.037109) q[1];
sx q[1];
rz(-2.6625405) q[1];
sx q[1];
rz(0.80161174) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1727337) q[3];
sx q[3];
rz(-1.549984) q[3];
sx q[3];
rz(1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(-3.1243096) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-0.52454138) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0073111) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(-1.4071464) q[0];
x q[1];
rz(1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(-1.7922572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9279889) q[1];
sx q[1];
rz(-3.0756746) q[1];
sx q[1];
rz(-0.35857486) q[1];
rz(-pi) q[2];
rz(-1.1509622) q[3];
sx q[3];
rz(-1.7063147) q[3];
sx q[3];
rz(2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-3.0526429) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(0.68960062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423858) q[0];
sx q[0];
rz(-2.8028574) q[0];
sx q[0];
rz(2.6329106) q[0];
rz(0.34822779) q[2];
sx q[2];
rz(-0.59154445) q[2];
sx q[2];
rz(2.0915742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1414813) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(0.4517171) q[1];
rz(-0.48456405) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(2.518667) q[2];
rz(-2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899807) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(-1.6437644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9153584) q[0];
sx q[0];
rz(-0.81904531) q[0];
sx q[0];
rz(0.35924964) q[0];
x q[1];
rz(3.0450902) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(-0.9347136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.076482) q[1];
sx q[1];
rz(-1.5434192) q[1];
sx q[1];
rz(2.4379424) q[1];
rz(-pi) q[2];
rz(-2.6260914) q[3];
sx q[3];
rz(-1.7121797) q[3];
sx q[3];
rz(0.055671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(0.90562138) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5934138) q[0];
sx q[0];
rz(-0.81028623) q[0];
sx q[0];
rz(-0.20472783) q[0];
rz(-0.058850364) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(2.5318052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.012949) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(-3.0055771) q[1];
rz(-1.3330323) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(0.8386855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0662213) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-0.92648363) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28618318) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(-2.1320027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44430915) q[0];
sx q[0];
rz(-0.59519207) q[0];
sx q[0];
rz(2.5213581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6206714) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(-1.6233363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91612591) q[1];
sx q[1];
rz(-1.2327317) q[1];
sx q[1];
rz(-2.613693) q[1];
x q[2];
rz(-2.1647251) q[3];
sx q[3];
rz(-2.2403324) q[3];
sx q[3];
rz(-2.9059448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0895773) q[0];
sx q[0];
rz(-1.1451632) q[0];
sx q[0];
rz(-0.94469597) q[0];
rz(0.66419454) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(-1.9265837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6726482) q[1];
sx q[1];
rz(-0.69744195) q[1];
sx q[1];
rz(1.3717321) q[1];
x q[2];
rz(0.41110699) q[3];
sx q[3];
rz(-1.4484222) q[3];
sx q[3];
rz(1.3947595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42177054) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(2.3262809) q[2];
rz(-1.0347962) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(2.8544193) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-2.8093991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824326) q[0];
sx q[0];
rz(-2.3258665) q[0];
sx q[0];
rz(-2.0072323) q[0];
x q[1];
rz(3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(-2.7547714) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13739535) q[1];
sx q[1];
rz(-1.157837) q[1];
sx q[1];
rz(1.8505627) q[1];
rz(-pi) q[2];
rz(-0.62854564) q[3];
sx q[3];
rz(-0.62918951) q[3];
sx q[3];
rz(-0.43720804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(-0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(2.8870781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33728889) q[0];
sx q[0];
rz(-2.8992607) q[0];
sx q[0];
rz(1.8810349) q[0];
rz(-0.48760957) q[2];
sx q[2];
rz(-1.3541823) q[2];
sx q[2];
rz(2.285514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75603112) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(0.24753333) q[1];
x q[2];
rz(-1.47154) q[3];
sx q[3];
rz(-1.1343065) q[3];
sx q[3];
rz(1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(2.9011177) q[2];
sx q[2];
rz(-1.0576116) q[2];
sx q[2];
rz(-0.53838421) q[2];
rz(-1.249282) q[3];
sx q[3];
rz(-2.5026863) q[3];
sx q[3];
rz(1.7659059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];