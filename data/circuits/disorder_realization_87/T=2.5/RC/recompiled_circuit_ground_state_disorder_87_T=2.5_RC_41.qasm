OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3992231) q[0];
sx q[0];
rz(-2.7288781) q[0];
sx q[0];
rz(-0.8362008) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158186) q[0];
sx q[0];
rz(-1.649722) q[0];
sx q[0];
rz(1.945334) q[0];
x q[1];
rz(0.37990976) q[2];
sx q[2];
rz(-1.7899982) q[2];
sx q[2];
rz(0.70218147) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5656149) q[1];
sx q[1];
rz(-1.4139785) q[1];
sx q[1];
rz(-2.6601936) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3775695) q[3];
sx q[3];
rz(-1.591103) q[3];
sx q[3];
rz(-1.8373674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8334373) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(-3.0270882) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(-2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(1.9812745) q[1];
sx q[1];
rz(-2.7313045) q[1];
sx q[1];
rz(-2.5210099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68408332) q[0];
sx q[0];
rz(-2.7433222) q[0];
sx q[0];
rz(1.0578367) q[0];
rz(-pi) q[1];
rz(-2.0257904) q[2];
sx q[2];
rz(-0.7209076) q[2];
sx q[2];
rz(-1.4269331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1329886) q[1];
sx q[1];
rz(-2.3594366) q[1];
sx q[1];
rz(-0.52495606) q[1];
x q[2];
rz(2.1866131) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(1.3835088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6985942) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(-1.5551785) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(-2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8860633) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(1.080876) q[0];
rz(2.4983662) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(-0.08531514) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236508) q[0];
sx q[0];
rz(-2.1408014) q[0];
sx q[0];
rz(2.355769) q[0];
rz(-2.8181452) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(-1.8307387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6416723) q[1];
sx q[1];
rz(-2.9292078) q[1];
sx q[1];
rz(-2.4714064) q[1];
rz(-3.078104) q[3];
sx q[3];
rz(-2.5194019) q[3];
sx q[3];
rz(0.49294642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9366511) q[2];
sx q[2];
rz(-0.0095657883) q[2];
sx q[2];
rz(2.9452475) q[2];
rz(0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333756) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-2.7561482) q[0];
rz(1.4133981) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(0.24994303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7628544) q[0];
sx q[0];
rz(-0.97134274) q[0];
sx q[0];
rz(-2.6361639) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.087935178) q[2];
sx q[2];
rz(-1.6307978) q[2];
sx q[2];
rz(-1.1052002) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0562607) q[1];
sx q[1];
rz(-1.0823432) q[1];
sx q[1];
rz(1.7621098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3479191) q[3];
sx q[3];
rz(-1.4303615) q[3];
sx q[3];
rz(0.92257231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6803153) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(1.8810898) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(-0.48511037) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9125875) q[0];
sx q[0];
rz(-1.8857755) q[0];
sx q[0];
rz(-2.7664685) q[0];
rz(-2.2274687) q[2];
sx q[2];
rz(-2.6466922) q[2];
sx q[2];
rz(-0.47593853) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5647392) q[1];
sx q[1];
rz(-2.4702063) q[1];
sx q[1];
rz(-1.4935843) q[1];
rz(-pi) q[2];
rz(1.5077293) q[3];
sx q[3];
rz(-1.1707889) q[3];
sx q[3];
rz(-2.0162752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89852077) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(2.2122993) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.9418075) q[0];
sx q[0];
rz(-2.6564044) q[0];
sx q[0];
rz(-2.3727681) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(-0.88968712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2425491) q[0];
sx q[0];
rz(-0.43936037) q[0];
sx q[0];
rz(2.0813294) q[0];
rz(-pi) q[1];
rz(3.1402428) q[2];
sx q[2];
rz(-1.6151516) q[2];
sx q[2];
rz(2.7173017) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1353768) q[1];
sx q[1];
rz(-1.6172503) q[1];
sx q[1];
rz(1.7170639) q[1];
rz(-1.6857288) q[3];
sx q[3];
rz(-2.0191666) q[3];
sx q[3];
rz(-1.2796206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5865405) q[2];
rz(1.8709315) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(0.13203013) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(1.9975115) q[0];
rz(-2.3106958) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(-2.3544618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2427776) q[0];
sx q[0];
rz(-1.4393596) q[0];
sx q[0];
rz(-1.4351033) q[0];
rz(-pi) q[1];
rz(0.93573715) q[2];
sx q[2];
rz(-2.3496186) q[2];
sx q[2];
rz(2.9284649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5097376) q[1];
sx q[1];
rz(-0.69586772) q[1];
sx q[1];
rz(2.3515755) q[1];
rz(-pi) q[2];
rz(2.5821286) q[3];
sx q[3];
rz(-1.2351523) q[3];
sx q[3];
rz(2.1223202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-3.1084295) q[2];
rz(1.5432594) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(-1.6395114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.892266) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(1.7484885) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(1.1005864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8580061) q[0];
sx q[0];
rz(-1.6467983) q[0];
sx q[0];
rz(-1.572524) q[0];
rz(0.69533657) q[2];
sx q[2];
rz(-2.0954663) q[2];
sx q[2];
rz(-2.6161043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7707899) q[1];
sx q[1];
rz(-0.90401559) q[1];
sx q[1];
rz(0.043170269) q[1];
rz(2.750179) q[3];
sx q[3];
rz(-1.4415223) q[3];
sx q[3];
rz(1.5551274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0988203) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(-1.5607321) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(-1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(-1.2358707) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(-0.12289563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.204958) q[0];
sx q[0];
rz(-2.7911148) q[0];
sx q[0];
rz(2.9913783) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2332623) q[2];
sx q[2];
rz(-0.35229063) q[2];
sx q[2];
rz(-0.34991821) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5220241) q[1];
sx q[1];
rz(-1.777473) q[1];
sx q[1];
rz(-3.0289438) q[1];
x q[2];
rz(-1.6121665) q[3];
sx q[3];
rz(-2.9968516) q[3];
sx q[3];
rz(2.5998757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5128532) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(-2.0780308) q[2];
rz(2.9641446) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43750957) q[0];
sx q[0];
rz(-0.45314416) q[0];
sx q[0];
rz(-1.7105239) q[0];
rz(2.5408632) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-0.61753714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73752943) q[0];
sx q[0];
rz(-1.6949708) q[0];
sx q[0];
rz(-1.7880718) q[0];
x q[1];
rz(2.3138758) q[2];
sx q[2];
rz(-1.1679782) q[2];
sx q[2];
rz(0.043443505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75467089) q[1];
sx q[1];
rz(-2.6082509) q[1];
sx q[1];
rz(-1.6120595) q[1];
x q[2];
rz(2.6838949) q[3];
sx q[3];
rz(-2.5182596) q[3];
sx q[3];
rz(-1.3184354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2403229) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(2.1642302) q[2];
rz(-1.0824925) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8025773) q[0];
sx q[0];
rz(-0.76887283) q[0];
sx q[0];
rz(-0.73706891) q[0];
rz(3.0958685) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(0.81188079) q[2];
sx q[2];
rz(-2.4627081) q[2];
sx q[2];
rz(-2.2344786) q[2];
rz(1.4010728) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
