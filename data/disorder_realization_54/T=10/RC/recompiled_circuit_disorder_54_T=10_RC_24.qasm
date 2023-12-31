OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(-2.4523003) q[0];
sx q[0];
rz(0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256975) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(-0.35981052) q[0];
rz(-pi) q[1];
rz(-1.2572631) q[2];
sx q[2];
rz(-1.5402113) q[2];
sx q[2];
rz(0.41536301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.658537) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(-2.2621821) q[1];
rz(-pi) q[2];
rz(-2.1756644) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(-1.0498429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7222897) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(-0.62227917) q[0];
rz(-2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-2.2252749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4692357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(2.310918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4257405) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(-1.2765826) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20827046) q[1];
sx q[1];
rz(-1.9083438) q[1];
sx q[1];
rz(-0.34668215) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.968859) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24094412) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(-2.1333372) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-2.6170513) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265822) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.7517356) q[0];
rz(-pi) q[1];
rz(2.8183297) q[2];
sx q[2];
rz(-2.4369536) q[2];
sx q[2];
rz(-1.3779674) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4265392) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(-3.0798562) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14822652) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(-1.1566597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(0.51070172) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(0.68960062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(-0.50868209) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56324048) q[2];
sx q[2];
rz(-13/(3*pi)) q[2];
sx q[2];
rz(-2.9134977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2526306) q[1];
sx q[1];
rz(-1.2446212) q[1];
sx q[1];
rz(-2.3701282) q[1];
rz(0.48456405) q[3];
sx q[3];
rz(-1.6996517) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(-2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899807) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(0.583453) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.6437644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0935055) q[0];
sx q[0];
rz(-1.3110647) q[0];
sx q[0];
rz(0.78608677) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3439212) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(-2.5728512) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46206611) q[1];
sx q[1];
rz(-2.4375009) q[1];
sx q[1];
rz(0.042298869) q[1];
rz(2.8605117) q[3];
sx q[3];
rz(-2.6087458) q[3];
sx q[3];
rz(1.7588774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3008266) q[0];
sx q[0];
rz(-2.3593785) q[0];
sx q[0];
rz(1.3602815) q[0];
x q[1];
rz(-1.0211421) q[2];
sx q[2];
rz(-1.6209941) q[2];
sx q[2];
rz(-2.2113138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.012949) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(-3.0055771) q[1];
rz(-pi) q[2];
x q[2];
rz(0.027408882) q[3];
sx q[3];
rz(-1.3331183) q[3];
sx q[3];
rz(2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0662213) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(-1.4525157) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-0.92648363) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28618318) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(2.563971) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(2.1320027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6972835) q[0];
sx q[0];
rz(-0.59519207) q[0];
sx q[0];
rz(-2.5213581) q[0];
rz(1.7989484) q[2];
sx q[2];
rz(-1.9462799) q[2];
sx q[2];
rz(1.0588888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2254667) q[1];
sx q[1];
rz(-1.9088609) q[1];
sx q[1];
rz(-2.613693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97686751) q[3];
sx q[3];
rz(-2.2403324) q[3];
sx q[3];
rz(-2.9059448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-0.2391267) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(2.696864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3702104) q[0];
sx q[0];
rz(-2.133773) q[0];
sx q[0];
rz(-2.6315106) q[0];
rz(-2.4773981) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(1.2150089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6726482) q[1];
sx q[1];
rz(-2.4441507) q[1];
sx q[1];
rz(-1.7698606) q[1];
rz(-pi) q[2];
rz(2.7304857) q[3];
sx q[3];
rz(-1.4484222) q[3];
sx q[3];
rz(1.7468332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-2.8544193) q[0];
rz(-2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(2.8093991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824326) q[0];
sx q[0];
rz(-2.3258665) q[0];
sx q[0];
rz(1.1343603) q[0];
rz(-0.8719445) q[2];
sx q[2];
rz(-1.5351864) q[2];
sx q[2];
rz(1.2138838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5933983) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(0.42773186) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6094749) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(-1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93280783) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(-2.3020111) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(0.65657842) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(0.97958952) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33728889) q[0];
sx q[0];
rz(-2.8992607) q[0];
sx q[0];
rz(-1.8810349) q[0];
rz(-1.326667) q[2];
sx q[2];
rz(-2.0460528) q[2];
sx q[2];
rz(2.313386) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6470681) q[1];
sx q[1];
rz(-1.2459323) q[1];
sx q[1];
rz(-1.4855794) q[1];
x q[2];
rz(-2.7032095) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(-0.59059483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-2.7813773) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(2.9011177) q[2];
sx q[2];
rz(-1.0576116) q[2];
sx q[2];
rz(-0.53838421) q[2];
rz(1.8923106) q[3];
sx q[3];
rz(-2.5026863) q[3];
sx q[3];
rz(1.7659059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
