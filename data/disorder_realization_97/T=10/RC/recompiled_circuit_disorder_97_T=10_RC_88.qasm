OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4176183) q[0];
sx q[0];
rz(-1.4899878) q[0];
sx q[0];
rz(2.2111501) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(-2.0342483) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270277) q[0];
sx q[0];
rz(-1.2455997) q[0];
sx q[0];
rz(0.70744275) q[0];
rz(-pi) q[1];
rz(-1.0385752) q[2];
sx q[2];
rz(-2.3468809) q[2];
sx q[2];
rz(0.46193916) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5109374) q[1];
sx q[1];
rz(-1.2927379) q[1];
sx q[1];
rz(3.0250938) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34875617) q[3];
sx q[3];
rz(-2.1766799) q[3];
sx q[3];
rz(-2.6989258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1399122) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(-1.2225479) q[0];
rz(-pi) q[1];
rz(1.0802644) q[2];
sx q[2];
rz(-1.3776407) q[2];
sx q[2];
rz(-2.1527388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0751794) q[1];
sx q[1];
rz(-1.315409) q[1];
sx q[1];
rz(0.65402072) q[1];
rz(-pi) q[2];
rz(0.18873429) q[3];
sx q[3];
rz(-2.0368324) q[3];
sx q[3];
rz(2.1621494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(0.98637995) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36689511) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-0.96631518) q[0];
rz(-2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(1.6945217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149949) q[0];
sx q[0];
rz(-1.5883755) q[0];
sx q[0];
rz(2.9265762) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8410728) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(-1.2031872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9008873) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(-2.4918873) q[1];
x q[2];
rz(-1.0560016) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(-1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(-2.7526061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8303191) q[0];
sx q[0];
rz(-1.3415601) q[0];
sx q[0];
rz(-2.8515408) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.94611) q[2];
sx q[2];
rz(-1.0921548) q[2];
sx q[2];
rz(3.0340956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71967857) q[1];
sx q[1];
rz(-2.419157) q[1];
sx q[1];
rz(2.4721485) q[1];
rz(-pi) q[2];
rz(-0.74794482) q[3];
sx q[3];
rz(-1.4191295) q[3];
sx q[3];
rz(-2.9880854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(-2.6848865) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(2.8592498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(-2.7815681) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(-2.6470851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8202782) q[0];
sx q[0];
rz(-1.8501794) q[0];
sx q[0];
rz(-1.7000291) q[0];
rz(-pi) q[1];
rz(0.039426609) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(2.6089422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33489409) q[1];
sx q[1];
rz(-2.1999199) q[1];
sx q[1];
rz(1.6449528) q[1];
rz(-pi) q[2];
rz(1.5962283) q[3];
sx q[3];
rz(-1.126822) q[3];
sx q[3];
rz(2.3412995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(-0.29850706) q[2];
rz(2.8295637) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(1.235199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4069177) q[0];
sx q[0];
rz(-1.9039246) q[0];
sx q[0];
rz(-1.8963277) q[0];
rz(-2.4371106) q[2];
sx q[2];
rz(-1.3085758) q[2];
sx q[2];
rz(-0.30202497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27855733) q[1];
sx q[1];
rz(-0.99579358) q[1];
sx q[1];
rz(1.0689736) q[1];
rz(-pi) q[2];
rz(-2.8060693) q[3];
sx q[3];
rz(-0.24970679) q[3];
sx q[3];
rz(-1.5229131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(-0.43169272) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(-3.0294763) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(2.0281866) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413234) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(-1.8510438) q[0];
rz(-pi) q[1];
rz(0.31539519) q[2];
sx q[2];
rz(-1.8176259) q[2];
sx q[2];
rz(0.95755267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0943162) q[1];
sx q[1];
rz(-2.2502406) q[1];
sx q[1];
rz(-1.9654771) q[1];
rz(1.9951622) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(-0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814608) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(-1.8966282) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891114) q[0];
sx q[0];
rz(-1.4915823) q[0];
sx q[0];
rz(-1.3080025) q[0];
rz(-pi) q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-0.94617832) q[2];
sx q[2];
rz(2.5059932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3371256) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(2.0130403) q[1];
rz(0.21945159) q[3];
sx q[3];
rz(-0.95153522) q[3];
sx q[3];
rz(2.0451562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.7804902) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(2.7499054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93145934) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(-1.5643442) q[0];
x q[1];
rz(1.7708771) q[2];
sx q[2];
rz(-2.0152425) q[2];
sx q[2];
rz(1.0619628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.081025) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(-2.702436) q[1];
x q[2];
rz(-0.12556062) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(0.76198602) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-15/(14*pi)) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(-0.57089943) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(2.9796519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2019661) q[0];
sx q[0];
rz(-1.2830178) q[0];
sx q[0];
rz(2.1956325) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6185206) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(-1.7392841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0429749) q[1];
sx q[1];
rz(-0.21511714) q[1];
sx q[1];
rz(-2.2153562) q[1];
x q[2];
rz(-3.083769) q[3];
sx q[3];
rz(-0.52383808) q[3];
sx q[3];
rz(0.1334838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.4282248) q[2];
rz(-1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(-1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006871) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(-1.7715001) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(-0.8926819) q[2];
sx q[2];
rz(-0.18845367) q[2];
sx q[2];
rz(-1.0343196) q[2];
rz(-0.42967038) q[3];
sx q[3];
rz(-1.771314) q[3];
sx q[3];
rz(0.41021456) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
