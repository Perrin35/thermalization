OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13745452) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(-0.11560346) q[0];
rz(1.9723188) q[2];
sx q[2];
rz(-1.2520773) q[2];
sx q[2];
rz(-0.67414325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9143608) q[1];
sx q[1];
rz(-1.1162317) q[1];
sx q[1];
rz(-2.582002) q[1];
rz(-pi) q[2];
rz(0.75833851) q[3];
sx q[3];
rz(-1.7539277) q[3];
sx q[3];
rz(-1.6299712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(2.3056727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056136925) q[0];
sx q[0];
rz(-1.5526999) q[0];
sx q[0];
rz(1.4319112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1199273) q[2];
sx q[2];
rz(-1.2163391) q[2];
sx q[2];
rz(-2.3159388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.086120124) q[1];
sx q[1];
rz(-1.2985843) q[1];
sx q[1];
rz(-0.91113669) q[1];
rz(-1.9637945) q[3];
sx q[3];
rz(-2.2713695) q[3];
sx q[3];
rz(-2.8515479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.345984) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99840435) q[0];
sx q[0];
rz(-2.8371053) q[0];
sx q[0];
rz(-1.2782774) q[0];
x q[1];
rz(1.1630467) q[2];
sx q[2];
rz(-1.3464658) q[2];
sx q[2];
rz(-0.11220223) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25698173) q[1];
sx q[1];
rz(-1.7092488) q[1];
sx q[1];
rz(2.715766) q[1];
rz(-pi) q[2];
rz(2.0655572) q[3];
sx q[3];
rz(-1.6238188) q[3];
sx q[3];
rz(0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6879038) q[0];
sx q[0];
rz(-2.9691302) q[0];
sx q[0];
rz(2.0989291) q[0];
x q[1];
rz(-1.8727826) q[2];
sx q[2];
rz(-2.3846845) q[2];
sx q[2];
rz(2.7243171) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7820218) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(-0.77918474) q[1];
x q[2];
rz(-2.5883834) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33870944) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(2.5849735) q[0];
rz(-3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(2.1690878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199729) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(0.58453154) q[0];
rz(0.82053484) q[2];
sx q[2];
rz(-1.8168601) q[2];
sx q[2];
rz(0.29124242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79341187) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-2.1972149) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95485165) q[3];
sx q[3];
rz(-1.6007489) q[3];
sx q[3];
rz(-1.1101462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(1.548432) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71972972) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(-2.0772207) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(-0.37429601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6573338) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(1.2249468) q[0];
rz(-pi) q[1];
rz(-2.6475545) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(1.5649232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1230871) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(-1.8803384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66773325) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(-2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7999254) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-3.0068586) q[0];
x q[1];
rz(-2.4285165) q[2];
sx q[2];
rz(-1.8250069) q[2];
sx q[2];
rz(-1.6377246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20565614) q[1];
sx q[1];
rz(-1.9034791) q[1];
sx q[1];
rz(-1.7722305) q[1];
x q[2];
rz(-3.0974) q[3];
sx q[3];
rz(-2.5857539) q[3];
sx q[3];
rz(-1.0122055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(1.3423052) q[0];
x q[1];
rz(-1.4555898) q[2];
sx q[2];
rz(-2.4552279) q[2];
sx q[2];
rz(2.3442868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97493193) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(1.5207661) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3109342) q[3];
sx q[3];
rz(-1.4514918) q[3];
sx q[3];
rz(-0.80727778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(2.3855709) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5091771) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.8639494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.02361) q[0];
sx q[0];
rz(-1.8007468) q[0];
sx q[0];
rz(0.50595565) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5759376) q[2];
sx q[2];
rz(-0.7753765) q[2];
sx q[2];
rz(-2.2765809) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2904418) q[1];
sx q[1];
rz(-2.6944707) q[1];
sx q[1];
rz(1.4852344) q[1];
x q[2];
rz(0.83417474) q[3];
sx q[3];
rz(-1.8636384) q[3];
sx q[3];
rz(1.0398231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8032802) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(-0.35783106) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(-3.0293368) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(-0.21044883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4608599) q[0];
sx q[0];
rz(-1.0675149) q[0];
sx q[0];
rz(-1.3637278) q[0];
rz(-1.9333282) q[2];
sx q[2];
rz(-0.4368442) q[2];
sx q[2];
rz(-2.3853962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12293832) q[1];
sx q[1];
rz(-1.8784874) q[1];
sx q[1];
rz(1.7880746) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8912192) q[3];
sx q[3];
rz(-1.8009406) q[3];
sx q[3];
rz(1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18008733) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(1.2735584) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(2.9755637) q[2];
sx q[2];
rz(-1.1062853) q[2];
sx q[2];
rz(1.7044978) q[2];
rz(-1.1917226) q[3];
sx q[3];
rz(-2.3754397) q[3];
sx q[3];
rz(0.55879186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
