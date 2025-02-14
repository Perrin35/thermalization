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
rz(0.089224815) q[0];
sx q[0];
rz(2.3951946) q[0];
sx q[0];
rz(8.7328773) q[0];
rz(0.46299419) q[1];
sx q[1];
rz(-1.6166592) q[1];
sx q[1];
rz(2.1114299) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2041516) q[0];
sx q[0];
rz(-2.5290956) q[0];
sx q[0];
rz(0.4093616) q[0];
rz(-pi) q[1];
rz(1.9857731) q[2];
sx q[2];
rz(-2.7105467) q[2];
sx q[2];
rz(2.8243714) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2166442) q[1];
sx q[1];
rz(-2.7949667) q[1];
sx q[1];
rz(1.9199172) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9943732) q[3];
sx q[3];
rz(-0.69425827) q[3];
sx q[3];
rz(-2.8822299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2618711) q[2];
sx q[2];
rz(-2.3679831) q[2];
sx q[2];
rz(-0.36641463) q[2];
rz(2.053818) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.825603) q[0];
sx q[0];
rz(-0.082254224) q[0];
sx q[0];
rz(-2.2218623) q[0];
rz(-2.3986744) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(1.119335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3920721) q[0];
sx q[0];
rz(-1.6320092) q[0];
sx q[0];
rz(-1.4257512) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.934147) q[2];
sx q[2];
rz(-0.72978696) q[2];
sx q[2];
rz(-1.5130254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97174469) q[1];
sx q[1];
rz(-1.3339284) q[1];
sx q[1];
rz(2.4011517) q[1];
x q[2];
rz(0.85823409) q[3];
sx q[3];
rz(-2.4284114) q[3];
sx q[3];
rz(-0.76480097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4982345) q[2];
sx q[2];
rz(-1.9469807) q[2];
sx q[2];
rz(0.052113459) q[2];
rz(1.107996) q[3];
sx q[3];
rz(-2.2823157) q[3];
sx q[3];
rz(3.0767379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010083) q[0];
sx q[0];
rz(-1.2315741) q[0];
sx q[0];
rz(-1.8517866) q[0];
rz(2.3884804) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(0.49711102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.755794) q[0];
sx q[0];
rz(-2.442474) q[0];
sx q[0];
rz(-0.6363477) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2154664) q[2];
sx q[2];
rz(-2.1417649) q[2];
sx q[2];
rz(2.1059011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5376007) q[1];
sx q[1];
rz(-1.5328755) q[1];
sx q[1];
rz(-2.7673519) q[1];
rz(1.6670973) q[3];
sx q[3];
rz(-0.9203921) q[3];
sx q[3];
rz(1.4092829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46811732) q[2];
sx q[2];
rz(-0.27286369) q[2];
sx q[2];
rz(-2.4172778) q[2];
rz(2.916548) q[3];
sx q[3];
rz(-1.8658172) q[3];
sx q[3];
rz(-2.3465274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.1555136) q[0];
sx q[0];
rz(-2.2561095) q[0];
sx q[0];
rz(0.28050637) q[0];
rz(-1.4836503) q[1];
sx q[1];
rz(-1.9397441) q[1];
sx q[1];
rz(-0.4712421) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126033) q[0];
sx q[0];
rz(-1.5935197) q[0];
sx q[0];
rz(1.2034334) q[0];
rz(-pi) q[1];
rz(-1.4578826) q[2];
sx q[2];
rz(-1.6364003) q[2];
sx q[2];
rz(2.0173617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71472834) q[1];
sx q[1];
rz(-2.8694723) q[1];
sx q[1];
rz(1.6741899) q[1];
x q[2];
rz(-0.80914167) q[3];
sx q[3];
rz(-2.005072) q[3];
sx q[3];
rz(-2.6971522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19807854) q[2];
sx q[2];
rz(-0.55800617) q[2];
sx q[2];
rz(-2.1264709) q[2];
rz(-0.73457581) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(2.6034897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839142) q[0];
sx q[0];
rz(-1.2175125) q[0];
sx q[0];
rz(0.96018803) q[0];
rz(-3.0958815) q[1];
sx q[1];
rz(-0.8784596) q[1];
sx q[1];
rz(3.1083621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3094745) q[0];
sx q[0];
rz(-0.72369472) q[0];
sx q[0];
rz(-2.3242818) q[0];
rz(1.8383547) q[2];
sx q[2];
rz(-1.5853527) q[2];
sx q[2];
rz(-0.47786682) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0000764) q[1];
sx q[1];
rz(-1.4404529) q[1];
sx q[1];
rz(-1.9059859) q[1];
rz(-0.44132061) q[3];
sx q[3];
rz(-1.1989824) q[3];
sx q[3];
rz(-1.1692493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9783322) q[2];
sx q[2];
rz(-0.99908081) q[2];
sx q[2];
rz(-0.31326374) q[2];
rz(2.7836109) q[3];
sx q[3];
rz(-2.101818) q[3];
sx q[3];
rz(-3.0430005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98704308) q[0];
sx q[0];
rz(-1.510226) q[0];
sx q[0];
rz(0.48598591) q[0];
rz(-1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(-0.33014306) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81583038) q[0];
sx q[0];
rz(-1.5410402) q[0];
sx q[0];
rz(1.8682006) q[0];
x q[1];
rz(-2.2171872) q[2];
sx q[2];
rz(-2.2644448) q[2];
sx q[2];
rz(-2.2104757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79336004) q[1];
sx q[1];
rz(-0.72472445) q[1];
sx q[1];
rz(0.7591922) q[1];
x q[2];
rz(2.5810235) q[3];
sx q[3];
rz(-0.40321002) q[3];
sx q[3];
rz(1.3399269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6351472) q[2];
sx q[2];
rz(-2.5250285) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(0.59456524) q[3];
sx q[3];
rz(-1.6584572) q[3];
sx q[3];
rz(2.4026292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67721382) q[0];
sx q[0];
rz(-2.3742299) q[0];
sx q[0];
rz(-0.58694029) q[0];
rz(-0.37239536) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(0.24758235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42668396) q[0];
sx q[0];
rz(-1.6021906) q[0];
sx q[0];
rz(2.035729) q[0];
rz(-2.5093497) q[2];
sx q[2];
rz(-1.768501) q[2];
sx q[2];
rz(-1.3121669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5953427) q[1];
sx q[1];
rz(-2.3997325) q[1];
sx q[1];
rz(1.0651903) q[1];
rz(-1.8195175) q[3];
sx q[3];
rz(-1.910672) q[3];
sx q[3];
rz(-0.23409493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26329142) q[2];
sx q[2];
rz(-1.2211439) q[2];
sx q[2];
rz(2.9325874) q[2];
rz(0.59669295) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(1.5359115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7335032) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(-2.7485513) q[0];
rz(-2.4709002) q[1];
sx q[1];
rz(-1.1436983) q[1];
sx q[1];
rz(1.6938946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9056975) q[0];
sx q[0];
rz(-1.6950771) q[0];
sx q[0];
rz(-1.0727364) q[0];
rz(-1.1426439) q[2];
sx q[2];
rz(-0.65946666) q[2];
sx q[2];
rz(3.1127549) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5911342) q[1];
sx q[1];
rz(-1.0626284) q[1];
sx q[1];
rz(2.152446) q[1];
rz(-pi) q[2];
rz(-1.0974698) q[3];
sx q[3];
rz(-0.83126634) q[3];
sx q[3];
rz(-0.39287469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4286917) q[2];
sx q[2];
rz(-2.28076) q[2];
sx q[2];
rz(-3.0248896) q[2];
rz(-2.1278925) q[3];
sx q[3];
rz(-1.9096392) q[3];
sx q[3];
rz(-1.8359756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75334615) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(1.4353132) q[0];
rz(2.1460136) q[1];
sx q[1];
rz(-1.8228056) q[1];
sx q[1];
rz(-0.54042712) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0198959) q[0];
sx q[0];
rz(-1.833124) q[0];
sx q[0];
rz(-2.645825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2378906) q[2];
sx q[2];
rz(-2.0533178) q[2];
sx q[2];
rz(-0.36090349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3167773) q[1];
sx q[1];
rz(-1.9476003) q[1];
sx q[1];
rz(-0.74700673) q[1];
x q[2];
rz(-2.5164817) q[3];
sx q[3];
rz(-1.2338936) q[3];
sx q[3];
rz(0.57527104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.871375) q[2];
sx q[2];
rz(-1.3048708) q[2];
sx q[2];
rz(-2.3649141) q[2];
rz(-3.1089879) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(-1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8107574) q[0];
sx q[0];
rz(-0.36892712) q[0];
sx q[0];
rz(2.708013) q[0];
rz(2.4724204) q[1];
sx q[1];
rz(-1.1829665) q[1];
sx q[1];
rz(2.7152667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0737555) q[0];
sx q[0];
rz(-1.3569419) q[0];
sx q[0];
rz(0.20238371) q[0];
x q[1];
rz(1.2157006) q[2];
sx q[2];
rz(-0.85596687) q[2];
sx q[2];
rz(1.0155884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1703699) q[1];
sx q[1];
rz(-2.055238) q[1];
sx q[1];
rz(1.1722527) q[1];
rz(2.0570404) q[3];
sx q[3];
rz(-2.1003688) q[3];
sx q[3];
rz(-2.2467006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3907889) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(2.0354347) q[2];
rz(3.028051) q[3];
sx q[3];
rz(-0.93324408) q[3];
sx q[3];
rz(-0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15585598) q[0];
sx q[0];
rz(-0.899213) q[0];
sx q[0];
rz(-0.56094299) q[0];
rz(2.7583495) q[1];
sx q[1];
rz(-2.0189197) q[1];
sx q[1];
rz(-0.22519208) q[1];
rz(-1.8763992) q[2];
sx q[2];
rz(-1.0355179) q[2];
sx q[2];
rz(-2.1260362) q[2];
rz(-1.0998955) q[3];
sx q[3];
rz(-2.3318137) q[3];
sx q[3];
rz(-1.6398738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
