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
rz(-1.6753766) q[0];
sx q[0];
rz(-2.197062) q[0];
sx q[0];
rz(2.9401927) q[0];
rz(-2.0693076) q[1];
sx q[1];
rz(-2.3787002) q[1];
sx q[1];
rz(1.720517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147707) q[0];
sx q[0];
rz(-2.1807007) q[0];
sx q[0];
rz(3.0746835) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95008738) q[2];
sx q[2];
rz(-0.7751152) q[2];
sx q[2];
rz(-2.895854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3609429) q[1];
sx q[1];
rz(-0.93027481) q[1];
sx q[1];
rz(-2.065413) q[1];
x q[2];
rz(-0.87405494) q[3];
sx q[3];
rz(-0.49421453) q[3];
sx q[3];
rz(-0.53411247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(-0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-0.24638677) q[0];
rz(-1.9465744) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.5066719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6794066) q[0];
sx q[0];
rz(-0.52216086) q[0];
sx q[0];
rz(2.0798111) q[0];
x q[1];
rz(3.0987034) q[2];
sx q[2];
rz(-1.2681307) q[2];
sx q[2];
rz(-1.3009225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8283702) q[1];
sx q[1];
rz(-1.5130318) q[1];
sx q[1];
rz(-1.934113) q[1];
x q[2];
rz(-1.0492658) q[3];
sx q[3];
rz(-1.9535011) q[3];
sx q[3];
rz(-2.692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40111497) q[2];
sx q[2];
rz(-2.1666708) q[2];
sx q[2];
rz(-1.9096036) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6983637) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(0.90782905) q[0];
rz(-2.5746131) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(0.10279113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55926052) q[0];
sx q[0];
rz(-1.0114504) q[0];
sx q[0];
rz(1.8190967) q[0];
x q[1];
rz(-1.6750828) q[2];
sx q[2];
rz(-1.490218) q[2];
sx q[2];
rz(-1.8957418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56270617) q[1];
sx q[1];
rz(-1.1670615) q[1];
sx q[1];
rz(-2.5078234) q[1];
rz(1.0166753) q[3];
sx q[3];
rz(-1.8401056) q[3];
sx q[3];
rz(2.179972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28321442) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.7041448) q[2];
rz(0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(-0.2984305) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0933541) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(-1.5123222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8522569) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(-2.1136978) q[0];
x q[1];
rz(-0.34660201) q[2];
sx q[2];
rz(-0.79814974) q[2];
sx q[2];
rz(-0.59890282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5315631) q[1];
sx q[1];
rz(-1.3852784) q[1];
sx q[1];
rz(-0.97037913) q[1];
rz(-pi) q[2];
rz(-1.063594) q[3];
sx q[3];
rz(-1.385687) q[3];
sx q[3];
rz(1.1066268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-2.4646087) q[2];
rz(1.8949159) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(-1.6251132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077654) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(-3.1166792) q[0];
rz(-2.086153) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(0.28344646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0830529) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(1.3015981) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7666498) q[2];
sx q[2];
rz(-1.699866) q[2];
sx q[2];
rz(-1.7650676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9495715) q[1];
sx q[1];
rz(-1.6013491) q[1];
sx q[1];
rz(-1.0699238) q[1];
x q[2];
rz(2.2437566) q[3];
sx q[3];
rz(-1.3011609) q[3];
sx q[3];
rz(1.8561038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8289566) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(-0.32583315) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(-1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40389898) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(-2.7666336) q[0];
rz(0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(-0.37014827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3309114) q[0];
sx q[0];
rz(-3.1157012) q[0];
sx q[0];
rz(-2.6736892) q[0];
rz(-pi) q[1];
rz(-3.0022718) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(-0.56688165) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1620336) q[1];
sx q[1];
rz(-2.0311072) q[1];
sx q[1];
rz(2.3737337) q[1];
rz(-pi) q[2];
rz(3.1305976) q[3];
sx q[3];
rz(-2.7137626) q[3];
sx q[3];
rz(-1.4233936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(-0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8884856) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-2.2174368) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.7402657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3649021) q[0];
sx q[0];
rz(-0.49743891) q[0];
sx q[0];
rz(-3.0203392) q[0];
x q[1];
rz(-2.9193814) q[2];
sx q[2];
rz(-2.4743818) q[2];
sx q[2];
rz(-1.8052342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3360916) q[1];
sx q[1];
rz(-0.84118836) q[1];
sx q[1];
rz(1.7859687) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21715607) q[3];
sx q[3];
rz(-2.0494866) q[3];
sx q[3];
rz(1.9793122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9132729) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(-1.9165967) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55815721) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(2.9971314) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(-2.2772148) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4557171) q[0];
sx q[0];
rz(-0.17016958) q[0];
sx q[0];
rz(1.8491114) q[0];
rz(-pi) q[1];
rz(-0.34274613) q[2];
sx q[2];
rz(-0.10048332) q[2];
sx q[2];
rz(-3.036694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64329735) q[1];
sx q[1];
rz(-1.3302696) q[1];
sx q[1];
rz(1.0161736) q[1];
x q[2];
rz(1.8501353) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(3.0077028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(2.5065191) q[2];
rz(-0.57279974) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(0.12426678) q[0];
rz(-0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(1.7100547) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99125049) q[0];
sx q[0];
rz(-1.7479154) q[0];
sx q[0];
rz(2.9464673) q[0];
rz(0.65406873) q[2];
sx q[2];
rz(-1.1032915) q[2];
sx q[2];
rz(-1.2142177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6874976) q[1];
sx q[1];
rz(-2.3907067) q[1];
sx q[1];
rz(-0.45529699) q[1];
x q[2];
rz(1.133119) q[3];
sx q[3];
rz(-2.056155) q[3];
sx q[3];
rz(0.35972586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.734599) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.2146568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491972) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(2.5469575) q[0];
rz(-2.1356964) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(2.2524021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85814637) q[0];
sx q[0];
rz(-2.915417) q[0];
sx q[0];
rz(-1.3429789) q[0];
rz(1.4290733) q[2];
sx q[2];
rz(-0.46438875) q[2];
sx q[2];
rz(2.1187256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8831625) q[1];
sx q[1];
rz(-0.11181242) q[1];
sx q[1];
rz(-2.1360141) q[1];
rz(-0.51639207) q[3];
sx q[3];
rz(-1.603873) q[3];
sx q[3];
rz(2.3865478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1903926) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(-0.54538837) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(1.4639328) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496898) q[0];
sx q[0];
rz(-0.90072537) q[0];
sx q[0];
rz(-1.8186722) q[0];
rz(-0.84359618) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(0.19914535) q[2];
sx q[2];
rz(-1.1391339) q[2];
sx q[2];
rz(1.7720411) q[2];
rz(-0.69915184) q[3];
sx q[3];
rz(-0.71157645) q[3];
sx q[3];
rz(-0.012307766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
