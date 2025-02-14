OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1516079) q[0];
sx q[0];
rz(-0.94010544) q[0];
sx q[0];
rz(0.54036933) q[0];
rz(0.49343935) q[1];
sx q[1];
rz(-2.4210338) q[1];
sx q[1];
rz(2.5277353) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709124) q[0];
sx q[0];
rz(-1.7550091) q[0];
sx q[0];
rz(2.5066916) q[0];
x q[1];
rz(3.1356642) q[2];
sx q[2];
rz(-1.8625755) q[2];
sx q[2];
rz(-2.8242982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7481374) q[1];
sx q[1];
rz(-2.2917213) q[1];
sx q[1];
rz(1.1925405) q[1];
rz(3.0560232) q[3];
sx q[3];
rz(-2.8935379) q[3];
sx q[3];
rz(0.55094592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8453688) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(-0.63429147) q[2];
rz(-0.93572179) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(2.3989357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82035404) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(-0.4441922) q[0];
rz(0.76002899) q[1];
sx q[1];
rz(-1.1385695) q[1];
sx q[1];
rz(-0.98145032) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652577) q[0];
sx q[0];
rz(-2.0786571) q[0];
sx q[0];
rz(1.2351551) q[0];
rz(0.18155988) q[2];
sx q[2];
rz(-1.4019792) q[2];
sx q[2];
rz(-3.0491875) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5983898) q[1];
sx q[1];
rz(-1.1335422) q[1];
sx q[1];
rz(1.6679126) q[1];
rz(-pi) q[2];
rz(-2.7090453) q[3];
sx q[3];
rz(-1.3290231) q[3];
sx q[3];
rz(1.4885308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11357073) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(-1.2051955) q[2];
rz(-2.6049854) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(-1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9179012) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(2.3028288) q[0];
rz(2.5054848) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(-3.1210693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3394525) q[0];
sx q[0];
rz(-2.1198556) q[0];
sx q[0];
rz(0.86717506) q[0];
x q[1];
rz(-2.5155847) q[2];
sx q[2];
rz(-1.9591718) q[2];
sx q[2];
rz(3.1386536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.518259) q[1];
sx q[1];
rz(-1.9837399) q[1];
sx q[1];
rz(-1.5625801) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2853363) q[3];
sx q[3];
rz(-1.721964) q[3];
sx q[3];
rz(-0.54394875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3402349) q[2];
sx q[2];
rz(-1.6131718) q[2];
sx q[2];
rz(2.197263) q[2];
rz(2.1186192) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870938) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.1573855) q[0];
rz(1.6429139) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(1.9532983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12637114) q[0];
sx q[0];
rz(-2.043521) q[0];
sx q[0];
rz(-2.0822099) q[0];
x q[1];
rz(1.6158197) q[2];
sx q[2];
rz(-1.0226215) q[2];
sx q[2];
rz(-0.5612824) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21056023) q[1];
sx q[1];
rz(-2.1347441) q[1];
sx q[1];
rz(2.2497247) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1340302) q[3];
sx q[3];
rz(-2.7705857) q[3];
sx q[3];
rz(1.8414258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1184065) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(0.6905306) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(-1.5007277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(-2.114356) q[0];
rz(-1.8401624) q[1];
sx q[1];
rz(-0.65168989) q[1];
sx q[1];
rz(0.32381907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831858) q[0];
sx q[0];
rz(-2.2290285) q[0];
sx q[0];
rz(2.6290274) q[0];
rz(2.1236046) q[2];
sx q[2];
rz(-1.5587285) q[2];
sx q[2];
rz(-0.34876212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4253772) q[1];
sx q[1];
rz(-2.1105123) q[1];
sx q[1];
rz(-1.9464689) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0824329) q[3];
sx q[3];
rz(-1.4338974) q[3];
sx q[3];
rz(0.22668992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7258437) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(-0.48247639) q[2];
rz(-1.9735362) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2060858) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(-0.8859984) q[0];
rz(-0.63367263) q[1];
sx q[1];
rz(-2.251667) q[1];
sx q[1];
rz(0.69127965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871498) q[0];
sx q[0];
rz(-1.6073174) q[0];
sx q[0];
rz(-2.4958688) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4497767) q[2];
sx q[2];
rz(-1.6199281) q[2];
sx q[2];
rz(-0.092158801) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.012612494) q[1];
sx q[1];
rz(-0.75769934) q[1];
sx q[1];
rz(1.6974259) q[1];
x q[2];
rz(-1.8227303) q[3];
sx q[3];
rz(-0.79753424) q[3];
sx q[3];
rz(-3.0424527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1401691) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(0.26958618) q[2];
rz(-2.1932898) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(-1.74291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7798994) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(1.0070739) q[0];
rz(-2.7359447) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(-2.5880623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29964511) q[0];
sx q[0];
rz(-1.3163693) q[0];
sx q[0];
rz(0.05924745) q[0];
rz(-pi) q[1];
rz(0.78041665) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(-1.9165906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1055401) q[1];
sx q[1];
rz(-2.3337488) q[1];
sx q[1];
rz(2.4921662) q[1];
x q[2];
rz(-0.38000472) q[3];
sx q[3];
rz(-1.5940574) q[3];
sx q[3];
rz(-1.8398566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3125399) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(1.9276169) q[2];
rz(-0.78643262) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(-0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9455652) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.8096402) q[0];
rz(0.61344433) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(0.44949284) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025279538) q[0];
sx q[0];
rz(-0.49508245) q[0];
sx q[0];
rz(1.7465003) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.055213) q[2];
sx q[2];
rz(-1.2360209) q[2];
sx q[2];
rz(2.2074062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0969857) q[1];
sx q[1];
rz(-1.7076483) q[1];
sx q[1];
rz(-0.78949331) q[1];
rz(1.6291796) q[3];
sx q[3];
rz(-1.3135785) q[3];
sx q[3];
rz(1.8622049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.31676644) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(-0.60834926) q[2];
rz(-1.1634722) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(-0.56330645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1083531) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(-0.60229993) q[0];
rz(-2.1741518) q[1];
sx q[1];
rz(-2.2106705) q[1];
sx q[1];
rz(-2.0379351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33537827) q[0];
sx q[0];
rz(-0.69604448) q[0];
sx q[0];
rz(2.2048414) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63656143) q[2];
sx q[2];
rz(-1.5590073) q[2];
sx q[2];
rz(1.4135041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59246588) q[1];
sx q[1];
rz(-1.5696206) q[1];
sx q[1];
rz(2.0219314) q[1];
rz(3.07396) q[3];
sx q[3];
rz(-2.1478232) q[3];
sx q[3];
rz(-2.2216376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.8359567) q[2];
sx q[2];
rz(-2.804011) q[2];
rz(0.28389367) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(-2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5481446) q[0];
sx q[0];
rz(-1.5080867) q[0];
sx q[0];
rz(-0.067807587) q[0];
rz(-1.9920805) q[1];
sx q[1];
rz(-1.5905453) q[1];
sx q[1];
rz(2.6670719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.041458) q[0];
sx q[0];
rz(-0.23915072) q[0];
sx q[0];
rz(-1.6275703) q[0];
rz(1.4350008) q[2];
sx q[2];
rz(-1.8095922) q[2];
sx q[2];
rz(-2.7416624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.671512) q[1];
sx q[1];
rz(-1.5776411) q[1];
sx q[1];
rz(1.6575302) q[1];
x q[2];
rz(-0.84375937) q[3];
sx q[3];
rz(-2.2402128) q[3];
sx q[3];
rz(-2.9694174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6600251) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(-2.068326) q[2];
rz(2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58707033) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(1.5837689) q[1];
sx q[1];
rz(-1.0777892) q[1];
sx q[1];
rz(-0.52660175) q[1];
rz(0.46109328) q[2];
sx q[2];
rz(-1.8965707) q[2];
sx q[2];
rz(-2.0988219) q[2];
rz(-0.14231331) q[3];
sx q[3];
rz(-1.489779) q[3];
sx q[3];
rz(-1.6708556) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
