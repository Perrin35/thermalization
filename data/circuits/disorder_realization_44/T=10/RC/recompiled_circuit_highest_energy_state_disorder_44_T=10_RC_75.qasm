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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(-1.1535147) q[0];
rz(1.310362) q[1];
sx q[1];
rz(-0.49340931) q[1];
sx q[1];
rz(-0.44841132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2758688) q[0];
sx q[0];
rz(-2.4910036) q[0];
sx q[0];
rz(-3.0537729) q[0];
rz(2.885778) q[2];
sx q[2];
rz(-1.4728776) q[2];
sx q[2];
rz(2.4061587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5858305) q[1];
sx q[1];
rz(-0.96802789) q[1];
sx q[1];
rz(2.4911111) q[1];
rz(0.5011933) q[3];
sx q[3];
rz(-1.6128916) q[3];
sx q[3];
rz(0.30748366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6260234) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-2.5173729) q[2];
rz(2.7624687) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17495951) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(2.6587291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0014526) q[0];
sx q[0];
rz(-1.4743544) q[0];
sx q[0];
rz(1.524855) q[0];
rz(2.9752983) q[2];
sx q[2];
rz(-1.4539084) q[2];
sx q[2];
rz(-0.4695732) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.711593) q[1];
sx q[1];
rz(-1.5614127) q[1];
sx q[1];
rz(1.1389772) q[1];
x q[2];
rz(0.24077529) q[3];
sx q[3];
rz(-2.8117883) q[3];
sx q[3];
rz(2.2036116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1172993) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(0.67000669) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(-2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4033177) q[0];
sx q[0];
rz(-2.6955695) q[0];
sx q[0];
rz(0.73905149) q[0];
rz(-0.96145472) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(-0.75698537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805646) q[0];
sx q[0];
rz(-1.3098728) q[0];
sx q[0];
rz(0.83103128) q[0];
rz(-pi) q[1];
rz(0.44481014) q[2];
sx q[2];
rz(-0.74193566) q[2];
sx q[2];
rz(0.43659376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.095159141) q[1];
sx q[1];
rz(-2.1951402) q[1];
sx q[1];
rz(-1.7364362) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0717172) q[3];
sx q[3];
rz(-1.8777851) q[3];
sx q[3];
rz(1.2941293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5028533) q[2];
sx q[2];
rz(-0.91705489) q[2];
sx q[2];
rz(0.70510954) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(-2.7307935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(-1.2257082) q[0];
rz(-1.2340087) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(3.0753678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561838) q[0];
sx q[0];
rz(-1.8161621) q[0];
sx q[0];
rz(-0.092459926) q[0];
x q[1];
rz(2.340256) q[2];
sx q[2];
rz(-1.3532012) q[2];
sx q[2];
rz(0.8720397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5410616) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(-0.97452428) q[1];
x q[2];
rz(2.8230571) q[3];
sx q[3];
rz(-2.7580166) q[3];
sx q[3];
rz(2.2936402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0951198) q[2];
sx q[2];
rz(-2.1472223) q[2];
sx q[2];
rz(0.44060102) q[2];
rz(-1.0062086) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42136583) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(1.5059858) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(1.45586) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57350006) q[0];
sx q[0];
rz(-1.716189) q[0];
sx q[0];
rz(-2.306434) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6728476) q[2];
sx q[2];
rz(-0.25368099) q[2];
sx q[2];
rz(-1.4616592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8090022) q[1];
sx q[1];
rz(-0.9192411) q[1];
sx q[1];
rz(-1.5963735) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7248575) q[3];
sx q[3];
rz(-2.1733781) q[3];
sx q[3];
rz(-1.3101206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(-2.2301162) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(0.13033303) q[0];
rz(-2.6538972) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(0.74388751) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7672507) q[0];
sx q[0];
rz(-3.0453186) q[0];
sx q[0];
rz(-2.1153347) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0024372) q[2];
sx q[2];
rz(-0.94991131) q[2];
sx q[2];
rz(0.54229743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1229314) q[1];
sx q[1];
rz(-2.5187188) q[1];
sx q[1];
rz(-2.0707002) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7508932) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(-0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0222212) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(-2.7458701) q[3];
sx q[3];
rz(-1.3362249) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845067) q[0];
sx q[0];
rz(-2.0261903) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1395124) q[0];
sx q[0];
rz(-2.2497228) q[0];
sx q[0];
rz(-2.9154791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.057633295) q[2];
sx q[2];
rz(-2.1890292) q[2];
sx q[2];
rz(1.0317486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8128374) q[1];
sx q[1];
rz(-2.8649674) q[1];
sx q[1];
rz(-3.0429716) q[1];
rz(-pi) q[2];
rz(0.32870558) q[3];
sx q[3];
rz(-1.8227326) q[3];
sx q[3];
rz(3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(1.3537815) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-2.8712809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29276174) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(2.8344179) q[0];
rz(-1.3245026) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(-2.2408392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7194289) q[0];
sx q[0];
rz(-1.3086645) q[0];
sx q[0];
rz(-1.5533226) q[0];
x q[1];
rz(-0.45489725) q[2];
sx q[2];
rz(-2.1730246) q[2];
sx q[2];
rz(-1.0554316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0010502) q[1];
sx q[1];
rz(-1.2549572) q[1];
sx q[1];
rz(1.2211826) q[1];
rz(-pi) q[2];
rz(0.19260223) q[3];
sx q[3];
rz(-0.64036548) q[3];
sx q[3];
rz(-2.7964724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(1.0963415) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480302) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(2.6322741) q[0];
rz(2.9604984) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(0.96955713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80705331) q[0];
sx q[0];
rz(-1.7248123) q[0];
sx q[0];
rz(2.3898983) q[0];
rz(0.050886919) q[2];
sx q[2];
rz(-2.0774088) q[2];
sx q[2];
rz(-1.8614872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6020376) q[1];
sx q[1];
rz(-1.6786715) q[1];
sx q[1];
rz(-2.8007568) q[1];
rz(-pi) q[2];
rz(1.9593616) q[3];
sx q[3];
rz(-0.7668743) q[3];
sx q[3];
rz(2.9076613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8934882) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(0.23165101) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-0.20802465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(-2.6864247) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(0.761935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021558048) q[0];
sx q[0];
rz(-1.2823449) q[0];
sx q[0];
rz(-0.7201654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3222515) q[2];
sx q[2];
rz(-2.0692673) q[2];
sx q[2];
rz(-0.33453951) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.476452) q[1];
sx q[1];
rz(-2.2287205) q[1];
sx q[1];
rz(-1.614078) q[1];
rz(-pi) q[2];
rz(-1.5860709) q[3];
sx q[3];
rz(-2.6741283) q[3];
sx q[3];
rz(1.8031507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.053293856) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(-0.90395149) q[2];
rz(-1.5229185) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(-1.9797549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.224613) q[0];
sx q[0];
rz(-1.1876748) q[0];
sx q[0];
rz(-1.351958) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(1.1197208) q[2];
sx q[2];
rz(-1.7985761) q[2];
sx q[2];
rz(0.41440415) q[2];
rz(-1.5736754) q[3];
sx q[3];
rz(-1.1408014) q[3];
sx q[3];
rz(-0.037527966) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
