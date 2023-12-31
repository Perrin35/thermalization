OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4364606) q[0];
sx q[0];
rz(-0.55186614) q[0];
sx q[0];
rz(-3.119757) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37786814) q[0];
sx q[0];
rz(-1.8452497) q[0];
sx q[0];
rz(0.26835693) q[0];
rz(-pi) q[1];
rz(0.46918842) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-0.56433041) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5948407) q[1];
sx q[1];
rz(-2.0205803) q[1];
sx q[1];
rz(2.2307598) q[1];
rz(2.4258852) q[3];
sx q[3];
rz(-1.9519836) q[3];
sx q[3];
rz(0.79431278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-2.2136097) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0052764) q[0];
sx q[0];
rz(-1.7668084) q[0];
sx q[0];
rz(0.85640237) q[0];
rz(1.6438269) q[2];
sx q[2];
rz(-1.4023997) q[2];
sx q[2];
rz(-1.8254691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36205081) q[1];
sx q[1];
rz(-0.78501399) q[1];
sx q[1];
rz(1.7464459) q[1];
rz(-0.41665839) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(0.43513402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4774287) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(-0.35350032) q[0];
rz(1.7180195) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(-1.148828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16776925) q[1];
sx q[1];
rz(-0.5127738) q[1];
sx q[1];
rz(-1.9085401) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76336236) q[3];
sx q[3];
rz(-1.8244787) q[3];
sx q[3];
rz(-2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-3.0505676) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19642297) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(0.26487574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0565856) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(2.4823275) q[0];
rz(-pi) q[1];
rz(0.27819602) q[2];
sx q[2];
rz(-1.2831266) q[2];
sx q[2];
rz(-2.1759335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.781573) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(-0.54846958) q[1];
x q[2];
rz(-0.78062765) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(-1.2300223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(2.662861) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.032741) q[0];
sx q[0];
rz(-1.4534338) q[0];
sx q[0];
rz(2.5928241) q[0];
rz(-1.9031992) q[2];
sx q[2];
rz(-1.9661511) q[2];
sx q[2];
rz(-2.7796641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.5336256) q[3];
sx q[3];
rz(0.93070785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7197363) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(-2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-2.9690202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1369551) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(-0.01491551) q[0];
x q[1];
rz(2.0667604) q[2];
sx q[2];
rz(-0.50903382) q[2];
sx q[2];
rz(2.3376655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8494106) q[1];
sx q[1];
rz(-1.1574355) q[1];
sx q[1];
rz(0.11489111) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8014088) q[3];
sx q[3];
rz(-2.1884544) q[3];
sx q[3];
rz(2.2889083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-0.75659928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7998357) q[0];
sx q[0];
rz(-0.17827398) q[0];
sx q[0];
rz(-1.0210277) q[0];
rz(-2.7034764) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(-0.24030906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36180624) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(-0.96868412) q[1];
x q[2];
rz(2.9499801) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38935223) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(-0.35895343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.093401508) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(1.2707368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5576396) q[1];
sx q[1];
rz(-1.156731) q[1];
sx q[1];
rz(-0.072903452) q[1];
x q[2];
rz(-2.057468) q[3];
sx q[3];
rz(-2.0351962) q[3];
sx q[3];
rz(1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-2.5323903) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163923) q[0];
sx q[0];
rz(-1.3915477) q[0];
sx q[0];
rz(0.068530131) q[0];
rz(-pi) q[1];
rz(-0.14752578) q[2];
sx q[2];
rz(-1.4359056) q[2];
sx q[2];
rz(0.26441661) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59233353) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(-1.1997644) q[1];
x q[2];
rz(2.14823) q[3];
sx q[3];
rz(-1.9462898) q[3];
sx q[3];
rz(1.7253699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436214) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.7451161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6122702) q[0];
sx q[0];
rz(-1.6565874) q[0];
sx q[0];
rz(-0.32353185) q[0];
rz(-pi) q[1];
rz(1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(-1.0422848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1435946) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(0.29923156) q[1];
rz(-pi) q[2];
rz(0.26225984) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(-3.0083169) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(2.9524654) q[2];
sx q[2];
rz(-2.1928939) q[2];
sx q[2];
rz(2.1415276) q[2];
rz(1.6205447) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
