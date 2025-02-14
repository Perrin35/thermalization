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
rz(2.4647291) q[0];
sx q[0];
rz(-1.530175) q[0];
sx q[0];
rz(2.8443008) q[0];
rz(0.77904207) q[1];
sx q[1];
rz(-0.160633) q[1];
sx q[1];
rz(1.4359441) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061062) q[0];
sx q[0];
rz(-1.0260884) q[0];
sx q[0];
rz(-2.5404055) q[0];
x q[1];
rz(3.0129635) q[2];
sx q[2];
rz(-1.8278215) q[2];
sx q[2];
rz(-0.84398735) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7117225) q[1];
sx q[1];
rz(-0.87583576) q[1];
sx q[1];
rz(1.9942787) q[1];
rz(-pi) q[2];
rz(1.7706031) q[3];
sx q[3];
rz(-1.8733896) q[3];
sx q[3];
rz(1.8963008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5433189) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-0.16779009) q[2];
rz(2.0134036) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096864916) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(-2.0665533) q[0];
rz(-1.3661522) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(1.8278106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9843861) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(-2.8048672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4046217) q[2];
sx q[2];
rz(-2.8923375) q[2];
sx q[2];
rz(-2.1971306) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35006501) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(-1.6863053) q[1];
rz(-pi) q[2];
rz(1.9017392) q[3];
sx q[3];
rz(-1.6972491) q[3];
sx q[3];
rz(0.053205333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4904867) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(0.42438486) q[2];
rz(-0.73244798) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(-1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(-2.6388229) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(-0.81519333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8250587) q[0];
sx q[0];
rz(-1.382917) q[0];
sx q[0];
rz(1.5231537) q[0];
x q[1];
rz(0.72702144) q[2];
sx q[2];
rz(-0.15910782) q[2];
sx q[2];
rz(0.66691676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5126172) q[1];
sx q[1];
rz(-1.830258) q[1];
sx q[1];
rz(-1.2995059) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.268154) q[3];
sx q[3];
rz(-1.1213574) q[3];
sx q[3];
rz(3.0118503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.076685585) q[2];
sx q[2];
rz(-1.6904597) q[2];
sx q[2];
rz(-2.0150851) q[2];
rz(1.3569776) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55946881) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(2.7816787) q[0];
rz(-2.9008046) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(2.7684033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2997595) q[0];
sx q[0];
rz(-1.4275121) q[0];
sx q[0];
rz(-3.0413607) q[0];
x q[1];
rz(1.8055438) q[2];
sx q[2];
rz(-2.1214607) q[2];
sx q[2];
rz(-2.5900813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2149005) q[1];
sx q[1];
rz(-2.5040031) q[1];
sx q[1];
rz(0.59507782) q[1];
rz(-pi) q[2];
rz(2.7591428) q[3];
sx q[3];
rz(-0.47690629) q[3];
sx q[3];
rz(-2.675569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.084426247) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(0.90844321) q[2];
rz(-1.069979) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(-0.98328868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(2.4519582) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(0.77401179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871983) q[0];
sx q[0];
rz(-2.3582705) q[0];
sx q[0];
rz(-0.67449595) q[0];
x q[1];
rz(3.1104964) q[2];
sx q[2];
rz(-2.3457639) q[2];
sx q[2];
rz(2.9253256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8929114) q[1];
sx q[1];
rz(-1.9982583) q[1];
sx q[1];
rz(-2.9315154) q[1];
x q[2];
rz(1.6004531) q[3];
sx q[3];
rz(-2.1184485) q[3];
sx q[3];
rz(2.1911007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8959877) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(0.73949933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6226115) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(-0.13825026) q[1];
sx q[1];
rz(-1.6743276) q[1];
sx q[1];
rz(1.5672055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7049887) q[0];
sx q[0];
rz(-1.2452176) q[0];
sx q[0];
rz(-2.7537936) q[0];
rz(-pi) q[1];
rz(0.60092314) q[2];
sx q[2];
rz(-0.4888566) q[2];
sx q[2];
rz(-0.86684858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3688599) q[1];
sx q[1];
rz(-2.8535758) q[1];
sx q[1];
rz(2.6898328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54924567) q[3];
sx q[3];
rz(-1.4016782) q[3];
sx q[3];
rz(1.0602407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77341998) q[2];
sx q[2];
rz(-1.0680826) q[2];
sx q[2];
rz(-2.0424021) q[2];
rz(2.1558971) q[3];
sx q[3];
rz(-1.516195) q[3];
sx q[3];
rz(-2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8428335) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(2.47593) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(0.98141247) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.736826) q[0];
sx q[0];
rz(-2.5676709) q[0];
sx q[0];
rz(1.4653652) q[0];
x q[1];
rz(-1.3150042) q[2];
sx q[2];
rz(-2.062254) q[2];
sx q[2];
rz(1.0813528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1011103) q[1];
sx q[1];
rz(-2.6553934) q[1];
sx q[1];
rz(-1.2771439) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3990226) q[3];
sx q[3];
rz(-2.4618759) q[3];
sx q[3];
rz(-1.8006067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0770646) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(2.7725753) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(-1.3212475) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(-0.33018026) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(0.59828573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56494026) q[0];
sx q[0];
rz(-0.9441174) q[0];
sx q[0];
rz(-1.0231072) q[0];
x q[1];
rz(1.1016229) q[2];
sx q[2];
rz(-1.507505) q[2];
sx q[2];
rz(-0.67055145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48409325) q[1];
sx q[1];
rz(-1.7998905) q[1];
sx q[1];
rz(0.34219663) q[1];
rz(-0.98294799) q[3];
sx q[3];
rz(-2.1390669) q[3];
sx q[3];
rz(0.98230386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-0.064621933) q[2];
rz(-1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(-1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794063) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(0.52126467) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(-2.1102139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9654664) q[0];
sx q[0];
rz(-2.3076008) q[0];
sx q[0];
rz(-1.6710207) q[0];
rz(-pi) q[1];
rz(1.6907482) q[2];
sx q[2];
rz(-1.861683) q[2];
sx q[2];
rz(-1.6294668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0887275) q[1];
sx q[1];
rz(-1.8590312) q[1];
sx q[1];
rz(-1.0239787) q[1];
rz(1.4724588) q[3];
sx q[3];
rz(-1.0957484) q[3];
sx q[3];
rz(2.6779384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3970268) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(0.087372027) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(-0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(1.1391621) q[0];
rz(-2.6423404) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(-1.8081236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718057) q[0];
sx q[0];
rz(-3.0020368) q[0];
sx q[0];
rz(0.7044756) q[0];
rz(1.8257481) q[2];
sx q[2];
rz(-1.6189908) q[2];
sx q[2];
rz(-0.15010897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.002961) q[1];
sx q[1];
rz(-1.3182782) q[1];
sx q[1];
rz(1.7934402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74204294) q[3];
sx q[3];
rz(-2.2931406) q[3];
sx q[3];
rz(-0.53398856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.501005) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(2.8583543) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3129811) q[0];
sx q[0];
rz(-0.77684488) q[0];
sx q[0];
rz(-2.0364398) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(1.3731643) q[2];
sx q[2];
rz(-1.0734954) q[2];
sx q[2];
rz(-1.3263477) q[2];
rz(0.5091359) q[3];
sx q[3];
rz(-0.73496277) q[3];
sx q[3];
rz(-3.1258813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
