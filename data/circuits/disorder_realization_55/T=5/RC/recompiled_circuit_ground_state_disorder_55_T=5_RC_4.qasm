OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9094698) q[0];
sx q[0];
rz(-2.1535518) q[0];
sx q[0];
rz(2.7890132) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-0.7847473) q[1];
sx q[1];
rz(0.027332505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77307149) q[0];
sx q[0];
rz(-1.7681244) q[0];
sx q[0];
rz(-0.34223603) q[0];
rz(-pi) q[1];
rz(0.17345239) q[2];
sx q[2];
rz(-1.7043742) q[2];
sx q[2];
rz(-2.0138149) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.057291659) q[1];
sx q[1];
rz(-2.9708869) q[1];
sx q[1];
rz(-2.956326) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23689278) q[3];
sx q[3];
rz(-1.7125355) q[3];
sx q[3];
rz(-0.6976034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0722384) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(1.5864774) q[2];
rz(0.72748264) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(-2.6684842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.469406) q[0];
sx q[0];
rz(-1.3453901) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-0.32497111) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(0.51345888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8425861) q[0];
sx q[0];
rz(-0.76548274) q[0];
sx q[0];
rz(2.1730636) q[0];
rz(-0.32294257) q[2];
sx q[2];
rz(-0.86951288) q[2];
sx q[2];
rz(-0.27837929) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8850234) q[1];
sx q[1];
rz(-2.3833439) q[1];
sx q[1];
rz(0.77039657) q[1];
rz(-2.6504032) q[3];
sx q[3];
rz(-2.0955387) q[3];
sx q[3];
rz(0.82651807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010633858) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(1.0884253) q[2];
rz(0.66793495) q[3];
sx q[3];
rz(-2.9974944) q[3];
sx q[3];
rz(1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82069355) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(1.8841085) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(0.65748293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456928) q[0];
sx q[0];
rz(-0.4687216) q[0];
sx q[0];
rz(-0.29229852) q[0];
x q[1];
rz(-2.3917603) q[2];
sx q[2];
rz(-0.74374357) q[2];
sx q[2];
rz(0.86803555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0850683) q[1];
sx q[1];
rz(-2.6303362) q[1];
sx q[1];
rz(-0.35447094) q[1];
rz(-pi) q[2];
x q[2];
rz(2.231953) q[3];
sx q[3];
rz(-1.8170905) q[3];
sx q[3];
rz(1.1681739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.377044) q[2];
sx q[2];
rz(-0.42082861) q[2];
sx q[2];
rz(1.2505038) q[2];
rz(1.4416384) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(-2.3120248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68469754) q[0];
sx q[0];
rz(-2.1853515) q[0];
sx q[0];
rz(-0.60234219) q[0];
rz(-0.96386987) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(-2.9197555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5103127) q[0];
sx q[0];
rz(-1.5820192) q[0];
sx q[0];
rz(1.5785286) q[0];
rz(-1.9811976) q[2];
sx q[2];
rz(-0.56804915) q[2];
sx q[2];
rz(-0.33524738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79961046) q[1];
sx q[1];
rz(-2.1253808) q[1];
sx q[1];
rz(-2.0444524) q[1];
x q[2];
rz(-2.9939566) q[3];
sx q[3];
rz(-2.6872928) q[3];
sx q[3];
rz(0.67377485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7092789) q[2];
sx q[2];
rz(-1.9807434) q[2];
sx q[2];
rz(-3.0529037) q[2];
rz(-2.6426219) q[3];
sx q[3];
rz(-1.6224529) q[3];
sx q[3];
rz(-3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(-0.29014507) q[1];
sx q[1];
rz(-2.2369308) q[1];
sx q[1];
rz(-1.9754999) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13577382) q[0];
sx q[0];
rz(-0.5899094) q[0];
sx q[0];
rz(-0.29340549) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8639216) q[2];
sx q[2];
rz(-1.8050965) q[2];
sx q[2];
rz(-1.9512343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2578967) q[1];
sx q[1];
rz(-1.9120815) q[1];
sx q[1];
rz(-2.8410556) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5968634) q[3];
sx q[3];
rz(-1.3192303) q[3];
sx q[3];
rz(2.9575916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0151998) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(1.8533206) q[2];
rz(1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.4748632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7159336) q[0];
sx q[0];
rz(-3.0379744) q[0];
sx q[0];
rz(-0.15956751) q[0];
rz(-2.5345934) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(-1.0608231) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.589405) q[0];
sx q[0];
rz(-0.97597835) q[0];
sx q[0];
rz(-2.908559) q[0];
x q[1];
rz(-2.4446746) q[2];
sx q[2];
rz(-2.7978511) q[2];
sx q[2];
rz(-0.70408193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86058206) q[1];
sx q[1];
rz(-1.4352903) q[1];
sx q[1];
rz(2.1675088) q[1];
x q[2];
rz(-2.5626866) q[3];
sx q[3];
rz(-1.488655) q[3];
sx q[3];
rz(-2.3711747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(0.38786495) q[2];
rz(-1.806949) q[3];
sx q[3];
rz(-2.4692061) q[3];
sx q[3];
rz(2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(2.2482596) q[0];
rz(0.54126254) q[1];
sx q[1];
rz(-2.2365139) q[1];
sx q[1];
rz(1.062324) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60266274) q[0];
sx q[0];
rz(-1.3100123) q[0];
sx q[0];
rz(2.3988924) q[0];
rz(1.5028033) q[2];
sx q[2];
rz(-0.74882245) q[2];
sx q[2];
rz(-0.48922005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7819643) q[1];
sx q[1];
rz(-2.5853695) q[1];
sx q[1];
rz(0.044597711) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19936187) q[3];
sx q[3];
rz(-2.3219675) q[3];
sx q[3];
rz(0.79341753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5541151) q[2];
sx q[2];
rz(-0.97439659) q[2];
sx q[2];
rz(-1.1678196) q[2];
rz(0.81985146) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(-2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.576936) q[0];
sx q[0];
rz(-2.7695203) q[0];
sx q[0];
rz(-1.8336953) q[0];
rz(-1.6851743) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(3.0110722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3878138) q[0];
sx q[0];
rz(-2.8243833) q[0];
sx q[0];
rz(1.3340201) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47696205) q[2];
sx q[2];
rz(-0.85627189) q[2];
sx q[2];
rz(0.80452308) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3106892) q[1];
sx q[1];
rz(-0.92055087) q[1];
sx q[1];
rz(1.6199153) q[1];
rz(2.4198615) q[3];
sx q[3];
rz(-2.2829307) q[3];
sx q[3];
rz(-0.16068383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9556344) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(1.4107417) q[2];
rz(0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(1.5131148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677583) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(0.14061418) q[0];
rz(-2.4070542) q[1];
sx q[1];
rz(-1.7081552) q[1];
sx q[1];
rz(0.75070423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561309) q[0];
sx q[0];
rz(-1.805725) q[0];
sx q[0];
rz(2.8402907) q[0];
rz(-1.9173938) q[2];
sx q[2];
rz(-1.3470413) q[2];
sx q[2];
rz(-0.72347298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.46671154) q[1];
sx q[1];
rz(-1.3916236) q[1];
sx q[1];
rz(-0.06249832) q[1];
x q[2];
rz(-1.4846913) q[3];
sx q[3];
rz(-0.70843452) q[3];
sx q[3];
rz(-1.751386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.004403) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(-2.3889551) q[2];
rz(0.33958069) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(-0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362157) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(2.5035653) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(-0.80150882) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0251533) q[0];
sx q[0];
rz(-2.354216) q[0];
sx q[0];
rz(1.8579432) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24157816) q[2];
sx q[2];
rz(-0.8525341) q[2];
sx q[2];
rz(3.0294543) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90097016) q[1];
sx q[1];
rz(-0.83604807) q[1];
sx q[1];
rz(3.088984) q[1];
rz(-pi) q[2];
rz(-1.784104) q[3];
sx q[3];
rz(-2.4856728) q[3];
sx q[3];
rz(1.3154958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30404299) q[2];
sx q[2];
rz(-1.3648698) q[2];
sx q[2];
rz(-1.1207646) q[2];
rz(-0.70901999) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(1.9161061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365702) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(2.2364521) q[1];
sx q[1];
rz(-0.97733472) q[1];
sx q[1];
rz(2.6505145) q[1];
rz(2.7837743) q[2];
sx q[2];
rz(-0.87643788) q[2];
sx q[2];
rz(-1.4446148) q[2];
rz(2.8560588) q[3];
sx q[3];
rz(-2.0750506) q[3];
sx q[3];
rz(-1.688523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
