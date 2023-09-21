OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(-0.021835672) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(2.3261855) q[0];
x q[1];
rz(2.6724042) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(2.5772622) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54675198) q[1];
sx q[1];
rz(-2.0205803) q[1];
sx q[1];
rz(0.91083281) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0825726) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-2.5773876) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3551536) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(1.2765221) q[0];
rz(-2.7361761) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(-2.237052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5337199) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-0.17287066) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0131301) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(-1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(2.7868328) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(0.43513402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1895771) q[0];
sx q[0];
rz(-2.7727685) q[0];
sx q[0];
rz(-2.8394305) q[0];
rz(-pi) q[1];
rz(1.7180195) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(-1.9927646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0355465) q[1];
sx q[1];
rz(-1.7340845) q[1];
sx q[1];
rz(1.0825023) q[1];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47539513) q[0];
sx q[0];
rz(-2.4817433) q[0];
sx q[0];
rz(-0.049113627) q[0];
rz(-pi) q[1];
rz(0.8226383) q[2];
sx q[2];
rz(-0.3974786) q[2];
sx q[2];
rz(-2.9646404) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(2.1125395) q[1];
rz(-2.3259006) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(-0.16251414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53341502) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.7081225) q[0];
rz(-pi) q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(-1.3410459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6528656) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-3.0609344) q[3];
sx q[3];
rz(-2.6622052) q[3];
sx q[3];
rz(-2.5731034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-2.3099242) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57219244) q[0];
sx q[0];
rz(-1.5844371) q[0];
sx q[0];
rz(-1.9872679) q[0];
rz(-0.2595915) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(-1.7829347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8166568) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(1.1549969) q[1];
rz(-0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(2.3849934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264383) q[0];
sx q[0];
rz(-1.7225791) q[0];
sx q[0];
rz(0.093869165) q[0];
x q[1];
rz(1.3691749) q[2];
sx q[2];
rz(-1.1668418) q[2];
sx q[2];
rz(2.902365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(1.8030333) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5114622) q[3];
sx q[3];
rz(-1.7620798) q[3];
sx q[3];
rz(-0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7635927) q[0];
sx q[0];
rz(-0.53056301) q[0];
sx q[0];
rz(-2.2647122) q[0];
x q[1];
rz(-2.9724389) q[2];
sx q[2];
rz(-0.094745853) q[2];
sx q[2];
rz(0.46846889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58395308) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(-0.072903452) q[1];
rz(-2.3905121) q[3];
sx q[3];
rz(-0.6595279) q[3];
sx q[3];
rz(-0.86576033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0489037) q[0];
sx q[0];
rz(-0.19177076) q[0];
sx q[0];
rz(1.9321241) q[0];
x q[1];
rz(2.9940669) q[2];
sx q[2];
rz(-1.4359056) q[2];
sx q[2];
rz(0.26441661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59233353) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(-1.9418282) q[1];
x q[2];
rz(2.196225) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.69797126) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012750082) q[0];
sx q[0];
rz(-1.8930952) q[0];
sx q[0];
rz(1.661257) q[0];
rz(2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(0.70891526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6203306) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(1.8950589) q[1];
x q[2];
rz(-0.20874899) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(0.97313125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-0.940154) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
rz(2.6079569) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
