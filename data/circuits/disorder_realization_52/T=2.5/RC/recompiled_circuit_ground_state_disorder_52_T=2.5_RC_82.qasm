OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(-2.2401016) q[0];
rz(-1.787552) q[1];
sx q[1];
rz(-1.6156337) q[1];
sx q[1];
rz(-1.807133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62518277) q[0];
sx q[0];
rz(-1.3526655) q[0];
sx q[0];
rz(1.5152009) q[0];
rz(-0.15057474) q[2];
sx q[2];
rz(-1.5158487) q[2];
sx q[2];
rz(-2.6274632) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1062231) q[1];
sx q[1];
rz(-1.5179885) q[1];
sx q[1];
rz(-1.04671) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6180575) q[3];
sx q[3];
rz(-1.8128554) q[3];
sx q[3];
rz(-0.74064613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91288519) q[2];
sx q[2];
rz(-0.097147377) q[2];
sx q[2];
rz(0.65931064) q[2];
rz(0.77624503) q[3];
sx q[3];
rz(-0.018748911) q[3];
sx q[3];
rz(-0.68500486) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9126251) q[0];
sx q[0];
rz(-1.9321059) q[0];
sx q[0];
rz(1.9483161) q[0];
rz(0.044366447) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(2.9150229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3831155) q[0];
sx q[0];
rz(-1.5742013) q[0];
sx q[0];
rz(0.0016511516) q[0];
x q[1];
rz(1.5878651) q[2];
sx q[2];
rz(-1.1917795) q[2];
sx q[2];
rz(-1.5755149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5331796) q[1];
sx q[1];
rz(-0.4371818) q[1];
sx q[1];
rz(2.9730148) q[1];
rz(-2.2983589) q[3];
sx q[3];
rz(-1.7403462) q[3];
sx q[3];
rz(2.7156626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6988397) q[2];
sx q[2];
rz(-1.5719465) q[2];
sx q[2];
rz(-1.6119831) q[2];
rz(-2.2085341) q[3];
sx q[3];
rz(-1.6589087) q[3];
sx q[3];
rz(0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642692) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(2.9413057) q[0];
rz(-3.1411723) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(0.013484152) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015756135) q[0];
sx q[0];
rz(-1.5061597) q[0];
sx q[0];
rz(-1.9958853) q[0];
x q[1];
rz(-3.0749226) q[2];
sx q[2];
rz(-1.613703) q[2];
sx q[2];
rz(1.5547475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2993907) q[1];
sx q[1];
rz(-2.9889571) q[1];
sx q[1];
rz(-1.0959036) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30463574) q[3];
sx q[3];
rz(-2.9713062) q[3];
sx q[3];
rz(-1.1089604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9424332) q[2];
sx q[2];
rz(-1.5853256) q[2];
sx q[2];
rz(1.5419434) q[2];
rz(-2.13983) q[3];
sx q[3];
rz(-0.21654138) q[3];
sx q[3];
rz(2.969363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7672985) q[0];
sx q[0];
rz(-0.16574398) q[0];
sx q[0];
rz(-2.7941008) q[0];
rz(2.5047498) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(-1.1905131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40909262) q[0];
sx q[0];
rz(-0.14483368) q[0];
sx q[0];
rz(1.1454789) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6830446) q[2];
sx q[2];
rz(-0.069321037) q[2];
sx q[2];
rz(-1.6905418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3733565) q[1];
sx q[1];
rz(-2.3563085) q[1];
sx q[1];
rz(0.26765243) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3509963) q[3];
sx q[3];
rz(-1.9270718) q[3];
sx q[3];
rz(1.364546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5845329) q[2];
sx q[2];
rz(-3.1019042) q[2];
sx q[2];
rz(-1.7812799) q[2];
rz(-1.6654061) q[3];
sx q[3];
rz(-1.5675661) q[3];
sx q[3];
rz(-2.7358957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0535468) q[0];
sx q[0];
rz(-2.4021554) q[0];
sx q[0];
rz(0.0060225688) q[0];
rz(-1.7209523) q[1];
sx q[1];
rz(-3.0907478) q[1];
sx q[1];
rz(3.0555225) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1981144) q[0];
sx q[0];
rz(-3.0331342) q[0];
sx q[0];
rz(-2.9861241) q[0];
x q[1];
rz(-3.086523) q[2];
sx q[2];
rz(-1.5447642) q[2];
sx q[2];
rz(-0.30773417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7300749) q[1];
sx q[1];
rz(-0.098654276) q[1];
sx q[1];
rz(-2.7559571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5102328) q[3];
sx q[3];
rz(-1.525536) q[3];
sx q[3];
rz(1.1332944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(1.7576199) q[2];
rz(0.39517394) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(1.1672195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60642099) q[0];
sx q[0];
rz(-0.21927729) q[0];
sx q[0];
rz(-1.0004591) q[0];
rz(-2.4011627) q[1];
sx q[1];
rz(-2.705997) q[1];
sx q[1];
rz(2.7620517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1943037) q[0];
sx q[0];
rz(-1.475388) q[0];
sx q[0];
rz(1.545115) q[0];
rz(-pi) q[1];
rz(0.47490317) q[2];
sx q[2];
rz(-1.8142209) q[2];
sx q[2];
rz(0.43192513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.740447) q[1];
sx q[1];
rz(-1.5441431) q[1];
sx q[1];
rz(1.0668287) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3336181) q[3];
sx q[3];
rz(-2.6814125) q[3];
sx q[3];
rz(1.7250012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.015811054) q[2];
sx q[2];
rz(-0.28818473) q[2];
sx q[2];
rz(-3.06456) q[2];
rz(0.020126255) q[3];
sx q[3];
rz(-3.088981) q[3];
sx q[3];
rz(0.7974112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064442) q[0];
sx q[0];
rz(-3.1290717) q[0];
sx q[0];
rz(1.7092108) q[0];
rz(-0.2969946) q[1];
sx q[1];
rz(-2.98731) q[1];
sx q[1];
rz(0.061554734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3307694) q[0];
sx q[0];
rz(-0.71414882) q[0];
sx q[0];
rz(-2.215395) q[0];
x q[1];
rz(-0.21451471) q[2];
sx q[2];
rz(-1.5763064) q[2];
sx q[2];
rz(-0.44098976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0377215) q[1];
sx q[1];
rz(-1.6391409) q[1];
sx q[1];
rz(-2.8496024) q[1];
rz(-0.43470862) q[3];
sx q[3];
rz(-1.5366597) q[3];
sx q[3];
rz(-0.0061090547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2986472) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(-1.418815) q[2];
rz(-1.8055387) q[3];
sx q[3];
rz(-3.1137443) q[3];
sx q[3];
rz(1.4068039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0272738) q[0];
sx q[0];
rz(-0.71684664) q[0];
sx q[0];
rz(1.5007716) q[0];
rz(1.1227135) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.2935125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1294848) q[0];
sx q[0];
rz(-0.92381682) q[0];
sx q[0];
rz(-0.39374521) q[0];
rz(-pi) q[1];
rz(1.229153) q[2];
sx q[2];
rz(-2.2356996) q[2];
sx q[2];
rz(-1.5835012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97630097) q[1];
sx q[1];
rz(-1.2982076) q[1];
sx q[1];
rz(-1.7458245) q[1];
rz(-pi) q[2];
x q[2];
rz(2.294306) q[3];
sx q[3];
rz(-1.022911) q[3];
sx q[3];
rz(-2.5112266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5571755) q[2];
sx q[2];
rz(-2.7807005) q[2];
sx q[2];
rz(-1.7822251) q[2];
rz(0.44698295) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(-2.9744448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3249224) q[0];
sx q[0];
rz(-1.8055547) q[0];
sx q[0];
rz(-2.0409806) q[0];
rz(-2.0657516) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(0.67695391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8487602) q[0];
sx q[0];
rz(-2.7037132) q[0];
sx q[0];
rz(0.48012244) q[0];
rz(-pi) q[1];
rz(-1.4482037) q[2];
sx q[2];
rz(-1.4559146) q[2];
sx q[2];
rz(2.5455395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7266287) q[1];
sx q[1];
rz(-1.5708956) q[1];
sx q[1];
rz(-1.5703771) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4550981) q[3];
sx q[3];
rz(-1.4346204) q[3];
sx q[3];
rz(-1.0973339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0000275) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(-2.4139717) q[2];
rz(-1.242312) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(1.9696994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2321371) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(-2.452028) q[0];
rz(-1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(0.15588674) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4094226) q[0];
sx q[0];
rz(-0.66930938) q[0];
sx q[0];
rz(-2.9271892) q[0];
rz(-pi) q[1];
rz(1.6901659) q[2];
sx q[2];
rz(-1.3210591) q[2];
sx q[2];
rz(-2.2636556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2624596) q[1];
sx q[1];
rz(-0.004460881) q[1];
sx q[1];
rz(-0.35426472) q[1];
x q[2];
rz(-1.7266375) q[3];
sx q[3];
rz(-1.4771802) q[3];
sx q[3];
rz(0.35767143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.186824) q[2];
sx q[2];
rz(-0.12207741) q[2];
sx q[2];
rz(0.99511498) q[2];
rz(0.22594813) q[3];
sx q[3];
rz(-3.093231) q[3];
sx q[3];
rz(0.82609716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36289287) q[0];
sx q[0];
rz(-2.1669372) q[0];
sx q[0];
rz(-1.7397276) q[0];
rz(-1.4317935) q[1];
sx q[1];
rz(-1.8451537) q[1];
sx q[1];
rz(0.61737212) q[1];
rz(2.4564248) q[2];
sx q[2];
rz(-1.4643659) q[2];
sx q[2];
rz(2.3902389) q[2];
rz(-1.1817839) q[3];
sx q[3];
rz(-1.6970194) q[3];
sx q[3];
rz(0.57803911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
