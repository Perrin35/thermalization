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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(0.65879917) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(0.08858362) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1155498) q[0];
sx q[0];
rz(-2.6691873) q[0];
sx q[0];
rz(-2.0055072) q[0];
rz(-pi) q[1];
rz(0.93841432) q[2];
sx q[2];
rz(-2.4749196) q[2];
sx q[2];
rz(3.0423903) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40868739) q[1];
sx q[1];
rz(-2.5026606) q[1];
sx q[1];
rz(1.8270135) q[1];
rz(-2.1351571) q[3];
sx q[3];
rz(-0.86216037) q[3];
sx q[3];
rz(2.9624727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52458557) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(2.6535772) q[2];
rz(2.6739142) q[3];
sx q[3];
rz(-2.6165163) q[3];
sx q[3];
rz(-2.974143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2117598) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(1.2540586) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(2.7235203) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195898) q[0];
sx q[0];
rz(-1.749335) q[0];
sx q[0];
rz(-2.5994632) q[0];
rz(-pi) q[1];
rz(0.16440161) q[2];
sx q[2];
rz(-2.4864956) q[2];
sx q[2];
rz(3.0384105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7217191) q[1];
sx q[1];
rz(-1.5054107) q[1];
sx q[1];
rz(-1.8871149) q[1];
rz(-2.5097519) q[3];
sx q[3];
rz(-0.92636331) q[3];
sx q[3];
rz(-1.0570248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40689251) q[2];
sx q[2];
rz(-2.0105346) q[2];
sx q[2];
rz(-0.1184173) q[2];
rz(-0.40220574) q[3];
sx q[3];
rz(-1.6536568) q[3];
sx q[3];
rz(-1.3291298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44620946) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(1.1545908) q[0];
rz(-1.062695) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(-1.9564995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9364734) q[0];
sx q[0];
rz(-1.2944049) q[0];
sx q[0];
rz(-2.2689681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0555004) q[2];
sx q[2];
rz(-1.764311) q[2];
sx q[2];
rz(-0.12963477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94471632) q[1];
sx q[1];
rz(-1.8507974) q[1];
sx q[1];
rz(1.9619322) q[1];
x q[2];
rz(2.5873379) q[3];
sx q[3];
rz(-0.46559428) q[3];
sx q[3];
rz(-1.8150282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2350754) q[2];
sx q[2];
rz(-0.86248988) q[2];
sx q[2];
rz(1.4441351) q[2];
rz(-0.92153543) q[3];
sx q[3];
rz(-0.60465616) q[3];
sx q[3];
rz(3.0864033) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7017355) q[0];
sx q[0];
rz(-3.0210962) q[0];
sx q[0];
rz(-0.078027092) q[0];
rz(-1.1174508) q[1];
sx q[1];
rz(-1.2470587) q[1];
sx q[1];
rz(1.4720346) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3465418) q[0];
sx q[0];
rz(-1.4637837) q[0];
sx q[0];
rz(-2.3404164) q[0];
x q[1];
rz(0.98791583) q[2];
sx q[2];
rz(-2.0327518) q[2];
sx q[2];
rz(-0.5505901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7546447) q[1];
sx q[1];
rz(-2.7803763) q[1];
sx q[1];
rz(-1.4903461) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4185798) q[3];
sx q[3];
rz(-1.8545056) q[3];
sx q[3];
rz(0.32645389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3052519) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(-0.15023896) q[2];
rz(0.59208208) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2109461) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(2.3760702) q[0];
rz(-0.80134478) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(-1.9498922) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649851) q[0];
sx q[0];
rz(-2.5132127) q[0];
sx q[0];
rz(0.5324441) q[0];
rz(-pi) q[1];
rz(0.52764738) q[2];
sx q[2];
rz(-1.9867267) q[2];
sx q[2];
rz(-1.6017583) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.050151976) q[1];
sx q[1];
rz(-1.5066083) q[1];
sx q[1];
rz(-2.2540171) q[1];
rz(-0.7804773) q[3];
sx q[3];
rz(-2.7679518) q[3];
sx q[3];
rz(-2.3821156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2661065) q[2];
sx q[2];
rz(-0.85662872) q[2];
sx q[2];
rz(2.6066656) q[2];
rz(-2.5180425) q[3];
sx q[3];
rz(-1.1112735) q[3];
sx q[3];
rz(-2.258544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542338) q[0];
sx q[0];
rz(-0.41340241) q[0];
sx q[0];
rz(-2.2075388) q[0];
rz(-1.2454698) q[1];
sx q[1];
rz(-1.5555236) q[1];
sx q[1];
rz(-2.5159871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69961415) q[0];
sx q[0];
rz(-1.7480339) q[0];
sx q[0];
rz(2.3236047) q[0];
x q[1];
rz(1.6883739) q[2];
sx q[2];
rz(-2.46799) q[2];
sx q[2];
rz(2.6308311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6924507) q[1];
sx q[1];
rz(-2.3896273) q[1];
sx q[1];
rz(-0.43159952) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2295093) q[3];
sx q[3];
rz(-0.28277031) q[3];
sx q[3];
rz(-1.2706437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.049456747) q[2];
sx q[2];
rz(-2.226604) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26559386) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.884961) q[0];
rz(1.9427293) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(1.2748324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0755646) q[0];
sx q[0];
rz(-2.4504553) q[0];
sx q[0];
rz(-2.4253571) q[0];
rz(-pi) q[1];
x q[1];
rz(2.149717) q[2];
sx q[2];
rz(-0.61805586) q[2];
sx q[2];
rz(1.8628054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3209838) q[1];
sx q[1];
rz(-0.57772355) q[1];
sx q[1];
rz(-2.251365) q[1];
rz(-1.253264) q[3];
sx q[3];
rz(-1.4503363) q[3];
sx q[3];
rz(2.5250556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4510497) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(-0.96735442) q[2];
rz(1.5585772) q[3];
sx q[3];
rz(-1.5019006) q[3];
sx q[3];
rz(-2.8398452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3485182) q[0];
sx q[0];
rz(-0.59978849) q[0];
sx q[0];
rz(-3.0317958) q[0];
rz(-1.173136) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(-2.0593624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.968348) q[0];
sx q[0];
rz(-2.3254407) q[0];
sx q[0];
rz(-1.2307172) q[0];
rz(-pi) q[1];
rz(1.3668109) q[2];
sx q[2];
rz(-0.49353853) q[2];
sx q[2];
rz(-0.94369027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46690643) q[1];
sx q[1];
rz(-1.9787346) q[1];
sx q[1];
rz(0.80927421) q[1];
rz(-pi) q[2];
rz(-2.6519351) q[3];
sx q[3];
rz(-0.69935971) q[3];
sx q[3];
rz(-2.7780864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75862306) q[2];
sx q[2];
rz(-2.0489645) q[2];
sx q[2];
rz(-0.36337241) q[2];
rz(2.162497) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(0.50160971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30524224) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(-2.7889732) q[0];
rz(1.7800219) q[1];
sx q[1];
rz(-1.1905328) q[1];
sx q[1];
rz(3.000066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1175431) q[0];
sx q[0];
rz(-2.96857) q[0];
sx q[0];
rz(0.19917147) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1098712) q[2];
sx q[2];
rz(-1.0264947) q[2];
sx q[2];
rz(-1.5371295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5302116) q[1];
sx q[1];
rz(-0.57032835) q[1];
sx q[1];
rz(0.71871106) q[1];
rz(2.448672) q[3];
sx q[3];
rz(-0.77465215) q[3];
sx q[3];
rz(-1.2415394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3129468) q[2];
sx q[2];
rz(-1.4742278) q[2];
sx q[2];
rz(-1.0860156) q[2];
rz(-0.18276754) q[3];
sx q[3];
rz(-1.026231) q[3];
sx q[3];
rz(2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(0.64424789) q[0];
rz(3.0746025) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(3.0279874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73175921) q[0];
sx q[0];
rz(-1.5749965) q[0];
sx q[0];
rz(-1.3467711) q[0];
rz(-pi) q[1];
rz(-1.6729987) q[2];
sx q[2];
rz(-1.6530286) q[2];
sx q[2];
rz(0.036957994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41109243) q[1];
sx q[1];
rz(-0.20721315) q[1];
sx q[1];
rz(-2.5403008) q[1];
rz(-1.213041) q[3];
sx q[3];
rz(-1.5309257) q[3];
sx q[3];
rz(-1.2522344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99864787) q[2];
sx q[2];
rz(-1.2544268) q[2];
sx q[2];
rz(2.3764853) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(0.83711886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.9255623) q[0];
sx q[0];
rz(-0.75556527) q[0];
sx q[0];
rz(2.8788993) q[0];
rz(-0.33949159) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(1.2267494) q[2];
sx q[2];
rz(-0.94469246) q[2];
sx q[2];
rz(-2.4641987) q[2];
rz(-0.30173652) q[3];
sx q[3];
rz(-1.269125) q[3];
sx q[3];
rz(2.453809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
