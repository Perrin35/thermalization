OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(0.23101097) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(-0.45431554) q[1];
sx q[1];
rz(1.2872202) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9480913) q[0];
sx q[0];
rz(-2.5768902) q[0];
sx q[0];
rz(-2.094784) q[0];
rz(-pi) q[1];
rz(-1.8051992) q[2];
sx q[2];
rz(-1.3318828) q[2];
sx q[2];
rz(0.51942458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8227777) q[1];
sx q[1];
rz(-0.19365573) q[1];
sx q[1];
rz(0.53656399) q[1];
rz(3.0662905) q[3];
sx q[3];
rz(-1.7455532) q[3];
sx q[3];
rz(0.14698262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(1.2688961) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(0.65688175) q[0];
rz(-2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(-2.2064256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2774076) q[0];
sx q[0];
rz(-1.6011229) q[0];
sx q[0];
rz(0.10459374) q[0];
rz(-pi) q[1];
rz(-0.25721154) q[2];
sx q[2];
rz(-1.6597444) q[2];
sx q[2];
rz(-0.48438977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3352901) q[1];
sx q[1];
rz(-0.15220255) q[1];
sx q[1];
rz(-2.8782513) q[1];
rz(-1.7788497) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(1.5294242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8344581) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(1.2163986) q[2];
rz(0.52792102) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(-1.886604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(2.5566027) q[0];
rz(3.0168369) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(1.4488719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1216269) q[0];
sx q[0];
rz(-3.1375855) q[0];
sx q[0];
rz(0.83929707) q[0];
rz(-pi) q[1];
rz(-1.5572433) q[2];
sx q[2];
rz(-1.6429735) q[2];
sx q[2];
rz(-0.37301317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3347149) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(0.54912864) q[1];
x q[2];
rz(1.307906) q[3];
sx q[3];
rz(-1.2440727) q[3];
sx q[3];
rz(-0.19858352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8831732) q[2];
sx q[2];
rz(-1.0018188) q[2];
sx q[2];
rz(1.6955356) q[2];
rz(2.3840733) q[3];
sx q[3];
rz(-1.5816553) q[3];
sx q[3];
rz(2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739173) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(-1.0876592) q[0];
rz(0.37519535) q[1];
sx q[1];
rz(-1.1253858) q[1];
sx q[1];
rz(0.96022022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.023237) q[0];
sx q[0];
rz(-0.26817061) q[0];
sx q[0];
rz(1.4203181) q[0];
x q[1];
rz(-0.029898568) q[2];
sx q[2];
rz(-1.6702067) q[2];
sx q[2];
rz(-1.4292029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8254777) q[1];
sx q[1];
rz(-2.0160455) q[1];
sx q[1];
rz(-0.47834217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33444114) q[3];
sx q[3];
rz(-2.8570647) q[3];
sx q[3];
rz(-2.4872125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(1.6440294) q[2];
rz(2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(-2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(-0.60428756) q[0];
rz(-1.9970278) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(-2.7191275) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2328932) q[0];
sx q[0];
rz(-2.9541203) q[0];
sx q[0];
rz(-2.0220246) q[0];
rz(0.79549148) q[2];
sx q[2];
rz(-0.73372148) q[2];
sx q[2];
rz(2.7249641) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3389272) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(-0.1130123) q[1];
rz(2.4147968) q[3];
sx q[3];
rz(-2.6495547) q[3];
sx q[3];
rz(0.43471042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8233238) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(-2.4853415) q[2];
rz(-1.5308135) q[3];
sx q[3];
rz(-1.8717513) q[3];
sx q[3];
rz(-2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4055279) q[0];
sx q[0];
rz(-2.1298213) q[0];
sx q[0];
rz(-0.62498012) q[0];
rz(2.3896353) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(-1.5350852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4706375) q[0];
sx q[0];
rz(-2.8482901) q[0];
sx q[0];
rz(-0.17099149) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4703896) q[2];
sx q[2];
rz(-0.83121383) q[2];
sx q[2];
rz(-0.20944706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9823204) q[1];
sx q[1];
rz(-2.1814924) q[1];
sx q[1];
rz(1.8507694) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5967136) q[3];
sx q[3];
rz(-2.0996465) q[3];
sx q[3];
rz(-1.8008055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9525166) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(-2.2124186) q[2];
rz(-0.82516986) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(-1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-0.80699054) q[0];
sx q[0];
rz(-1.3744542) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-0.62320954) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40996288) q[0];
sx q[0];
rz(-1.4193168) q[0];
sx q[0];
rz(-1.234647) q[0];
x q[1];
rz(0.34131949) q[2];
sx q[2];
rz(-0.53487294) q[2];
sx q[2];
rz(0.90082263) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0125) q[1];
sx q[1];
rz(-0.18408891) q[1];
sx q[1];
rz(0.91839183) q[1];
rz(-pi) q[2];
rz(-2.3939952) q[3];
sx q[3];
rz(-1.4211402) q[3];
sx q[3];
rz(-0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0509384) q[2];
sx q[2];
rz(-2.3263558) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(0.022631571) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8959494) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(2.1242712) q[0];
rz(1.3941049) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(-0.62754935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16495569) q[0];
sx q[0];
rz(-2.4058488) q[0];
sx q[0];
rz(0.95860211) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6553788) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(-1.6135474) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4269451) q[1];
sx q[1];
rz(-1.3332187) q[1];
sx q[1];
rz(-1.3690884) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6672809) q[3];
sx q[3];
rz(-1.8758869) q[3];
sx q[3];
rz(-0.69661372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5158186) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(0.41075692) q[2];
rz(0.050203236) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(-3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.582616) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(-2.8726752) q[0];
rz(-2.8315663) q[1];
sx q[1];
rz(-2.0545484) q[1];
sx q[1];
rz(-0.1300098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271255) q[0];
sx q[0];
rz(-1.7916745) q[0];
sx q[0];
rz(1.2400024) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7902044) q[2];
sx q[2];
rz(-1.2341502) q[2];
sx q[2];
rz(-1.2996246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0669032) q[1];
sx q[1];
rz(-0.58614158) q[1];
sx q[1];
rz(-0.2939923) q[1];
x q[2];
rz(0.49874108) q[3];
sx q[3];
rz(-0.95512701) q[3];
sx q[3];
rz(-1.9174674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2076063) q[2];
sx q[2];
rz(-0.76675582) q[2];
sx q[2];
rz(0.096573528) q[2];
rz(1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0132975) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(-0.53303322) q[0];
rz(-3.10532) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(1.4520377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851327) q[0];
sx q[0];
rz(-1.0074573) q[0];
sx q[0];
rz(-2.6908532) q[0];
x q[1];
rz(2.8218837) q[2];
sx q[2];
rz(-2.8065348) q[2];
sx q[2];
rz(0.31949319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.888537) q[1];
sx q[1];
rz(-1.6007242) q[1];
sx q[1];
rz(-2.7821343) q[1];
x q[2];
rz(2.96569) q[3];
sx q[3];
rz(-0.47187343) q[3];
sx q[3];
rz(2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71296802) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(0.0326322) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3451155) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(-2.3445917) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(2.6673139) q[2];
sx q[2];
rz(-2.5534292) q[2];
sx q[2];
rz(-3.0234887) q[2];
rz(-1.6975523) q[3];
sx q[3];
rz(-2.3777665) q[3];
sx q[3];
rz(-0.39118097) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
