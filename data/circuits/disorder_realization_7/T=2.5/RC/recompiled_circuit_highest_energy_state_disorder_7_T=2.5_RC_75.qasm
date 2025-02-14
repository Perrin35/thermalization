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
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(-1.4359441) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4592847) q[0];
sx q[0];
rz(-2.3536682) q[0];
sx q[0];
rz(0.81972229) q[0];
rz(-1.1167657) q[2];
sx q[2];
rz(-0.28677702) q[2];
sx q[2];
rz(1.3146626) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85992009) q[1];
sx q[1];
rz(-1.8919196) q[1];
sx q[1];
rz(-0.74076498) q[1];
rz(-0.30835704) q[3];
sx q[3];
rz(-1.7614109) q[3];
sx q[3];
rz(-2.7558143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5982738) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(0.16779009) q[2];
rz(2.0134036) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(1.1892085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0447277) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(-1.0750394) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(1.8278106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2883817) q[0];
sx q[0];
rz(-2.7409005) q[0];
sx q[0];
rz(-2.543422) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4046217) q[2];
sx q[2];
rz(-0.24925512) q[2];
sx q[2];
rz(2.1971306) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35006501) q[1];
sx q[1];
rz(-2.413837) q[1];
sx q[1];
rz(1.4552874) q[1];
x q[2];
rz(1.9017392) q[3];
sx q[3];
rz(-1.6972491) q[3];
sx q[3];
rz(-3.0883873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4904867) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(2.7172078) q[2];
rz(2.4091447) q[3];
sx q[3];
rz(-1.4832352) q[3];
sx q[3];
rz(1.6781767) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80377793) q[0];
sx q[0];
rz(-0.71500272) q[0];
sx q[0];
rz(-0.5027698) q[0];
rz(-0.94364199) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(0.81519333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40476) q[0];
sx q[0];
rz(-1.6175999) q[0];
sx q[0];
rz(2.9535049) q[0];
rz(-pi) q[1];
rz(-1.6770467) q[2];
sx q[2];
rz(-1.4521404) q[2];
sx q[2];
rz(3.075171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1285308) q[1];
sx q[1];
rz(-1.3088041) q[1];
sx q[1];
rz(-2.872741) q[1];
rz(0.87343862) q[3];
sx q[3];
rz(-1.1213574) q[3];
sx q[3];
rz(-0.12974236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.076685585) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(-2.0150851) q[2];
rz(-1.7846151) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55946881) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(2.7816787) q[0];
rz(-2.9008046) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(2.7684033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561962) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(-1.7147934) q[0];
x q[1];
rz(0.56314205) q[2];
sx q[2];
rz(-1.7703319) q[2];
sx q[2];
rz(-2.2468001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9149941) q[1];
sx q[1];
rz(-1.0553331) q[1];
sx q[1];
rz(-1.9644009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3803102) q[3];
sx q[3];
rz(-1.130874) q[3];
sx q[3];
rz(-0.89118496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.084426247) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(2.2331494) q[2];
rz(-2.0716136) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(0.98328868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1835566) q[0];
sx q[0];
rz(-0.89360845) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(-2.4519582) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(0.77401179) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000074) q[0];
sx q[0];
rz(-0.98709269) q[0];
sx q[0];
rz(-2.1271749) q[0];
rz(-2.3460057) q[2];
sx q[2];
rz(-1.5485816) q[2];
sx q[2];
rz(1.3327679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23394984) q[1];
sx q[1];
rz(-1.3798668) q[1];
sx q[1];
rz(1.1348866) q[1];
rz(-3.0930107) q[3];
sx q[3];
rz(-0.54837275) q[3];
sx q[3];
rz(2.2480132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8959877) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(2.4020933) q[2];
rz(-1.6351522) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6226115) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(0.13825026) q[1];
sx q[1];
rz(-1.6743276) q[1];
sx q[1];
rz(-1.5672055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.137321) q[0];
sx q[0];
rz(-1.204365) q[0];
sx q[0];
rz(-1.2211114) q[0];
x q[1];
rz(1.862941) q[2];
sx q[2];
rz(-1.1730447) q[2];
sx q[2];
rz(-1.6145371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3688599) q[1];
sx q[1];
rz(-0.28801685) q[1];
sx q[1];
rz(2.6898328) q[1];
rz(-pi) q[2];
rz(2.592347) q[3];
sx q[3];
rz(-1.7399144) q[3];
sx q[3];
rz(-2.0813519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77341998) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(2.0424021) q[2];
rz(2.1558971) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(2.3327935) q[0];
rz(-2.47593) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(-0.98141247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.736826) q[0];
sx q[0];
rz(-0.5739218) q[0];
sx q[0];
rz(1.6762275) q[0];
rz(-1.8265884) q[2];
sx q[2];
rz(-1.0793387) q[2];
sx q[2];
rz(1.0813528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1011103) q[1];
sx q[1];
rz(-2.6553934) q[1];
sx q[1];
rz(-1.2771439) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0043169) q[3];
sx q[3];
rz(-2.2386754) q[3];
sx q[3];
rz(1.560451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0770646) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(0.36901739) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(-1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3082665) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(-0.33018026) q[0];
rz(-0.45686832) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(2.5433069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56494026) q[0];
sx q[0];
rz(-0.9441174) q[0];
sx q[0];
rz(-2.1184854) q[0];
x q[1];
rz(3.0706579) q[2];
sx q[2];
rz(-1.1026376) q[2];
sx q[2];
rz(2.2092961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48409325) q[1];
sx q[1];
rz(-1.3417021) q[1];
sx q[1];
rz(-0.34219663) q[1];
rz(-pi) q[2];
rz(2.4870858) q[3];
sx q[3];
rz(-1.084436) q[3];
sx q[3];
rz(0.24412046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8699708) q[2];
sx q[2];
rz(-0.87646708) q[2];
sx q[2];
rz(0.064621933) q[2];
rz(-1.4501075) q[3];
sx q[3];
rz(-0.40018574) q[3];
sx q[3];
rz(1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794063) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(-2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(2.1102139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8143896) q[0];
sx q[0];
rz(-1.4966244) q[0];
sx q[0];
rz(0.73930862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7613369) q[2];
sx q[2];
rz(-2.8275942) q[2];
sx q[2];
rz(1.2316201) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4881812) q[1];
sx q[1];
rz(-2.0926884) q[1];
sx q[1];
rz(0.33409358) q[1];
rz(-1.6691339) q[3];
sx q[3];
rz(-1.0957484) q[3];
sx q[3];
rz(2.6779384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7445659) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(-0.087372027) q[2];
rz(-1.2157863) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66607296) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(-1.1391621) q[0];
rz(0.49925223) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(1.3334691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718057) q[0];
sx q[0];
rz(-0.13955586) q[0];
sx q[0];
rz(-2.4371171) q[0];
rz(-pi) q[1];
rz(1.8257481) q[2];
sx q[2];
rz(-1.5226018) q[2];
sx q[2];
rz(0.15010897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1386316) q[1];
sx q[1];
rz(-1.3182782) q[1];
sx q[1];
rz(-1.7934402) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74204294) q[3];
sx q[3];
rz(-2.2931406) q[3];
sx q[3];
rz(-2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(-2.8583543) q[2];
rz(1.9404274) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(3.037187) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(-1.7684284) q[2];
sx q[2];
rz(-1.0734954) q[2];
sx q[2];
rz(-1.3263477) q[2];
rz(1.9857882) q[3];
sx q[3];
rz(-2.1963051) q[3];
sx q[3];
rz(2.5122234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
