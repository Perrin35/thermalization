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
rz(-1.9395186) q[0];
sx q[0];
rz(-1.0459067) q[0];
sx q[0];
rz(-2.3910971) q[0];
rz(0.22076386) q[1];
sx q[1];
rz(-2.0808487) q[1];
sx q[1];
rz(-1.1836675) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7322217) q[0];
sx q[0];
rz(-0.82287649) q[0];
sx q[0];
rz(-2.5636682) q[0];
rz(-pi) q[1];
rz(-0.92249845) q[2];
sx q[2];
rz(-2.9335409) q[2];
sx q[2];
rz(-1.3626984) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69931251) q[1];
sx q[1];
rz(-0.89095014) q[1];
sx q[1];
rz(-1.5857497) q[1];
x q[2];
rz(1.8791844) q[3];
sx q[3];
rz(-1.7636136) q[3];
sx q[3];
rz(-2.8578025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.815328) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(0.092078837) q[2];
rz(2.0103256) q[3];
sx q[3];
rz(-1.0695894) q[3];
sx q[3];
rz(0.85589516) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49694127) q[0];
sx q[0];
rz(-0.48321378) q[0];
sx q[0];
rz(2.605873) q[0];
rz(3.1247395) q[1];
sx q[1];
rz(-2.2622175) q[1];
sx q[1];
rz(-0.0171612) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494228) q[0];
sx q[0];
rz(-2.0606406) q[0];
sx q[0];
rz(-2.7754098) q[0];
x q[1];
rz(1.503741) q[2];
sx q[2];
rz(-2.0366327) q[2];
sx q[2];
rz(-0.28595823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70046798) q[1];
sx q[1];
rz(-1.4135612) q[1];
sx q[1];
rz(-2.9067405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6616749) q[3];
sx q[3];
rz(-1.7044997) q[3];
sx q[3];
rz(2.9667678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0422334) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(-0.24822203) q[2];
rz(-2.21375) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(-0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27702734) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(-0.16265854) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-0.22480741) q[1];
sx q[1];
rz(0.83388296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57788173) q[0];
sx q[0];
rz(-2.1378008) q[0];
sx q[0];
rz(-1.4253229) q[0];
x q[1];
rz(2.7835186) q[2];
sx q[2];
rz(-0.34293567) q[2];
sx q[2];
rz(-1.6274522) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0625588) q[1];
sx q[1];
rz(-1.9818881) q[1];
sx q[1];
rz(-2.316019) q[1];
x q[2];
rz(2.7008204) q[3];
sx q[3];
rz(-1.6218054) q[3];
sx q[3];
rz(0.15404242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0791846) q[2];
sx q[2];
rz(-1.1391613) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(1.0583813) q[3];
sx q[3];
rz(-0.94401413) q[3];
sx q[3];
rz(0.03037608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-2.3069032) q[0];
sx q[0];
rz(1.7133065) q[0];
rz(0.76438534) q[1];
sx q[1];
rz(-1.9995707) q[1];
sx q[1];
rz(1.3216602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40250889) q[0];
sx q[0];
rz(-2.5989977) q[0];
sx q[0];
rz(0.71904166) q[0];
x q[1];
rz(-1.4218421) q[2];
sx q[2];
rz(-0.33518727) q[2];
sx q[2];
rz(0.4544979) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5668754) q[1];
sx q[1];
rz(-1.2657026) q[1];
sx q[1];
rz(1.2724691) q[1];
rz(-pi) q[2];
rz(1.0760078) q[3];
sx q[3];
rz(-0.77829276) q[3];
sx q[3];
rz(-0.28631223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4798212) q[2];
sx q[2];
rz(-1.7281374) q[2];
sx q[2];
rz(-2.856355) q[2];
rz(-1.8789004) q[3];
sx q[3];
rz(-1.9243536) q[3];
sx q[3];
rz(-1.7181905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30094639) q[0];
sx q[0];
rz(-2.3941289) q[0];
sx q[0];
rz(-2.2593011) q[0];
rz(0.67059416) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(0.98463279) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0142954) q[0];
sx q[0];
rz(-1.6335575) q[0];
sx q[0];
rz(1.602229) q[0];
rz(-1.6419471) q[2];
sx q[2];
rz(-0.54279581) q[2];
sx q[2];
rz(2.0031479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2217) q[1];
sx q[1];
rz(-1.0439149) q[1];
sx q[1];
rz(-2.5008898) q[1];
x q[2];
rz(2.4299942) q[3];
sx q[3];
rz(-0.67005537) q[3];
sx q[3];
rz(-1.3942476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8661477) q[2];
sx q[2];
rz(-0.53469849) q[2];
sx q[2];
rz(2.5246485) q[2];
rz(-0.75211891) q[3];
sx q[3];
rz(-1.813846) q[3];
sx q[3];
rz(0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1521456) q[0];
sx q[0];
rz(-2.4247657) q[0];
sx q[0];
rz(0.77051198) q[0];
rz(2.41467) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(2.4368584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7186747) q[0];
sx q[0];
rz(-0.98121907) q[0];
sx q[0];
rz(-2.4797863) q[0];
rz(-2.1628468) q[2];
sx q[2];
rz(-0.60475588) q[2];
sx q[2];
rz(0.7515242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7721036) q[1];
sx q[1];
rz(-1.2210969) q[1];
sx q[1];
rz(1.5353126) q[1];
x q[2];
rz(-2.8877657) q[3];
sx q[3];
rz(-2.1269264) q[3];
sx q[3];
rz(-2.0230011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3720234) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(-0.43323576) q[2];
rz(-0.20089928) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(2.6065629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83605003) q[0];
sx q[0];
rz(-1.1490281) q[0];
sx q[0];
rz(2.7591163) q[0];
rz(-1.0410694) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(0.44630757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8167717) q[0];
sx q[0];
rz(-1.3123625) q[0];
sx q[0];
rz(1.6245317) q[0];
rz(2.7591797) q[2];
sx q[2];
rz(-1.7837886) q[2];
sx q[2];
rz(-1.8960003) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81627264) q[1];
sx q[1];
rz(-2.153646) q[1];
sx q[1];
rz(0.42709777) q[1];
rz(-pi) q[2];
rz(-2.2931523) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(0.27001303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81991601) q[2];
sx q[2];
rz(-1.1175464) q[2];
sx q[2];
rz(2.4134911) q[2];
rz(-0.32650945) q[3];
sx q[3];
rz(-1.510334) q[3];
sx q[3];
rz(-0.85730332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(-1.8439199) q[0];
rz(-2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(-0.40850684) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9169642) q[0];
sx q[0];
rz(-3.1273807) q[0];
sx q[0];
rz(-1.5111708) q[0];
rz(0.82168545) q[2];
sx q[2];
rz(-2.5251881) q[2];
sx q[2];
rz(-2.6954755) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6031987) q[1];
sx q[1];
rz(-1.6848752) q[1];
sx q[1];
rz(-1.7054547) q[1];
rz(-0.37398963) q[3];
sx q[3];
rz(-1.1701823) q[3];
sx q[3];
rz(3.0246322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11477509) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(-2.3180023) q[2];
rz(-0.94281998) q[3];
sx q[3];
rz(-1.4225682) q[3];
sx q[3];
rz(1.6387117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59458643) q[0];
sx q[0];
rz(-2.6972045) q[0];
sx q[0];
rz(1.8894926) q[0];
rz(1.3638672) q[1];
sx q[1];
rz(-1.5003279) q[1];
sx q[1];
rz(-1.8069256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8060243) q[0];
sx q[0];
rz(-0.3055521) q[0];
sx q[0];
rz(-0.1211042) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7450856) q[2];
sx q[2];
rz(-0.93610033) q[2];
sx q[2];
rz(-1.8746992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.081055911) q[1];
sx q[1];
rz(-1.6439652) q[1];
sx q[1];
rz(-0.45159486) q[1];
rz(-1.1030508) q[3];
sx q[3];
rz(-1.333101) q[3];
sx q[3];
rz(-2.2464744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021726457) q[2];
sx q[2];
rz(-0.29325565) q[2];
sx q[2];
rz(2.8450656) q[2];
rz(-1.6622539) q[3];
sx q[3];
rz(-1.6300423) q[3];
sx q[3];
rz(2.9843946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.5088365) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(-2.346709) q[0];
rz(0.55566135) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(-1.6937675) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7894365) q[0];
sx q[0];
rz(-1.1934501) q[0];
sx q[0];
rz(2.3579602) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5083639) q[2];
sx q[2];
rz(-1.0737906) q[2];
sx q[2];
rz(1.4344429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8965079) q[1];
sx q[1];
rz(-1.0963529) q[1];
sx q[1];
rz(-2.9175379) q[1];
rz(-pi) q[2];
rz(0.85526222) q[3];
sx q[3];
rz(-2.4088833) q[3];
sx q[3];
rz(-0.39855188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8453703) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(1.4307865) q[2];
rz(-0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(-2.9770765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13567781) q[0];
sx q[0];
rz(-0.68886859) q[0];
sx q[0];
rz(-1.9266358) q[0];
rz(-1.4665435) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(2.3382414) q[2];
sx q[2];
rz(-2.4016082) q[2];
sx q[2];
rz(-0.27863816) q[2];
rz(-1.0115276) q[3];
sx q[3];
rz(-0.20023919) q[3];
sx q[3];
rz(-1.7773624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
