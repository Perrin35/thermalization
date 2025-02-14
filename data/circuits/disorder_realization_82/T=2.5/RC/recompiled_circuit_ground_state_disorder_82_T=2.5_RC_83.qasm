OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(2.3400657) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(3.0473696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202336) q[0];
sx q[0];
rz(-0.84478837) q[0];
sx q[0];
rz(2.0840706) q[0];
x q[1];
rz(0.030709312) q[2];
sx q[2];
rz(-1.2982663) q[2];
sx q[2];
rz(0.8050134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6871145) q[1];
sx q[1];
rz(-1.2839437) q[1];
sx q[1];
rz(-1.0832562) q[1];
x q[2];
rz(0.15259907) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(-1.2945557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(1.5113277) q[2];
rz(0.18621914) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(2.490624) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52361012) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(3.004916) q[0];
rz(-3.0467721) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(0.76900855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13327293) q[0];
sx q[0];
rz(-1.1964487) q[0];
sx q[0];
rz(-1.7167313) q[0];
x q[1];
rz(-1.9273754) q[2];
sx q[2];
rz(-0.90889034) q[2];
sx q[2];
rz(-1.2939158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0817854) q[1];
sx q[1];
rz(-1.4084653) q[1];
sx q[1];
rz(-0.74682208) q[1];
x q[2];
rz(1.7058825) q[3];
sx q[3];
rz(-1.6052433) q[3];
sx q[3];
rz(-2.5273539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39917055) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(0.4906022) q[3];
sx q[3];
rz(-2.6503745) q[3];
sx q[3];
rz(-2.365999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2630149) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(-1.0947134) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(-3.1375695) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72832472) q[0];
sx q[0];
rz(-0.85094977) q[0];
sx q[0];
rz(-2.9664449) q[0];
x q[1];
rz(-2.268774) q[2];
sx q[2];
rz(-0.96257639) q[2];
sx q[2];
rz(1.5814511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5044671) q[1];
sx q[1];
rz(-1.5404623) q[1];
sx q[1];
rz(-0.25719579) q[1];
rz(-pi) q[2];
rz(-0.28272776) q[3];
sx q[3];
rz(-0.73388956) q[3];
sx q[3];
rz(2.3343488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5502988) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(0.44814056) q[2];
rz(2.6552933) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(-2.015131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22911856) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(-2.4718156) q[0];
rz(-1.9122596) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(2.5965447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81789595) q[0];
sx q[0];
rz(-2.4113089) q[0];
sx q[0];
rz(0.72636168) q[0];
x q[1];
rz(2.3178905) q[2];
sx q[2];
rz(-0.81497008) q[2];
sx q[2];
rz(-2.4832249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8432359) q[1];
sx q[1];
rz(-1.3562172) q[1];
sx q[1];
rz(-0.061857865) q[1];
x q[2];
rz(-2.6757984) q[3];
sx q[3];
rz(-2.287103) q[3];
sx q[3];
rz(2.4972625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(-1.2083017) q[3];
sx q[3];
rz(-2.8390563) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912882) q[0];
sx q[0];
rz(-1.2606324) q[0];
sx q[0];
rz(0.0038797832) q[0];
rz(2.8344391) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(2.602813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2870432) q[0];
sx q[0];
rz(-1.3417305) q[0];
sx q[0];
rz(1.4255217) q[0];
x q[1];
rz(-2.712115) q[2];
sx q[2];
rz(-1.4444077) q[2];
sx q[2];
rz(2.1548558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6262349) q[1];
sx q[1];
rz(-1.9222676) q[1];
sx q[1];
rz(-0.61351794) q[1];
rz(-pi) q[2];
rz(1.1191551) q[3];
sx q[3];
rz(-0.8913826) q[3];
sx q[3];
rz(-1.7070626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6946081) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(-1.2896607) q[2];
rz(1.5223711) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30037844) q[0];
sx q[0];
rz(-0.9391681) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(-2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(-0.57714677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4704551) q[0];
sx q[0];
rz(-1.5744493) q[0];
sx q[0];
rz(1.5994344) q[0];
rz(-1.7359176) q[2];
sx q[2];
rz(-2.151474) q[2];
sx q[2];
rz(0.48027793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0967193) q[1];
sx q[1];
rz(-1.4920939) q[1];
sx q[1];
rz(-3.031879) q[1];
x q[2];
rz(0.34586819) q[3];
sx q[3];
rz(-0.78639275) q[3];
sx q[3];
rz(-2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(1.6183759) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(-2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(2.1286569) q[0];
sx q[0];
rz(-2.1067297) q[0];
sx q[0];
rz(2.3701684) q[0];
rz(2.6029288) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(-2.0753986) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8583657) q[0];
sx q[0];
rz(-1.3182501) q[0];
sx q[0];
rz(3.0622803) q[0];
x q[1];
rz(1.1723521) q[2];
sx q[2];
rz(-0.5486998) q[2];
sx q[2];
rz(-2.0946225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33246189) q[1];
sx q[1];
rz(-1.5839424) q[1];
sx q[1];
rz(-2.9244208) q[1];
rz(-pi) q[2];
rz(-1.8366424) q[3];
sx q[3];
rz(-1.2472594) q[3];
sx q[3];
rz(-2.8697447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44925877) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(2.6991357) q[2];
rz(0.19065204) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.2451179) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(-1.3099439) q[0];
rz(-0.39778057) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(1.3936874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8373972) q[0];
sx q[0];
rz(-2.1050354) q[0];
sx q[0];
rz(-2.6287352) q[0];
x q[1];
rz(-2.2289235) q[2];
sx q[2];
rz(-0.37458146) q[2];
sx q[2];
rz(-0.41824579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2893148) q[1];
sx q[1];
rz(-0.56486579) q[1];
sx q[1];
rz(-0.58128618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.769479) q[3];
sx q[3];
rz(-1.4666838) q[3];
sx q[3];
rz(-0.35752359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(3.0681211) q[2];
rz(-0.69463378) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-2.1139483) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46591127) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(2.1821816) q[0];
rz(2.361946) q[1];
sx q[1];
rz(-1.0028853) q[1];
sx q[1];
rz(-2.8693105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.773943) q[0];
sx q[0];
rz(-0.01423562) q[0];
sx q[0];
rz(-1.1971426) q[0];
rz(-pi) q[1];
rz(-1.2113038) q[2];
sx q[2];
rz(-2.4036416) q[2];
sx q[2];
rz(-2.8518845) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.694937) q[1];
sx q[1];
rz(-2.4866101) q[1];
sx q[1];
rz(0.86465624) q[1];
rz(-2.394478) q[3];
sx q[3];
rz(-1.1076339) q[3];
sx q[3];
rz(-2.2833191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73725629) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(2.9340202) q[2];
rz(0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614457) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(2.0848059) q[0];
rz(0.54569221) q[1];
sx q[1];
rz(-2.1634407) q[1];
sx q[1];
rz(-0.32870764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0595039) q[0];
sx q[0];
rz(-2.5631813) q[0];
sx q[0];
rz(-2.2057981) q[0];
x q[1];
rz(0.64735554) q[2];
sx q[2];
rz(-0.090687625) q[2];
sx q[2];
rz(0.91435223) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4663965) q[1];
sx q[1];
rz(-1.5874784) q[1];
sx q[1];
rz(-1.2287324) q[1];
rz(-pi) q[2];
rz(2.8601951) q[3];
sx q[3];
rz(-0.79165047) q[3];
sx q[3];
rz(-2.7782914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9823965) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(-2.2176149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15976739) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-2.0682206) q[2];
sx q[2];
rz(-2.2523027) q[2];
sx q[2];
rz(-0.61423617) q[2];
rz(-0.0088469395) q[3];
sx q[3];
rz(-2.736241) q[3];
sx q[3];
rz(-0.74373087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
