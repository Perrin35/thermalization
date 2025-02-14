OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57920116) q[0];
sx q[0];
rz(-0.76051036) q[0];
sx q[0];
rz(1.0150681) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88841146) q[0];
sx q[0];
rz(-1.0602573) q[0];
sx q[0];
rz(-0.84485139) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2767645) q[2];
sx q[2];
rz(-2.6061686) q[2];
sx q[2];
rz(-0.12685093) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1972292) q[1];
sx q[1];
rz(-2.5355593) q[1];
sx q[1];
rz(0.89971772) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7389033) q[3];
sx q[3];
rz(-2.1283177) q[3];
sx q[3];
rz(1.3489189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(2.4386621) q[2];
rz(0.85567307) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(-1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913476) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(-2.3968089) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(2.1994798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73977913) q[0];
sx q[0];
rz(-1.7727858) q[0];
sx q[0];
rz(0.76633472) q[0];
x q[1];
rz(3.0622185) q[2];
sx q[2];
rz(-0.81981507) q[2];
sx q[2];
rz(-2.306983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86650634) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(-1.7656209) q[1];
rz(-2.5553988) q[3];
sx q[3];
rz(-1.6209416) q[3];
sx q[3];
rz(1.8184838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88910237) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(-1.4787632) q[3];
sx q[3];
rz(-2.3507599) q[3];
sx q[3];
rz(-1.7553294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11109322) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(-1.7362562) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(-2.2415846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4858682) q[0];
sx q[0];
rz(-1.2085755) q[0];
sx q[0];
rz(0.78804228) q[0];
x q[1];
rz(0.58061231) q[2];
sx q[2];
rz(-1.0913278) q[2];
sx q[2];
rz(-1.1632077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1088037) q[1];
sx q[1];
rz(-0.76251635) q[1];
sx q[1];
rz(-2.8897048) q[1];
rz(-2.4803512) q[3];
sx q[3];
rz(-2.0783278) q[3];
sx q[3];
rz(2.762261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46093837) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(2.1920965) q[2];
rz(-0.62120581) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.1533823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42156521) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(3.0928639) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(2.8299832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64612113) q[0];
sx q[0];
rz(-2.5866383) q[0];
sx q[0];
rz(1.036219) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0428869) q[2];
sx q[2];
rz(-1.8709057) q[2];
sx q[2];
rz(-1.2712511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78098044) q[1];
sx q[1];
rz(-0.9223088) q[1];
sx q[1];
rz(2.8070252) q[1];
rz(-pi) q[2];
rz(2.0120088) q[3];
sx q[3];
rz(-1.3008683) q[3];
sx q[3];
rz(-2.2306311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17778808) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(1.6507899) q[2];
rz(-0.35495159) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(-1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(1.2667013) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(1.5142534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26402347) q[0];
sx q[0];
rz(-2.2248587) q[0];
sx q[0];
rz(2.7436849) q[0];
rz(-1.7909796) q[2];
sx q[2];
rz(-1.4823969) q[2];
sx q[2];
rz(-2.2184586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7092412) q[1];
sx q[1];
rz(-2.5951457) q[1];
sx q[1];
rz(1.6691471) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84228869) q[3];
sx q[3];
rz(-1.5986048) q[3];
sx q[3];
rz(-2.9942346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9194455) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(0.97833943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670052) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(0.40147716) q[0];
rz(-1.1270771) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(2.3289767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81174034) q[0];
sx q[0];
rz(-0.63416687) q[0];
sx q[0];
rz(-1.5638173) q[0];
rz(1.8046384) q[2];
sx q[2];
rz(-1.8385781) q[2];
sx q[2];
rz(0.074169548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9795831) q[1];
sx q[1];
rz(-2.1623731) q[1];
sx q[1];
rz(-1.817837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98315696) q[3];
sx q[3];
rz(-0.91598375) q[3];
sx q[3];
rz(2.2838796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(-2.9511792) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82454005) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(-0.12271605) q[0];
rz(1.1931194) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(0.76748031) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5106414) q[0];
sx q[0];
rz(-0.92483339) q[0];
sx q[0];
rz(0.064152282) q[0];
x q[1];
rz(-0.20920472) q[2];
sx q[2];
rz(-0.38953094) q[2];
sx q[2];
rz(-2.4845633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2042522) q[1];
sx q[1];
rz(-0.91751999) q[1];
sx q[1];
rz(-0.055257992) q[1];
x q[2];
rz(1.5753059) q[3];
sx q[3];
rz(-1.9770375) q[3];
sx q[3];
rz(0.81886983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(1.5896612) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81812304) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(2.7662011) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(-2.4129131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.323384) q[0];
sx q[0];
rz(-2.2012156) q[0];
sx q[0];
rz(0.66108836) q[0];
rz(-pi) q[1];
rz(-2.5251472) q[2];
sx q[2];
rz(-2.8551276) q[2];
sx q[2];
rz(-1.1760515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3672626) q[1];
sx q[1];
rz(-2.3891797) q[1];
sx q[1];
rz(2.1062327) q[1];
x q[2];
rz(-0.94447926) q[3];
sx q[3];
rz(-2.2221643) q[3];
sx q[3];
rz(-2.8508972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(1.4061617) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(-1.5555443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9091699) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(2.4160093) q[0];
rz(-1.3369417) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(-2.5673089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5919898) q[0];
sx q[0];
rz(-1.6832325) q[0];
sx q[0];
rz(1.7312584) q[0];
rz(2.1663675) q[2];
sx q[2];
rz(-2.4831366) q[2];
sx q[2];
rz(2.0582046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6114823) q[1];
sx q[1];
rz(-0.84646314) q[1];
sx q[1];
rz(-0.15165374) q[1];
rz(-pi) q[2];
rz(2.3575555) q[3];
sx q[3];
rz(-2.0052387) q[3];
sx q[3];
rz(1.8798352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9743222) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(-2.4794225) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(1.3242599) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26226703) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.7437438) q[0];
rz(2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(1.8096583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71026826) q[0];
sx q[0];
rz(-1.6013751) q[0];
sx q[0];
rz(1.7215183) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8990614) q[2];
sx q[2];
rz(-1.1672535) q[2];
sx q[2];
rz(-1.5833441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.012318) q[1];
sx q[1];
rz(-2.8135317) q[1];
sx q[1];
rz(-2.3543958) q[1];
rz(0.9040971) q[3];
sx q[3];
rz(-0.81579627) q[3];
sx q[3];
rz(-2.286269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3146882) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(0.17987128) q[2];
rz(-1.3966857) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(-2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(1.9942453) q[2];
sx q[2];
rz(-2.1622787) q[2];
sx q[2];
rz(2.5522243) q[2];
rz(2.497429) q[3];
sx q[3];
rz(-2.5732891) q[3];
sx q[3];
rz(2.7960232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
