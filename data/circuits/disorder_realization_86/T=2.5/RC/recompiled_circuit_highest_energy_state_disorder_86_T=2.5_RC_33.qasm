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
rz(-2.4966709) q[0];
sx q[0];
rz(-2.6725197) q[0];
sx q[0];
rz(0.99035779) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(-0.85964179) q[1];
sx q[1];
rz(0.8492066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17371236) q[0];
sx q[0];
rz(-0.033941887) q[0];
sx q[0];
rz(-1.4871661) q[0];
rz(-0.87090308) q[2];
sx q[2];
rz(-1.8352301) q[2];
sx q[2];
rz(1.3178133) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1942595) q[1];
sx q[1];
rz(-2.4521356) q[1];
sx q[1];
rz(1.9973151) q[1];
x q[2];
rz(0.11470397) q[3];
sx q[3];
rz(-1.0513115) q[3];
sx q[3];
rz(1.8709559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2363756) q[2];
sx q[2];
rz(-1.6340916) q[2];
sx q[2];
rz(-2.6017453) q[2];
rz(1.5763464) q[3];
sx q[3];
rz(-0.40894517) q[3];
sx q[3];
rz(1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.0419615) q[0];
sx q[0];
rz(-2.4653682) q[0];
sx q[0];
rz(1.3270295) q[0];
rz(-0.78500336) q[1];
sx q[1];
rz(-1.7186586) q[1];
sx q[1];
rz(2.6410417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7206524) q[0];
sx q[0];
rz(-0.55188239) q[0];
sx q[0];
rz(-3.0938074) q[0];
x q[1];
rz(2.0346257) q[2];
sx q[2];
rz(-0.40574771) q[2];
sx q[2];
rz(-1.4357937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9119279) q[1];
sx q[1];
rz(-1.2357792) q[1];
sx q[1];
rz(-2.204621) q[1];
x q[2];
rz(2.331892) q[3];
sx q[3];
rz(-1.2375583) q[3];
sx q[3];
rz(2.7908195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0850204) q[2];
sx q[2];
rz(-0.73740059) q[2];
sx q[2];
rz(-2.6596587) q[2];
rz(-2.7815172) q[3];
sx q[3];
rz(-1.1432546) q[3];
sx q[3];
rz(2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29927403) q[0];
sx q[0];
rz(-2.0330918) q[0];
sx q[0];
rz(0.47766787) q[0];
rz(-2.281669) q[1];
sx q[1];
rz(-1.0483024) q[1];
sx q[1];
rz(2.8544676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52567083) q[0];
sx q[0];
rz(-1.6369372) q[0];
sx q[0];
rz(1.4445027) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7874009) q[2];
sx q[2];
rz(-1.2239211) q[2];
sx q[2];
rz(-1.5225449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9707253) q[1];
sx q[1];
rz(-1.5620462) q[1];
sx q[1];
rz(2.2544202) q[1];
rz(0.085137376) q[3];
sx q[3];
rz(-1.1907902) q[3];
sx q[3];
rz(3.0875077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3890248) q[2];
sx q[2];
rz(-1.3501046) q[2];
sx q[2];
rz(0.019850578) q[2];
rz(-0.94414532) q[3];
sx q[3];
rz(-0.70725924) q[3];
sx q[3];
rz(-0.39785644) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95553628) q[0];
sx q[0];
rz(-2.5113386) q[0];
sx q[0];
rz(1.8408884) q[0];
rz(-2.6553254) q[1];
sx q[1];
rz(-1.174289) q[1];
sx q[1];
rz(-3.1062612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0606421) q[0];
sx q[0];
rz(-1.5968906) q[0];
sx q[0];
rz(1.1423201) q[0];
x q[1];
rz(-2.6258518) q[2];
sx q[2];
rz(-2.6922142) q[2];
sx q[2];
rz(-0.55381227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9641387) q[1];
sx q[1];
rz(-1.5651917) q[1];
sx q[1];
rz(0.023375794) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0159303) q[3];
sx q[3];
rz(-1.0175704) q[3];
sx q[3];
rz(2.826626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47762927) q[2];
sx q[2];
rz(-0.56325459) q[2];
sx q[2];
rz(-2.8832054) q[2];
rz(-2.431331) q[3];
sx q[3];
rz(-0.60451549) q[3];
sx q[3];
rz(-1.5593504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.6579987) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(-2.4972231) q[0];
rz(-0.98980347) q[1];
sx q[1];
rz(-0.80392307) q[1];
sx q[1];
rz(2.03233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8744226) q[0];
sx q[0];
rz(-0.56550282) q[0];
sx q[0];
rz(-1.1995951) q[0];
x q[1];
rz(-2.0844056) q[2];
sx q[2];
rz(-1.536278) q[2];
sx q[2];
rz(-2.2396954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4586045) q[1];
sx q[1];
rz(-2.3420077) q[1];
sx q[1];
rz(1.7728642) q[1];
rz(2.3886834) q[3];
sx q[3];
rz(-1.782182) q[3];
sx q[3];
rz(1.1803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41205078) q[2];
sx q[2];
rz(-1.5770301) q[2];
sx q[2];
rz(0.61005074) q[2];
rz(-3.0610906) q[3];
sx q[3];
rz(-2.9952315) q[3];
sx q[3];
rz(-0.0065053594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215111) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(1.6625846) q[0];
rz(1.9526019) q[1];
sx q[1];
rz(-2.7956796) q[1];
sx q[1];
rz(1.4422653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449397) q[0];
sx q[0];
rz(-1.4126987) q[0];
sx q[0];
rz(2.197236) q[0];
x q[1];
rz(-1.6369319) q[2];
sx q[2];
rz(-2.1512262) q[2];
sx q[2];
rz(2.5709267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4870781) q[1];
sx q[1];
rz(-2.2527825) q[1];
sx q[1];
rz(2.712516) q[1];
rz(-pi) q[2];
rz(-0.54338065) q[3];
sx q[3];
rz(-2.5347049) q[3];
sx q[3];
rz(3.1306289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6815765) q[2];
sx q[2];
rz(-2.4676393) q[2];
sx q[2];
rz(-2.4665534) q[2];
rz(2.8071844) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(2.7736751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48208958) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(3.0159045) q[0];
rz(-2.9478759) q[1];
sx q[1];
rz(-2.2592762) q[1];
sx q[1];
rz(-1.151459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2053368) q[0];
sx q[0];
rz(-1.0837306) q[0];
sx q[0];
rz(2.4508935) q[0];
rz(2.861768) q[2];
sx q[2];
rz(-1.2132267) q[2];
sx q[2];
rz(-1.2841061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0380504) q[1];
sx q[1];
rz(-1.0443496) q[1];
sx q[1];
rz(-2.0851233) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6748333) q[3];
sx q[3];
rz(-1.950436) q[3];
sx q[3];
rz(-1.4745431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2842399) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(-2.8947158) q[2];
rz(-0.85865584) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(-0.27913678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.7849279) q[0];
sx q[0];
rz(-2.6981638) q[0];
sx q[0];
rz(2.8163633) q[0];
rz(2.6204956) q[1];
sx q[1];
rz(-2.5083713) q[1];
sx q[1];
rz(0.31347832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8756008) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(2.4317047) q[0];
x q[1];
rz(1.8584537) q[2];
sx q[2];
rz(-0.23459841) q[2];
sx q[2];
rz(1.2755659) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4904546) q[1];
sx q[1];
rz(-0.50781194) q[1];
sx q[1];
rz(2.3976624) q[1];
x q[2];
rz(1.9133592) q[3];
sx q[3];
rz(-2.5940764) q[3];
sx q[3];
rz(1.8208131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43361214) q[2];
sx q[2];
rz(-1.9233476) q[2];
sx q[2];
rz(1.8617967) q[2];
rz(-0.72569877) q[3];
sx q[3];
rz(-1.9819219) q[3];
sx q[3];
rz(-2.3316135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5028266) q[0];
sx q[0];
rz(-1.9483197) q[0];
sx q[0];
rz(-2.9916812) q[0];
rz(-1.8984849) q[1];
sx q[1];
rz(-1.5877692) q[1];
sx q[1];
rz(1.4814203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96326665) q[0];
sx q[0];
rz(-1.5218922) q[0];
sx q[0];
rz(0.045128926) q[0];
x q[1];
rz(-1.8114575) q[2];
sx q[2];
rz(-1.5207371) q[2];
sx q[2];
rz(-0.71739774) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2592377) q[1];
sx q[1];
rz(-1.9505525) q[1];
sx q[1];
rz(-1.8024551) q[1];
rz(-pi) q[2];
rz(-3.0143901) q[3];
sx q[3];
rz(-1.3281999) q[3];
sx q[3];
rz(0.85240817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3159065) q[2];
sx q[2];
rz(-1.3797727) q[2];
sx q[2];
rz(-2.1659577) q[2];
rz(1.1286831) q[3];
sx q[3];
rz(-0.55207878) q[3];
sx q[3];
rz(-2.8868207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9752556) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(2.9794203) q[0];
rz(-1.8099248) q[1];
sx q[1];
rz(-0.63260308) q[1];
sx q[1];
rz(-3.0955637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23102681) q[0];
sx q[0];
rz(-1.5648769) q[0];
sx q[0];
rz(-1.5886515) q[0];
rz(-pi) q[1];
rz(-1.1835021) q[2];
sx q[2];
rz(-2.0047054) q[2];
sx q[2];
rz(-2.6604685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6190336) q[1];
sx q[1];
rz(-1.7184034) q[1];
sx q[1];
rz(-1.9760494) q[1];
rz(0.33933731) q[3];
sx q[3];
rz(-1.4095613) q[3];
sx q[3];
rz(-1.9766903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8056246) q[2];
sx q[2];
rz(-2.7569572) q[2];
sx q[2];
rz(1.1495205) q[2];
rz(-0.51291054) q[3];
sx q[3];
rz(-1.7549113) q[3];
sx q[3];
rz(2.8368867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49160663) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(0.06123771) q[1];
sx q[1];
rz(-1.9495268) q[1];
sx q[1];
rz(-1.3757642) q[1];
rz(2.1381151) q[2];
sx q[2];
rz(-2.2879231) q[2];
sx q[2];
rz(-3.1262457) q[2];
rz(-0.41027222) q[3];
sx q[3];
rz(-2.3193852) q[3];
sx q[3];
rz(1.074765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
