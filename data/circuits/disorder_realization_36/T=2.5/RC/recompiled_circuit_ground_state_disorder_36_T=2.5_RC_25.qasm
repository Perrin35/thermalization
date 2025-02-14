OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26389709) q[0];
sx q[0];
rz(-1.0385624) q[0];
sx q[0];
rz(-0.39499083) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(-1.7698987) q[1];
sx q[1];
rz(2.3514907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34694865) q[0];
sx q[0];
rz(-2.145549) q[0];
sx q[0];
rz(-1.6084421) q[0];
rz(-pi) q[1];
rz(-2.6955881) q[2];
sx q[2];
rz(-2.3141907) q[2];
sx q[2];
rz(0.80660404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1829738) q[1];
sx q[1];
rz(-1.0369025) q[1];
sx q[1];
rz(0.16927232) q[1];
rz(-pi) q[2];
rz(-0.67050855) q[3];
sx q[3];
rz(-2.5250348) q[3];
sx q[3];
rz(0.77810615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6140952) q[2];
sx q[2];
rz(-0.082753269) q[2];
sx q[2];
rz(2.5772074) q[2];
rz(0.91974059) q[3];
sx q[3];
rz(-0.65110937) q[3];
sx q[3];
rz(-1.9981599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1564002) q[0];
sx q[0];
rz(-0.99034482) q[0];
sx q[0];
rz(-1.1457957) q[0];
rz(-2.1339553) q[1];
sx q[1];
rz(-2.2538908) q[1];
sx q[1];
rz(-3.0767344) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52043693) q[0];
sx q[0];
rz(-1.1722992) q[0];
sx q[0];
rz(-0.99834672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5186132) q[2];
sx q[2];
rz(-0.88277915) q[2];
sx q[2];
rz(2.1114897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25003281) q[1];
sx q[1];
rz(-1.3866826) q[1];
sx q[1];
rz(1.5450983) q[1];
rz(-pi) q[2];
rz(-0.50527827) q[3];
sx q[3];
rz(-2.0483096) q[3];
sx q[3];
rz(1.2569515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8811983) q[2];
sx q[2];
rz(-3.03646) q[2];
sx q[2];
rz(-0.51327389) q[2];
rz(1.5567635) q[3];
sx q[3];
rz(-1.9624174) q[3];
sx q[3];
rz(1.5796327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.2960812) q[0];
sx q[0];
rz(-3.0873612) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(-0.71146479) q[1];
sx q[1];
rz(-2.3399946) q[1];
sx q[1];
rz(2.5149288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.959528) q[0];
sx q[0];
rz(-1.7668889) q[0];
sx q[0];
rz(0.0053868731) q[0];
rz(-2.9128252) q[2];
sx q[2];
rz(-1.1950777) q[2];
sx q[2];
rz(-1.8919593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16594812) q[1];
sx q[1];
rz(-0.86417809) q[1];
sx q[1];
rz(-0.23166512) q[1];
rz(-pi) q[2];
x q[2];
rz(2.061743) q[3];
sx q[3];
rz(-1.5643334) q[3];
sx q[3];
rz(0.53538509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5034379) q[2];
sx q[2];
rz(-1.902709) q[2];
sx q[2];
rz(-0.16866355) q[2];
rz(-2.9851798) q[3];
sx q[3];
rz(-2.5037239) q[3];
sx q[3];
rz(-2.8580581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6969358) q[0];
sx q[0];
rz(-1.1035656) q[0];
sx q[0];
rz(0.37100938) q[0];
rz(-1.3320529) q[1];
sx q[1];
rz(-1.5468372) q[1];
sx q[1];
rz(0.36088774) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661117) q[0];
sx q[0];
rz(-1.4752441) q[0];
sx q[0];
rz(-0.086009228) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0195186) q[2];
sx q[2];
rz(-1.586682) q[2];
sx q[2];
rz(2.0636914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8743048) q[1];
sx q[1];
rz(-2.1308769) q[1];
sx q[1];
rz(0.31927424) q[1];
x q[2];
rz(-0.64334433) q[3];
sx q[3];
rz(-2.3549101) q[3];
sx q[3];
rz(-2.1718882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0767625) q[2];
sx q[2];
rz(-2.7702489) q[2];
sx q[2];
rz(-1.4999464) q[2];
rz(1.0547981) q[3];
sx q[3];
rz(-2.0513556) q[3];
sx q[3];
rz(-1.1928026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.763279) q[0];
sx q[0];
rz(-1.9739456) q[0];
sx q[0];
rz(1.2953229) q[0];
rz(-0.12652215) q[1];
sx q[1];
rz(-0.69252068) q[1];
sx q[1];
rz(-1.893938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0737138) q[0];
sx q[0];
rz(-0.31737161) q[0];
sx q[0];
rz(-0.48010357) q[0];
x q[1];
rz(1.1791897) q[2];
sx q[2];
rz(-0.84184781) q[2];
sx q[2];
rz(-0.1192418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3646255) q[1];
sx q[1];
rz(-2.7086621) q[1];
sx q[1];
rz(2.6294699) q[1];
rz(-pi) q[2];
rz(0.18670795) q[3];
sx q[3];
rz(-1.2934341) q[3];
sx q[3];
rz(1.7565146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36914545) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(-2.6920953) q[2];
rz(-0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(1.0615758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.19675572) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(0.34360316) q[0];
rz(0.65394872) q[1];
sx q[1];
rz(-0.75075722) q[1];
sx q[1];
rz(-0.59069815) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4784929) q[0];
sx q[0];
rz(-2.0556774) q[0];
sx q[0];
rz(2.8829178) q[0];
rz(-pi) q[1];
rz(0.78738625) q[2];
sx q[2];
rz(-0.91286589) q[2];
sx q[2];
rz(1.7869365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6178083) q[1];
sx q[1];
rz(-1.7016546) q[1];
sx q[1];
rz(-0.010556825) q[1];
rz(-0.69575255) q[3];
sx q[3];
rz(-2.1348456) q[3];
sx q[3];
rz(-2.6524909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5609453) q[2];
sx q[2];
rz(-2.2717768) q[2];
sx q[2];
rz(-0.84442863) q[2];
rz(2.1010418) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(1.9677264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.7548673) q[0];
sx q[0];
rz(-0.52229053) q[0];
sx q[0];
rz(-1.0871543) q[0];
rz(2.6548751) q[1];
sx q[1];
rz(-1.4954647) q[1];
sx q[1];
rz(1.6513991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1468069) q[0];
sx q[0];
rz(-1.4658966) q[0];
sx q[0];
rz(0.86284967) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1370239) q[2];
sx q[2];
rz(-0.90551584) q[2];
sx q[2];
rz(-0.69349223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.381272) q[1];
sx q[1];
rz(-2.5650592) q[1];
sx q[1];
rz(-0.74200304) q[1];
rz(-0.20743117) q[3];
sx q[3];
rz(-1.279502) q[3];
sx q[3];
rz(-0.32793159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7955486) q[2];
sx q[2];
rz(-0.86980021) q[2];
sx q[2];
rz(2.8022433) q[2];
rz(2.7782373) q[3];
sx q[3];
rz(-0.27000579) q[3];
sx q[3];
rz(-1.5800765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5245847) q[0];
sx q[0];
rz(-2.0211077) q[0];
sx q[0];
rz(0.28550112) q[0];
rz(0.78954804) q[1];
sx q[1];
rz(-1.6098846) q[1];
sx q[1];
rz(1.550536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0302736) q[0];
sx q[0];
rz(-2.7919214) q[0];
sx q[0];
rz(0.5446148) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3197081) q[2];
sx q[2];
rz(-2.6007915) q[2];
sx q[2];
rz(-0.021440949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62666368) q[1];
sx q[1];
rz(-2.0928934) q[1];
sx q[1];
rz(1.0842647) q[1];
x q[2];
rz(0.66724369) q[3];
sx q[3];
rz(-1.5074456) q[3];
sx q[3];
rz(-2.1433307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8153814) q[2];
sx q[2];
rz(-2.9823494) q[2];
sx q[2];
rz(-2.2747269) q[2];
rz(0.73831144) q[3];
sx q[3];
rz(-1.5574995) q[3];
sx q[3];
rz(-0.90832925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0758122) q[0];
sx q[0];
rz(-2.3920238) q[0];
sx q[0];
rz(-0.69920364) q[0];
rz(0.54222822) q[1];
sx q[1];
rz(-1.7991853) q[1];
sx q[1];
rz(-0.29621616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2767773) q[0];
sx q[0];
rz(-0.84057284) q[0];
sx q[0];
rz(3.0481799) q[0];
x q[1];
rz(2.2454295) q[2];
sx q[2];
rz(-2.2632209) q[2];
sx q[2];
rz(1.2330173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6449127) q[1];
sx q[1];
rz(-1.0154004) q[1];
sx q[1];
rz(-2.4128561) q[1];
rz(-pi) q[2];
x q[2];
rz(0.086478905) q[3];
sx q[3];
rz(-1.1942549) q[3];
sx q[3];
rz(1.2264156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88973033) q[2];
sx q[2];
rz(-0.44318649) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(3.1098747) q[3];
sx q[3];
rz(-0.40516502) q[3];
sx q[3];
rz(-0.051137663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48145914) q[0];
sx q[0];
rz(-2.0429459) q[0];
sx q[0];
rz(0.81137401) q[0];
rz(-2.7157281) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(1.3070377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872914) q[0];
sx q[0];
rz(-1.5736921) q[0];
sx q[0];
rz(1.5737246) q[0];
rz(1.1063031) q[2];
sx q[2];
rz(-1.5896322) q[2];
sx q[2];
rz(2.1003758) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9237808) q[1];
sx q[1];
rz(-2.5805223) q[1];
sx q[1];
rz(-0.63528676) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47932239) q[3];
sx q[3];
rz(-0.36277825) q[3];
sx q[3];
rz(-0.58982217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9797152) q[2];
sx q[2];
rz(-2.3780509) q[2];
sx q[2];
rz(0.36583501) q[2];
rz(-0.93519917) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(2.759815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314421) q[0];
sx q[0];
rz(-1.3642385) q[0];
sx q[0];
rz(1.4862899) q[0];
rz(-1.4797795) q[1];
sx q[1];
rz(-1.9098837) q[1];
sx q[1];
rz(2.2178537) q[1];
rz(2.1187028) q[2];
sx q[2];
rz(-1.5199502) q[2];
sx q[2];
rz(0.71141031) q[2];
rz(-1.4167841) q[3];
sx q[3];
rz(-1.3091675) q[3];
sx q[3];
rz(-3.0021734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
