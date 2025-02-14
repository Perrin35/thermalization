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
rz(2.1030302) q[0];
sx q[0];
rz(9.8197688) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(-1.7698987) q[1];
sx q[1];
rz(-0.79010195) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2443197) q[0];
sx q[0];
rz(-1.5392015) q[0];
sx q[0];
rz(0.57507609) q[0];
x q[1];
rz(2.3656125) q[2];
sx q[2];
rz(-1.2476412) q[2];
sx q[2];
rz(2.6903642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2826281) q[1];
sx q[1];
rz(-2.583995) q[1];
sx q[1];
rz(-1.8484115) q[1];
rz(2.6346968) q[3];
sx q[3];
rz(-1.9383176) q[3];
sx q[3];
rz(-1.3669922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52749741) q[2];
sx q[2];
rz(-3.0588394) q[2];
sx q[2];
rz(0.56438524) q[2];
rz(-2.2218521) q[3];
sx q[3];
rz(-2.4904833) q[3];
sx q[3];
rz(1.9981599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98519242) q[0];
sx q[0];
rz(-0.99034482) q[0];
sx q[0];
rz(-1.9957969) q[0];
rz(2.1339553) q[1];
sx q[1];
rz(-2.2538908) q[1];
sx q[1];
rz(3.0767344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52043693) q[0];
sx q[0];
rz(-1.1722992) q[0];
sx q[0];
rz(2.1432459) q[0];
x q[1];
rz(-2.5186132) q[2];
sx q[2];
rz(-2.2588135) q[2];
sx q[2];
rz(1.0301029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7520719) q[1];
sx q[1];
rz(-0.18587843) q[1];
sx q[1];
rz(-3.0044772) q[1];
x q[2];
rz(-0.50527827) q[3];
sx q[3];
rz(-1.0932831) q[3];
sx q[3];
rz(-1.2569515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2603944) q[2];
sx q[2];
rz(-3.03646) q[2];
sx q[2];
rz(-2.6283188) q[2];
rz(-1.5567635) q[3];
sx q[3];
rz(-1.1791752) q[3];
sx q[3];
rz(1.5796327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2960812) q[0];
sx q[0];
rz(-0.054231461) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(0.71146479) q[1];
sx q[1];
rz(-2.3399946) q[1];
sx q[1];
rz(-2.5149288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7539106) q[0];
sx q[0];
rz(-1.5655127) q[0];
sx q[0];
rz(1.374701) q[0];
x q[1];
rz(-2.0925631) q[2];
sx q[2];
rz(-0.43704068) q[2];
sx q[2];
rz(-2.4573987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8887254) q[1];
sx q[1];
rz(-1.3952726) q[1];
sx q[1];
rz(2.2908126) q[1];
rz(3.1342642) q[3];
sx q[3];
rz(-2.0617319) q[3];
sx q[3];
rz(-2.1027264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5034379) q[2];
sx q[2];
rz(-1.902709) q[2];
sx q[2];
rz(2.9729291) q[2];
rz(-0.1564129) q[3];
sx q[3];
rz(-0.63786879) q[3];
sx q[3];
rz(0.28353459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6969358) q[0];
sx q[0];
rz(-1.1035656) q[0];
sx q[0];
rz(-0.37100938) q[0];
rz(-1.8095398) q[1];
sx q[1];
rz(-1.5468372) q[1];
sx q[1];
rz(-0.36088774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7545032) q[0];
sx q[0];
rz(-1.4851804) q[0];
sx q[0];
rz(-1.4748917) q[0];
x q[1];
rz(3.0118594) q[2];
sx q[2];
rz(-0.12309821) q[2];
sx q[2];
rz(-2.77746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8743048) q[1];
sx q[1];
rz(-2.1308769) q[1];
sx q[1];
rz(-0.31927424) q[1];
rz(-pi) q[2];
rz(2.4982483) q[3];
sx q[3];
rz(-0.78668252) q[3];
sx q[3];
rz(-0.96970448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0648301) q[2];
sx q[2];
rz(-0.37134376) q[2];
sx q[2];
rz(-1.6416462) q[2];
rz(-1.0547981) q[3];
sx q[3];
rz(-2.0513556) q[3];
sx q[3];
rz(-1.9487901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.763279) q[0];
sx q[0];
rz(-1.9739456) q[0];
sx q[0];
rz(-1.2953229) q[0];
rz(0.12652215) q[1];
sx q[1];
rz(-2.449072) q[1];
sx q[1];
rz(-1.893938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0737138) q[0];
sx q[0];
rz(-0.31737161) q[0];
sx q[0];
rz(0.48010357) q[0];
x q[1];
rz(-1.1791897) q[2];
sx q[2];
rz(-2.2997448) q[2];
sx q[2];
rz(-0.1192418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3314455) q[1];
sx q[1];
rz(-1.9451912) q[1];
sx q[1];
rz(1.7935171) q[1];
x q[2];
rz(1.8528102) q[3];
sx q[3];
rz(-1.3913034) q[3];
sx q[3];
rz(-0.13403758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36914545) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(-2.6920953) q[2];
rz(0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(-1.0615758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19675572) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(0.34360316) q[0];
rz(-2.4876439) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(0.59069815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1794248) q[0];
sx q[0];
rz(-0.54467595) q[0];
sx q[0];
rz(-1.1187798) q[0];
rz(-pi) q[1];
rz(-0.82876708) q[2];
sx q[2];
rz(-0.97835079) q[2];
sx q[2];
rz(2.3784021) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52378435) q[1];
sx q[1];
rz(-1.439938) q[1];
sx q[1];
rz(0.010556825) q[1];
rz(-pi) q[2];
rz(-0.88149397) q[3];
sx q[3];
rz(-2.1432264) q[3];
sx q[3];
rz(1.6400169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5609453) q[2];
sx q[2];
rz(-0.86981589) q[2];
sx q[2];
rz(2.297164) q[2];
rz(2.1010418) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(-1.1738663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867253) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(-2.0544384) q[0];
rz(2.6548751) q[1];
sx q[1];
rz(-1.4954647) q[1];
sx q[1];
rz(1.6513991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1468069) q[0];
sx q[0];
rz(-1.4658966) q[0];
sx q[0];
rz(2.278743) q[0];
rz(-pi) q[1];
rz(-2.1370239) q[2];
sx q[2];
rz(-0.90551584) q[2];
sx q[2];
rz(-0.69349223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9864038) q[1];
sx q[1];
rz(-1.1935368) q[1];
sx q[1];
rz(2.6946486) q[1];
rz(-pi) q[2];
rz(2.1726726) q[3];
sx q[3];
rz(-2.7857097) q[3];
sx q[3];
rz(2.8371135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7955486) q[2];
sx q[2];
rz(-2.2717924) q[2];
sx q[2];
rz(-0.33934936) q[2];
rz(0.36335534) q[3];
sx q[3];
rz(-0.27000579) q[3];
sx q[3];
rz(1.5800765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5245847) q[0];
sx q[0];
rz(-2.0211077) q[0];
sx q[0];
rz(2.8560915) q[0];
rz(-0.78954804) q[1];
sx q[1];
rz(-1.6098846) q[1];
sx q[1];
rz(-1.550536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94211468) q[0];
sx q[0];
rz(-1.7492332) q[0];
sx q[0];
rz(-2.8392544) q[0];
rz(-1.8218846) q[2];
sx q[2];
rz(-0.54080117) q[2];
sx q[2];
rz(0.021440949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68622491) q[1];
sx q[1];
rz(-1.9880725) q[1];
sx q[1];
rz(-0.57699202) q[1];
rz(-3.0394394) q[3];
sx q[3];
rz(-0.66978598) q[3];
sx q[3];
rz(-0.49234351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32621128) q[2];
sx q[2];
rz(-0.15924328) q[2];
sx q[2];
rz(2.2747269) q[2];
rz(-2.4032812) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(-2.2332634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0758122) q[0];
sx q[0];
rz(-2.3920238) q[0];
sx q[0];
rz(-0.69920364) q[0];
rz(2.5993644) q[1];
sx q[1];
rz(-1.3424073) q[1];
sx q[1];
rz(2.8453765) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86481536) q[0];
sx q[0];
rz(-2.3010198) q[0];
sx q[0];
rz(0.093412799) q[0];
x q[1];
rz(0.64546236) q[2];
sx q[2];
rz(-0.92593595) q[2];
sx q[2];
rz(-1.0114111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6275639) q[1];
sx q[1];
rz(-0.96935287) q[1];
sx q[1];
rz(0.87694962) q[1];
rz(-pi) q[2];
rz(3.0551137) q[3];
sx q[3];
rz(-1.1942549) q[3];
sx q[3];
rz(-1.2264156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2518623) q[2];
sx q[2];
rz(-0.44318649) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(-3.1098747) q[3];
sx q[3];
rz(-2.7364276) q[3];
sx q[3];
rz(3.090455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.6601335) q[0];
sx q[0];
rz(-1.0986468) q[0];
sx q[0];
rz(-2.3302186) q[0];
rz(2.7157281) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(-1.3070377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872914) q[0];
sx q[0];
rz(-1.5679006) q[0];
sx q[0];
rz(1.5737246) q[0];
rz(-pi) q[1];
rz(0.021067384) q[2];
sx q[2];
rz(-1.1063919) q[2];
sx q[2];
rz(-0.5201425) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9237808) q[1];
sx q[1];
rz(-2.5805223) q[1];
sx q[1];
rz(2.5063059) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32486954) q[3];
sx q[3];
rz(-1.7351955) q[3];
sx q[3];
rz(-1.4332958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1618774) q[2];
sx q[2];
rz(-0.76354176) q[2];
sx q[2];
rz(-0.36583501) q[2];
rz(-2.2063935) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(0.38177761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314421) q[0];
sx q[0];
rz(-1.3642385) q[0];
sx q[0];
rz(1.4862899) q[0];
rz(-1.6618132) q[1];
sx q[1];
rz(-1.231709) q[1];
sx q[1];
rz(-0.92373893) q[1];
rz(2.1187028) q[2];
sx q[2];
rz(-1.5199502) q[2];
sx q[2];
rz(0.71141031) q[2];
rz(1.4167841) q[3];
sx q[3];
rz(-1.8324251) q[3];
sx q[3];
rz(0.13941924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
