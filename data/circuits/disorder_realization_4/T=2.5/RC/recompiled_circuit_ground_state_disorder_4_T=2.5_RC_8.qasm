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
rz(5.5226749) q[0];
sx q[0];
rz(10.439846) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8684611) q[0];
sx q[0];
rz(-2.1885311) q[0];
sx q[0];
rz(-2.4988089) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7741987) q[2];
sx q[2];
rz(-1.9695373) q[2];
sx q[2];
rz(-0.90786394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0933669) q[1];
sx q[1];
rz(-1.9328572) q[1];
sx q[1];
rz(1.0735379) q[1];
rz(2.8794199) q[3];
sx q[3];
rz(-2.5618563) q[3];
sx q[3];
rz(1.6593103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(0.70293054) q[2];
rz(0.85567307) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4502451) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(2.3968089) q[0];
rz(1.567747) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(2.1994798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5013393) q[0];
sx q[0];
rz(-2.3177409) q[0];
sx q[0];
rz(1.8477316) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079374119) q[2];
sx q[2];
rz(-2.3217776) q[2];
sx q[2];
rz(-0.83460966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2750863) q[1];
sx q[1];
rz(-2.3672291) q[1];
sx q[1];
rz(1.3759717) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5553988) q[3];
sx q[3];
rz(-1.6209416) q[3];
sx q[3];
rz(1.8184838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(-1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(1.7553294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11109322) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(1.4053364) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(-0.90000802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42381313) q[0];
sx q[0];
rz(-2.2909145) q[0];
sx q[0];
rz(0.49085842) q[0];
rz(-pi) q[1];
rz(2.3829557) q[2];
sx q[2];
rz(-0.73497811) q[2];
sx q[2];
rz(2.1211565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0327889) q[1];
sx q[1];
rz(-0.76251635) q[1];
sx q[1];
rz(-0.25188781) q[1];
rz(-pi) q[2];
rz(0.7358968) q[3];
sx q[3];
rz(-2.3319339) q[3];
sx q[3];
rz(-2.5084605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6806543) q[2];
sx q[2];
rz(-2.9260981) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42156521) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(-3.0928639) q[0];
rz(-0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-0.31160942) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8870114) q[0];
sx q[0];
rz(-2.0413646) q[0];
sx q[0];
rz(-0.30593095) q[0];
rz(-pi) q[1];
rz(2.7976102) q[2];
sx q[2];
rz(-1.0687573) q[2];
sx q[2];
rz(-3.0127522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8383465) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(-1.9796728) q[1];
rz(1.1295839) q[3];
sx q[3];
rz(-1.8407243) q[3];
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
rz(-0.21169855) q[2];
sx q[2];
rz(1.4908028) q[2];
rz(2.7866411) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(-1.8748913) q[0];
rz(3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(1.6273392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851885) q[0];
sx q[0];
rz(-1.8833816) q[0];
sx q[0];
rz(0.87707918) q[0];
rz(-pi) q[1];
rz(-3.051018) q[2];
sx q[2];
rz(-1.7901058) q[2];
sx q[2];
rz(-2.5136869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5942638) q[1];
sx q[1];
rz(-1.0272861) q[1];
sx q[1];
rz(0.059652358) q[1];
x q[2];
rz(-2.299304) q[3];
sx q[3];
rz(-1.5429879) q[3];
sx q[3];
rz(0.14735809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2221471) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.9661281) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(-2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670052) q[0];
sx q[0];
rz(-2.8103204) q[0];
sx q[0];
rz(0.40147716) q[0];
rz(-1.1270771) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(-2.3289767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881587) q[0];
sx q[0];
rz(-1.5666612) q[0];
sx q[0];
rz(0.93664108) q[0];
rz(1.8046384) q[2];
sx q[2];
rz(-1.3030145) q[2];
sx q[2];
rz(-0.074169548) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9795831) q[1];
sx q[1];
rz(-0.97921959) q[1];
sx q[1];
rz(-1.3237557) q[1];
x q[2];
rz(2.5162302) q[3];
sx q[3];
rz(-0.84983045) q[3];
sx q[3];
rz(1.4537077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9408985) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(1.5563439) q[2];
rz(-2.9511792) q[3];
sx q[3];
rz(-2.3868581) q[3];
sx q[3];
rz(1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82454005) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(-1.9484733) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(-2.3741123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4043264) q[0];
sx q[0];
rz(-2.4929059) q[0];
sx q[0];
rz(1.4859597) q[0];
rz(2.9323879) q[2];
sx q[2];
rz(-2.7520617) q[2];
sx q[2];
rz(2.4845633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9373405) q[1];
sx q[1];
rz(-2.2240727) q[1];
sx q[1];
rz(0.055257992) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40624491) q[3];
sx q[3];
rz(-1.5666538) q[3];
sx q[3];
rz(2.3914482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(-1.5896612) q[2];
rz(-1.6520366) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(-1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
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
rz(0.72867957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4010854) q[0];
sx q[0];
rz(-0.87941636) q[0];
sx q[0];
rz(-0.87134937) q[0];
rz(-pi) q[1];
rz(-1.7394786) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(-0.53984914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0496484) q[1];
sx q[1];
rz(-2.1990805) q[1];
sx q[1];
rz(2.6960084) q[1];
rz(0.94447926) q[3];
sx q[3];
rz(-2.2221643) q[3];
sx q[3];
rz(-0.29069549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(-0.7473839) q[2];
rz(-1.4061617) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9091699) q[0];
sx q[0];
rz(-0.94038832) q[0];
sx q[0];
rz(-0.7255834) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(0.57428378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6273514) q[0];
sx q[0];
rz(-0.19565576) q[0];
sx q[0];
rz(0.95558863) q[0];
rz(0.97522511) q[2];
sx q[2];
rz(-2.4831366) q[2];
sx q[2];
rz(-2.0582046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.060238801) q[1];
sx q[1];
rz(-1.4574086) q[1];
sx q[1];
rz(-2.3008623) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5602407) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(2.4331828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9743222) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(0.66217011) q[2];
rz(-2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(-1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(-1.3978488) q[0];
rz(-2.0715879) q[1];
sx q[1];
rz(-2.013423) q[1];
sx q[1];
rz(-1.3319344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764212) q[0];
sx q[0];
rz(-1.4201453) q[0];
sx q[0];
rz(-3.1106635) q[0];
rz(1.1564163) q[2];
sx q[2];
rz(-1.7935026) q[2];
sx q[2];
rz(-3.0321995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3146914) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(-1.3342085) q[1];
rz(0.9040971) q[3];
sx q[3];
rz(-0.81579627) q[3];
sx q[3];
rz(-2.286269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3146882) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(2.9617214) q[2];
rz(1.7449069) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(-2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.3708645) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(2.6188359) q[1];
sx q[1];
rz(-1.1154543) q[1];
sx q[1];
rz(1.9784068) q[1];
rz(-1.9942453) q[2];
sx q[2];
rz(-0.97931391) q[2];
sx q[2];
rz(-0.58936832) q[2];
rz(-0.64416364) q[3];
sx q[3];
rz(-2.5732891) q[3];
sx q[3];
rz(2.7960232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
