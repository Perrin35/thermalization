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
rz(3.562869) q[1];
sx q[1];
rz(8.1864551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9559798) q[0];
sx q[0];
rz(-2.2816814) q[0];
sx q[0];
rz(0.86998633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9948436) q[2];
sx q[2];
rz(-1.2334261) q[2];
sx q[2];
rz(0.81126311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4290966) q[1];
sx q[1];
rz(-2.0331906) q[1];
sx q[1];
rz(-2.7346947) q[1];
x q[2];
rz(-2.8794199) q[3];
sx q[3];
rz(-2.5618563) q[3];
sx q[3];
rz(-1.6593103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51332659) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-2.4386621) q[2];
rz(2.2859196) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6913476) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(-2.3968089) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(-0.94211284) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73977913) q[0];
sx q[0];
rz(-1.3688068) q[0];
sx q[0];
rz(2.3752579) q[0];
x q[1];
rz(-1.4860542) q[2];
sx q[2];
rz(-0.7543482) q[2];
sx q[2];
rz(-2.1909383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59717885) q[1];
sx q[1];
rz(-2.326818) q[1];
sx q[1];
rz(2.9544403) q[1];
rz(-2.5553988) q[3];
sx q[3];
rz(-1.5206511) q[3];
sx q[3];
rz(1.3231089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2524903) q[2];
sx q[2];
rz(-2.1847051) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(-1.7553294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304994) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(-1.4053364) q[0];
rz(-1.7449215) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(0.90000802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8843061) q[0];
sx q[0];
rz(-0.84601706) q[0];
sx q[0];
rz(-1.0777362) q[0];
x q[1];
rz(2.3829557) q[2];
sx q[2];
rz(-2.4066145) q[2];
sx q[2];
rz(-2.1211565) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69086466) q[1];
sx q[1];
rz(-0.83800661) q[1];
sx q[1];
rz(-1.8045252) q[1];
rz(-2.4803512) q[3];
sx q[3];
rz(-1.0632648) q[3];
sx q[3];
rz(0.37933168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6806543) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(0.94949618) q[2];
rz(0.62120581) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200274) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(-3.0928639) q[0];
rz(-2.2817634) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(0.31160942) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1216868) q[2];
sx q[2];
rz(-0.60014987) q[2];
sx q[2];
rz(-0.7687591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58285634) q[1];
sx q[1];
rz(-1.3060102) q[1];
sx q[1];
rz(2.2469421) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0120088) q[3];
sx q[3];
rz(-1.3008683) q[3];
sx q[3];
rz(2.2306311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9638046) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(1.6507899) q[2];
rz(0.35495159) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.8420334) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0519003) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(1.2667013) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(-1.5142534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26402347) q[0];
sx q[0];
rz(-2.2248587) q[0];
sx q[0];
rz(-2.7436849) q[0];
rz(-pi) q[1];
rz(0.090574663) q[2];
sx q[2];
rz(-1.7901058) q[2];
sx q[2];
rz(0.62790576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5473289) q[1];
sx q[1];
rz(-2.1143066) q[1];
sx q[1];
rz(-0.059652358) q[1];
rz(-pi) q[2];
rz(-2.299304) q[3];
sx q[3];
rz(-1.5986048) q[3];
sx q[3];
rz(2.9942346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(-0.90421024) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(0.40147716) q[0];
rz(-2.0145156) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(0.81261596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75343392) q[0];
sx q[0];
rz(-1.5666612) q[0];
sx q[0];
rz(-2.2049516) q[0];
x q[1];
rz(1.8046384) q[2];
sx q[2];
rz(-1.8385781) q[2];
sx q[2];
rz(0.074169548) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1620096) q[1];
sx q[1];
rz(-0.97921959) q[1];
sx q[1];
rz(1.817837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62536247) q[3];
sx q[3];
rz(-0.84983045) q[3];
sx q[3];
rz(-1.687885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(1.5563439) q[2];
rz(-2.9511792) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3170526) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(3.0188766) q[0];
rz(1.1931194) q[1];
sx q[1];
rz(-2.2331608) q[1];
sx q[1];
rz(-0.76748031) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4043264) q[0];
sx q[0];
rz(-0.64868673) q[0];
sx q[0];
rz(1.655633) q[0];
rz(-0.38184719) q[2];
sx q[2];
rz(-1.4918461) q[2];
sx q[2];
rz(-0.71982924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7414416) q[1];
sx q[1];
rz(-1.6146683) q[1];
sx q[1];
rz(0.9167826) q[1];
rz(-1.5662868) q[3];
sx q[3];
rz(-1.9770375) q[3];
sx q[3];
rz(-2.3227228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(1.5896612) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81812304) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-2.7662011) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(-0.72867957) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4010854) q[0];
sx q[0];
rz(-0.87941636) q[0];
sx q[0];
rz(-0.87134937) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.402114) q[2];
sx q[2];
rz(-1.3381492) q[2];
sx q[2];
rz(2.6017435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77433005) q[1];
sx q[1];
rz(-2.3891797) q[1];
sx q[1];
rz(2.1062327) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3865984) q[3];
sx q[3];
rz(-2.0557311) q[3];
sx q[3];
rz(0.86673966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(-1.4061617) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9091699) q[0];
sx q[0];
rz(-0.94038832) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.888211) q[1];
sx q[1];
rz(-0.57428378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6273514) q[0];
sx q[0];
rz(-2.9459369) q[0];
sx q[0];
rz(2.186004) q[0];
rz(2.1663675) q[2];
sx q[2];
rz(-0.65845602) q[2];
sx q[2];
rz(1.083388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6114823) q[1];
sx q[1];
rz(-0.84646314) q[1];
sx q[1];
rz(2.9899389) q[1];
rz(-pi) q[2];
rz(2.5602407) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(-0.70840981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(1.3978488) q[0];
rz(-1.0700048) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(-1.3319344) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86517143) q[0];
sx q[0];
rz(-1.7214473) q[0];
sx q[0];
rz(0.030929203) q[0];
rz(-pi) q[1];
rz(-2.0831624) q[2];
sx q[2];
rz(-2.67423) q[2];
sx q[2];
rz(2.1455554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9402658) q[1];
sx q[1];
rz(-1.801071) q[1];
sx q[1];
rz(0.23576945) q[1];
rz(-pi) q[2];
rz(0.58140786) q[3];
sx q[3];
rz(-2.1801342) q[3];
sx q[3];
rz(-1.7096303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3146882) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(0.17987128) q[2];
rz(-1.7449069) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(1.1473473) q[2];
sx q[2];
rz(-0.97931391) q[2];
sx q[2];
rz(-0.58936832) q[2];
rz(-1.9369851) q[3];
sx q[3];
rz(-1.1259176) q[3];
sx q[3];
rz(2.068145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
