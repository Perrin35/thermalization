OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(-1.0150681) q[0];
rz(1.9782344) q[1];
sx q[1];
rz(-0.42127633) q[1];
sx q[1];
rz(-1.9032698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8684611) q[0];
sx q[0];
rz(-0.95306153) q[0];
sx q[0];
rz(-0.6427838) q[0];
rz(-pi) q[1];
rz(-1.9948436) q[2];
sx q[2];
rz(-1.2334261) q[2];
sx q[2];
rz(-0.81126311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1972292) q[1];
sx q[1];
rz(-2.5355593) q[1];
sx q[1];
rz(-0.89971772) q[1];
x q[2];
rz(2.8794199) q[3];
sx q[3];
rz(-2.5618563) q[3];
sx q[3];
rz(-1.4822823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-2.4386621) q[2];
rz(0.85567307) q[3];
sx q[3];
rz(-2.296505) q[3];
sx q[3];
rz(1.2969016) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4502451) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(-2.3968089) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(-0.94211284) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73977913) q[0];
sx q[0];
rz(-1.7727858) q[0];
sx q[0];
rz(-2.3752579) q[0];
rz(-1.4860542) q[2];
sx q[2];
rz(-0.7543482) q[2];
sx q[2];
rz(0.95065439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86650634) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(-1.3759717) q[1];
rz(-pi) q[2];
rz(-3.0511176) q[3];
sx q[3];
rz(-2.5535085) q[3];
sx q[3];
rz(-2.9692441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2524903) q[2];
sx q[2];
rz(-2.1847051) q[2];
sx q[2];
rz(1.0955742) q[2];
rz(-1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(-1.3862632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(3.0304994) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(1.4053364) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(0.90000802) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6557245) q[0];
sx q[0];
rz(-1.9330171) q[0];
sx q[0];
rz(-0.78804228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58061231) q[2];
sx q[2];
rz(-2.0502649) q[2];
sx q[2];
rz(1.978385) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1088037) q[1];
sx q[1];
rz(-0.76251635) q[1];
sx q[1];
rz(-2.8897048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4803512) q[3];
sx q[3];
rz(-2.0783278) q[3];
sx q[3];
rz(-2.762261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6806543) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(2.1920965) q[2];
rz(-2.5203868) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(1.1533823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42156521) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(-0.048728745) q[0];
rz(2.2817634) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-0.31160942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4954715) q[0];
sx q[0];
rz(-0.55495431) q[0];
sx q[0];
rz(2.1053736) q[0];
rz(-2.1216868) q[2];
sx q[2];
rz(-2.5414428) q[2];
sx q[2];
rz(-2.3728336) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8383465) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.1619199) q[1];
x q[2];
rz(-2.8446557) q[3];
sx q[3];
rz(-1.9949759) q[3];
sx q[3];
rz(-0.78510982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9638046) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(-1.4908028) q[2];
rz(2.7866411) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(1.8420334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896924) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(1.2667013) q[0];
rz(-0.1162687) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(1.6273392) q[1];
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
rz(-2.5136869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0872495) q[1];
sx q[1];
rz(-1.5197481) q[1];
sx q[1];
rz(-2.1150949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.037260696) q[3];
sx q[3];
rz(-2.2989591) q[3];
sx q[3];
rz(1.3986349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(-1.9661281) q[2];
rz(-0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(0.97833943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(-2.7401155) q[0];
rz(-2.0145156) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(-0.81261596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82040374) q[0];
sx q[0];
rz(-0.93664737) q[0];
sx q[0];
rz(3.1364596) q[0];
rz(-pi) q[1];
x q[1];
rz(2.866686) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(-1.7079086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26906313) q[1];
sx q[1];
rz(-1.7751964) q[1];
sx q[1];
rz(0.60592954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98315696) q[3];
sx q[3];
rz(-0.91598375) q[3];
sx q[3];
rz(-0.85771307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9408985) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(1.5563439) q[2];
rz(0.19041348) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3170526) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(1.9484733) q[1];
sx q[1];
rz(-2.2331608) q[1];
sx q[1];
rz(-2.3741123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5106414) q[0];
sx q[0];
rz(-0.92483339) q[0];
sx q[0];
rz(-0.064152282) q[0];
rz(-pi) q[1];
rz(0.20920472) q[2];
sx q[2];
rz(-0.38953094) q[2];
sx q[2];
rz(2.4845633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7414416) q[1];
sx q[1];
rz(-1.6146683) q[1];
sx q[1];
rz(-2.2248101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7353477) q[3];
sx q[3];
rz(-1.5749388) q[3];
sx q[3];
rz(-2.3914482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3976589) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(1.5519315) q[2];
rz(1.6520366) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(1.1154729) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81812304) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-2.7662011) q[0];
rz(-2.2581532) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7405072) q[0];
sx q[0];
rz(-2.2621763) q[0];
sx q[0];
rz(0.87134937) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23587464) q[2];
sx q[2];
rz(-1.4067003) q[2];
sx q[2];
rz(0.99170384) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3672626) q[1];
sx q[1];
rz(-2.3891797) q[1];
sx q[1];
rz(-1.03536) q[1];
x q[2];
rz(-2.3865984) q[3];
sx q[3];
rz(-2.0557311) q[3];
sx q[3];
rz(-2.274853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(1.735431) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9091699) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(2.4160093) q[0];
rz(1.3369417) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(-0.57428378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030356) q[0];
sx q[0];
rz(-1.4113562) q[0];
sx q[0];
rz(3.0277059) q[0];
x q[1];
rz(2.1404187) q[2];
sx q[2];
rz(-1.9211848) q[2];
sx q[2];
rz(-2.1621665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5301104) q[1];
sx q[1];
rz(-2.2951295) q[1];
sx q[1];
rz(0.15165374) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3575555) q[3];
sx q[3];
rz(-1.136354) q[3];
sx q[3];
rz(1.8798352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(2.4794225) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.013423) q[1];
sx q[1];
rz(1.3319344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71026826) q[0];
sx q[0];
rz(-1.5402176) q[0];
sx q[0];
rz(1.7215183) q[0];
rz(-2.0831624) q[2];
sx q[2];
rz(-0.46736267) q[2];
sx q[2];
rz(0.99603727) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1292746) q[1];
sx q[1];
rz(-0.32806096) q[1];
sx q[1];
rz(-0.78719689) q[1];
rz(-pi) q[2];
rz(0.9040971) q[3];
sx q[3];
rz(-2.3257964) q[3];
sx q[3];
rz(-0.85532361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3146882) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(-0.17987128) q[2];
rz(-1.7449069) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708645) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(-2.6188359) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(1.1473473) q[2];
sx q[2];
rz(-0.97931391) q[2];
sx q[2];
rz(-0.58936832) q[2];
rz(1.9369851) q[3];
sx q[3];
rz(-2.0156751) q[3];
sx q[3];
rz(-1.0734476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
