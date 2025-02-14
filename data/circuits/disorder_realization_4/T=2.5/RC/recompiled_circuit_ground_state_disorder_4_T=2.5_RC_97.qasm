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
rz(3.562869) q[1];
sx q[1];
rz(8.1864551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8684611) q[0];
sx q[0];
rz(-2.1885311) q[0];
sx q[0];
rz(0.6427838) q[0];
x q[1];
rz(-0.86482817) q[2];
sx q[2];
rz(-2.6061686) q[2];
sx q[2];
rz(0.12685093) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4290966) q[1];
sx q[1];
rz(-1.1084021) q[1];
sx q[1];
rz(0.40689792) q[1];
rz(-1.7389033) q[3];
sx q[3];
rz(-1.013275) q[3];
sx q[3];
rz(-1.7926737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51332659) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(2.4386621) q[2];
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
x q[2];
rz(-pi/2) q[3];
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
rz(-1.4502451) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(2.3968089) q[0];
rz(1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(0.94211284) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4018135) q[0];
sx q[0];
rz(-1.3688068) q[0];
sx q[0];
rz(2.3752579) q[0];
rz(-pi) q[1];
rz(-3.0622185) q[2];
sx q[2];
rz(-2.3217776) q[2];
sx q[2];
rz(0.83460966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5444138) q[1];
sx q[1];
rz(-2.326818) q[1];
sx q[1];
rz(-0.18715231) q[1];
x q[2];
rz(-3.0511176) q[3];
sx q[3];
rz(-0.58808413) q[3];
sx q[3];
rz(-0.17234853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(-1.6628294) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11109322) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(-1.4053364) q[0];
rz(-1.7449215) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(-0.90000802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4858682) q[0];
sx q[0];
rz(-1.9330171) q[0];
sx q[0];
rz(-0.78804228) q[0];
x q[1];
rz(0.75863691) q[2];
sx q[2];
rz(-2.4066145) q[2];
sx q[2];
rz(2.1211565) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72199539) q[1];
sx q[1];
rz(-1.743814) q[1];
sx q[1];
rz(-0.74651511) q[1];
rz(0.7358968) q[3];
sx q[3];
rz(-0.80965878) q[3];
sx q[3];
rz(2.5084605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46093837) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(0.94949618) q[2];
rz(-2.5203868) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7200274) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(3.0928639) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(-0.31160942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8870114) q[0];
sx q[0];
rz(-1.100228) q[0];
sx q[0];
rz(-0.30593095) q[0];
rz(0.34398244) q[2];
sx q[2];
rz(-2.0728353) q[2];
sx q[2];
rz(0.12884049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5587363) q[1];
sx q[1];
rz(-1.3060102) q[1];
sx q[1];
rz(-2.2469421) q[1];
rz(-pi) q[2];
rz(-0.29693691) q[3];
sx q[3];
rz(-1.9949759) q[3];
sx q[3];
rz(-2.3564828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17778808) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(1.6507899) q[2];
rz(-2.7866411) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(-1.8420334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0519003) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(-1.2667013) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(-1.5142534) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86858803) q[0];
sx q[0];
rz(-2.3914861) q[0];
sx q[0];
rz(2.0387893) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1852988) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(-2.1182107) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5473289) q[1];
sx q[1];
rz(-1.0272861) q[1];
sx q[1];
rz(-0.059652358) q[1];
rz(-pi) q[2];
rz(1.5290401) q[3];
sx q[3];
rz(-0.72894086) q[3];
sx q[3];
rz(1.6870013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2221471) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670052) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(2.7401155) q[0];
rz(-2.0145156) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(-2.3289767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81174034) q[0];
sx q[0];
rz(-0.63416687) q[0];
sx q[0];
rz(1.5777753) q[0];
rz(-pi) q[1];
rz(-0.70131371) q[2];
sx q[2];
rz(-0.35365401) q[2];
sx q[2];
rz(2.3344699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5548812) q[1];
sx q[1];
rz(-0.63535345) q[1];
sx q[1];
rz(2.79252) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3964368) q[3];
sx q[3];
rz(-2.0260915) q[3];
sx q[3];
rz(2.8139092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(2.9511792) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.93369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(2.3741123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90119367) q[0];
sx q[0];
rz(-1.6220105) q[0];
sx q[0];
rz(0.9238433) q[0];
rz(-pi) q[1];
rz(-1.4857471) q[2];
sx q[2];
rz(-1.1901996) q[2];
sx q[2];
rz(-2.2589661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2042522) q[1];
sx q[1];
rz(-0.91751999) q[1];
sx q[1];
rz(3.0863347) q[1];
rz(-pi) q[2];
rz(-1.5753059) q[3];
sx q[3];
rz(-1.1645551) q[3];
sx q[3];
rz(-2.3227228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(1.5896612) q[2];
rz(-1.6520366) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(-1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234696) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(2.7662011) q[0];
rz(-2.2581532) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(2.4129131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7405072) q[0];
sx q[0];
rz(-0.87941636) q[0];
sx q[0];
rz(-2.2702433) q[0];
rz(-1.402114) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(-2.6017435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9364215) q[1];
sx q[1];
rz(-1.2146307) q[1];
sx q[1];
rz(0.89295279) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75499423) q[3];
sx q[3];
rz(-2.0557311) q[3];
sx q[3];
rz(-0.86673966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(-1.735431) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(-1.5555443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9091699) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(2.4160093) q[0];
rz(-1.8046509) q[1];
sx q[1];
rz(-1.888211) q[1];
sx q[1];
rz(-2.5673089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5142412) q[0];
sx q[0];
rz(-2.9459369) q[0];
sx q[0];
rz(2.186004) q[0];
rz(0.4094643) q[2];
sx q[2];
rz(-2.1019693) q[2];
sx q[2];
rz(-2.766618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7567819) q[1];
sx q[1];
rz(-0.73721209) q[1];
sx q[1];
rz(-1.739915) q[1];
x q[2];
rz(0.78403715) q[3];
sx q[3];
rz(-1.136354) q[3];
sx q[3];
rz(-1.2617574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1672704) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(-0.66217011) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(-1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(-1.3978488) q[0];
rz(-1.0700048) q[1];
sx q[1];
rz(-2.013423) q[1];
sx q[1];
rz(1.3319344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66187132) q[0];
sx q[0];
rz(-2.9878231) q[0];
sx q[0];
rz(1.3698306) q[0];
rz(-pi) q[1];
rz(0.24253129) q[2];
sx q[2];
rz(-1.9743391) q[2];
sx q[2];
rz(-1.5582486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20132682) q[1];
sx q[1];
rz(-1.3405217) q[1];
sx q[1];
rz(0.23576945) q[1];
x q[2];
rz(-2.5601848) q[3];
sx q[3];
rz(-2.1801342) q[3];
sx q[3];
rz(1.4319624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(-0.17987128) q[2];
rz(-1.7449069) q[3];
sx q[3];
rz(-1.4254009) q[3];
sx q[3];
rz(0.21656187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.77072813) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(0.54900563) q[2];
sx q[2];
rz(-0.71239757) q[2];
sx q[2];
rz(-3.0511643) q[2];
rz(-0.47209817) q[3];
sx q[3];
rz(-1.8999128) q[3];
sx q[3];
rz(0.66090665) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
