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
rz(-1.3136366) q[0];
sx q[0];
rz(-3.0744636) q[0];
sx q[0];
rz(3.1188174) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(0.013962362) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3168105) q[0];
sx q[0];
rz(-0.75763065) q[0];
sx q[0];
rz(3.0824721) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3910782) q[2];
sx q[2];
rz(-0.57518164) q[2];
sx q[2];
rz(1.2746948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96367946) q[1];
sx q[1];
rz(-1.3909893) q[1];
sx q[1];
rz(-0.20188971) q[1];
rz(-pi) q[2];
rz(0.45043378) q[3];
sx q[3];
rz(-1.2518468) q[3];
sx q[3];
rz(-1.3802841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17249168) q[2];
sx q[2];
rz(-1.3250985) q[2];
sx q[2];
rz(1.3410404) q[2];
rz(-1.4260882) q[3];
sx q[3];
rz(-2.0243578) q[3];
sx q[3];
rz(-2.8360227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7476927) q[0];
sx q[0];
rz(-1.1423528) q[0];
sx q[0];
rz(2.4161412) q[0];
rz(-2.2254288) q[1];
sx q[1];
rz(-2.3326645) q[1];
sx q[1];
rz(2.5221672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9206132) q[0];
sx q[0];
rz(-1.0695463) q[0];
sx q[0];
rz(-0.46904012) q[0];
rz(-pi) q[1];
rz(-2.1690363) q[2];
sx q[2];
rz(-2.2248039) q[2];
sx q[2];
rz(-1.1197661) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5350772) q[1];
sx q[1];
rz(-1.9961352) q[1];
sx q[1];
rz(-0.40014793) q[1];
rz(-2.0158012) q[3];
sx q[3];
rz(-1.3135305) q[3];
sx q[3];
rz(0.98350888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3118423) q[2];
sx q[2];
rz(-2.1205015) q[2];
sx q[2];
rz(-2.1799083) q[2];
rz(2.035615) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(2.7963514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2894534) q[0];
sx q[0];
rz(-2.0134605) q[0];
sx q[0];
rz(-1.7732675) q[0];
rz(0.013956919) q[1];
sx q[1];
rz(-2.0266666) q[1];
sx q[1];
rz(1.4143292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7507259) q[0];
sx q[0];
rz(-1.8818151) q[0];
sx q[0];
rz(3.0469131) q[0];
rz(-1.2017085) q[2];
sx q[2];
rz(-1.6371583) q[2];
sx q[2];
rz(-1.8770799) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77905475) q[1];
sx q[1];
rz(-1.0666872) q[1];
sx q[1];
rz(0.49493954) q[1];
rz(-pi) q[2];
rz(2.7501864) q[3];
sx q[3];
rz(-1.3308755) q[3];
sx q[3];
rz(-0.6894043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76325512) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(0.21667996) q[2];
rz(-2.2583151) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(2.0047552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20275177) q[0];
sx q[0];
rz(-2.1348248) q[0];
sx q[0];
rz(0.83813611) q[0];
rz(-2.5996767) q[1];
sx q[1];
rz(-2.7653265) q[1];
sx q[1];
rz(0.045104973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7411557) q[0];
sx q[0];
rz(-0.61624762) q[0];
sx q[0];
rz(-0.23016302) q[0];
rz(-pi) q[1];
rz(0.25928478) q[2];
sx q[2];
rz(-2.0392937) q[2];
sx q[2];
rz(1.1069654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20238961) q[1];
sx q[1];
rz(-1.09859) q[1];
sx q[1];
rz(3.0665728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4724422) q[3];
sx q[3];
rz(-1.0788267) q[3];
sx q[3];
rz(1.6823925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1230459) q[2];
sx q[2];
rz(-0.16877731) q[2];
sx q[2];
rz(0.65400845) q[2];
rz(0.6319913) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(-0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.80394799) q[0];
sx q[0];
rz(-0.73886442) q[0];
sx q[0];
rz(1.7919354) q[0];
rz(0.82756394) q[1];
sx q[1];
rz(-2.5065828) q[1];
sx q[1];
rz(2.6571224) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0837528) q[0];
sx q[0];
rz(-2.6415402) q[0];
sx q[0];
rz(-2.3127348) q[0];
x q[1];
rz(1.7155649) q[2];
sx q[2];
rz(-2.5826404) q[2];
sx q[2];
rz(1.1048855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4836241) q[1];
sx q[1];
rz(-1.4438902) q[1];
sx q[1];
rz(-2.7183141) q[1];
rz(-pi) q[2];
rz(2.2665759) q[3];
sx q[3];
rz(-1.9364803) q[3];
sx q[3];
rz(-2.65923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95939246) q[2];
sx q[2];
rz(-1.0687989) q[2];
sx q[2];
rz(-0.63641316) q[2];
rz(3.0427129) q[3];
sx q[3];
rz(-1.585588) q[3];
sx q[3];
rz(1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8168378) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.9689993) q[0];
rz(1.532754) q[1];
sx q[1];
rz(-2.1295261) q[1];
sx q[1];
rz(0.00072678725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4945723) q[0];
sx q[0];
rz(-1.9969517) q[0];
sx q[0];
rz(2.1151353) q[0];
x q[1];
rz(-0.41623314) q[2];
sx q[2];
rz(-2.9011167) q[2];
sx q[2];
rz(-1.3593624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24405133) q[1];
sx q[1];
rz(-1.6836299) q[1];
sx q[1];
rz(-0.33399069) q[1];
rz(-pi) q[2];
rz(-0.29720593) q[3];
sx q[3];
rz(-0.5999705) q[3];
sx q[3];
rz(1.6855437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8751004) q[2];
sx q[2];
rz(-2.2537587) q[2];
sx q[2];
rz(0.20981851) q[2];
rz(0.80875129) q[3];
sx q[3];
rz(-2.5918312) q[3];
sx q[3];
rz(-0.89890283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631735) q[0];
sx q[0];
rz(-1.2043948) q[0];
sx q[0];
rz(-3.0294982) q[0];
rz(-3.0922999) q[1];
sx q[1];
rz(-0.53105989) q[1];
sx q[1];
rz(-1.8045527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508731) q[0];
sx q[0];
rz(-2.4639122) q[0];
sx q[0];
rz(-1.2356133) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12814204) q[2];
sx q[2];
rz(-1.0661835) q[2];
sx q[2];
rz(-2.5487633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9328135) q[1];
sx q[1];
rz(-1.5808388) q[1];
sx q[1];
rz(2.7947763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2483761) q[3];
sx q[3];
rz(-2.5395576) q[3];
sx q[3];
rz(2.4436434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2992531) q[2];
sx q[2];
rz(-2.3167593) q[2];
sx q[2];
rz(2.4315368) q[2];
rz(3.1402785) q[3];
sx q[3];
rz(-0.90130663) q[3];
sx q[3];
rz(2.5406751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32427078) q[0];
sx q[0];
rz(-1.3794427) q[0];
sx q[0];
rz(-0.58260179) q[0];
rz(2.461589) q[1];
sx q[1];
rz(-2.1425207) q[1];
sx q[1];
rz(-1.2635788) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1156968) q[0];
sx q[0];
rz(-1.8531593) q[0];
sx q[0];
rz(2.6656194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9367909) q[2];
sx q[2];
rz(-1.3582841) q[2];
sx q[2];
rz(0.16016618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7192051) q[1];
sx q[1];
rz(-2.4610956) q[1];
sx q[1];
rz(2.8533555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6035147) q[3];
sx q[3];
rz(-1.9944291) q[3];
sx q[3];
rz(-0.7684497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6363643) q[2];
sx q[2];
rz(-1.4437081) q[2];
sx q[2];
rz(-0.12602028) q[2];
rz(-1.5382918) q[3];
sx q[3];
rz(-2.3647629) q[3];
sx q[3];
rz(0.34415027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6525604) q[0];
sx q[0];
rz(-1.5131938) q[0];
sx q[0];
rz(-1.9788347) q[0];
rz(-1.3105505) q[1];
sx q[1];
rz(-0.73859155) q[1];
sx q[1];
rz(-1.7576677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726295) q[0];
sx q[0];
rz(-1.2314737) q[0];
sx q[0];
rz(-3.0940381) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6114386) q[2];
sx q[2];
rz(-2.9994476) q[2];
sx q[2];
rz(1.7696385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3466665) q[1];
sx q[1];
rz(-0.62267471) q[1];
sx q[1];
rz(-1.7706448) q[1];
rz(-3.0988068) q[3];
sx q[3];
rz(-1.7093911) q[3];
sx q[3];
rz(-0.78857546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80013529) q[2];
sx q[2];
rz(-1.2722445) q[2];
sx q[2];
rz(0.47194353) q[2];
rz(2.3927472) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(-0.81764618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6227459) q[0];
sx q[0];
rz(-0.8569583) q[0];
sx q[0];
rz(-2.6527606) q[0];
rz(-0.34293109) q[1];
sx q[1];
rz(-2.2551426) q[1];
sx q[1];
rz(2.0874646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3913522) q[0];
sx q[0];
rz(-3.0778193) q[0];
sx q[0];
rz(-0.94185726) q[0];
rz(-1.2230049) q[2];
sx q[2];
rz(-0.985983) q[2];
sx q[2];
rz(-2.0690167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78275241) q[1];
sx q[1];
rz(-0.61780518) q[1];
sx q[1];
rz(2.0899523) q[1];
rz(-pi) q[2];
rz(-3.0929186) q[3];
sx q[3];
rz(-1.0992388) q[3];
sx q[3];
rz(0.31192985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1663345) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(-3.0496822) q[2];
rz(-0.46936834) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(3.0288467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54409201) q[0];
sx q[0];
rz(-1.3213897) q[0];
sx q[0];
rz(-1.8505009) q[0];
rz(0.76580936) q[1];
sx q[1];
rz(-1.3594834) q[1];
sx q[1];
rz(1.8654738) q[1];
rz(-2.886404) q[2];
sx q[2];
rz(-1.219426) q[2];
sx q[2];
rz(0.30795369) q[2];
rz(2.6041531) q[3];
sx q[3];
rz(-1.7949355) q[3];
sx q[3];
rz(-1.8209281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
