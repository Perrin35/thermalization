OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9165186) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(-2.201447) q[0];
rz(-pi) q[1];
rz(0.93809442) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(0.11703581) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7251016) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(-1.7822669) q[1];
rz(2.3834121) q[3];
sx q[3];
rz(-1.4361793) q[3];
sx q[3];
rz(-1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663548) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-2.2251341) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(0.8786456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40936138) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(-2.2344927) q[0];
rz(1.8354561) q[2];
sx q[2];
rz(-2.172643) q[2];
sx q[2];
rz(-2.2167609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2616927) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(-2.9628423) q[1];
rz(-pi) q[2];
rz(1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-2.0464954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(0.15163264) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(2.0203967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(1.0870766) q[0];
x q[1];
rz(-1.7240702) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(-0.085445554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.725863) q[1];
sx q[1];
rz(-1.2408857) q[1];
sx q[1];
rz(-1.8833453) q[1];
x q[2];
rz(-1.4533914) q[3];
sx q[3];
rz(-1.0798287) q[3];
sx q[3];
rz(-1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(0.78805584) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(0.11638164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7096036) q[0];
sx q[0];
rz(-2.5076712) q[0];
sx q[0];
rz(3.0984127) q[0];
rz(-pi) q[1];
rz(-2.0747244) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(2.4525814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0105646) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(-2.4478711) q[1];
rz(-pi) q[2];
rz(-2.651752) q[3];
sx q[3];
rz(-1.5252938) q[3];
sx q[3];
rz(-0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085885) q[0];
sx q[0];
rz(-2.5143393) q[0];
sx q[0];
rz(-1.6968615) q[0];
x q[1];
rz(-2.2448036) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(1.0587495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8685638) q[1];
sx q[1];
rz(-1.7595353) q[1];
sx q[1];
rz(-0.36004685) q[1];
x q[2];
rz(0.43907459) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(-2.4174158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5500568) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2287184) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(0.74761439) q[0];
rz(-pi) q[1];
rz(-2.4762819) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(0.88924185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1604662) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(1.6756945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19595512) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(-2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4326101) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(1.9463584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8853828) q[1];
sx q[1];
rz(-1.7269644) q[1];
sx q[1];
rz(3.0728818) q[1];
rz(-2.1738449) q[3];
sx q[3];
rz(-1.839404) q[3];
sx q[3];
rz(-2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.9205836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894665) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(0.37129398) q[0];
rz(-pi) q[1];
rz(-2.9646655) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(-2.4307876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3386821) q[1];
sx q[1];
rz(-1.7065587) q[1];
sx q[1];
rz(1.6822097) q[1];
x q[2];
rz(1.2113167) q[3];
sx q[3];
rz(-2.3593138) q[3];
sx q[3];
rz(-0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.5195297) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9916358) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.3791023) q[0];
x q[1];
rz(-0.33340402) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(2.8508027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77482624) q[1];
sx q[1];
rz(-2.3407196) q[1];
sx q[1];
rz(1.1714539) q[1];
rz(-0.59154193) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(-1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126497) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(1.2055231) q[0];
x q[1];
rz(2.0613725) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(2.0928004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0298437) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(0.44058056) q[1];
rz(1.2530008) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(-0.12249891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-2.649879) q[2];
sx q[2];
rz(-1.9946873) q[2];
sx q[2];
rz(-0.36110525) q[2];
rz(0.58581523) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
