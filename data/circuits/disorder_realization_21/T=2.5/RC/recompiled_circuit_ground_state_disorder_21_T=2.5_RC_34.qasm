OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(3.6977036) q[0];
sx q[0];
rz(7.2365427) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(-2.1929725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9069288) q[0];
sx q[0];
rz(-2.1000104) q[0];
sx q[0];
rz(-1.0906903) q[0];
rz(2.9041102) q[2];
sx q[2];
rz(-1.5735224) q[2];
sx q[2];
rz(-3.1196371) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1114914) q[1];
sx q[1];
rz(-1.9989357) q[1];
sx q[1];
rz(-0.74049048) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0417651) q[3];
sx q[3];
rz(-2.5210125) q[3];
sx q[3];
rz(2.4247501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3775776) q[2];
sx q[2];
rz(-0.074187584) q[2];
sx q[2];
rz(-0.78262502) q[2];
rz(-3.022656) q[3];
sx q[3];
rz(-1.0340034) q[3];
sx q[3];
rz(2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9453732) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(2.8837606) q[0];
rz(-0.091015426) q[1];
sx q[1];
rz(-2.0423404) q[1];
sx q[1];
rz(1.6450504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89423075) q[0];
sx q[0];
rz(-1.6391488) q[0];
sx q[0];
rz(-2.1240881) q[0];
rz(2.7690954) q[2];
sx q[2];
rz(-1.8939233) q[2];
sx q[2];
rz(-1.1748287) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5241201) q[1];
sx q[1];
rz(-1.6008428) q[1];
sx q[1];
rz(-1.9829795) q[1];
rz(-pi) q[2];
rz(-1.1608382) q[3];
sx q[3];
rz(-0.94469584) q[3];
sx q[3];
rz(-2.5253341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1220793) q[2];
sx q[2];
rz(-1.2747108) q[2];
sx q[2];
rz(-0.40461928) q[2];
rz(0.98958611) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(-2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.1125546) q[0];
sx q[0];
rz(-2.9974388) q[0];
sx q[0];
rz(-0.83830225) q[0];
rz(0.84838947) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(-0.99260509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557195) q[0];
sx q[0];
rz(-0.8536754) q[0];
sx q[0];
rz(1.3846057) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5136593) q[2];
sx q[2];
rz(-1.9046202) q[2];
sx q[2];
rz(-1.305507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4731208) q[1];
sx q[1];
rz(-1.8278012) q[1];
sx q[1];
rz(0.42196749) q[1];
rz(-2.8512673) q[3];
sx q[3];
rz(-1.464352) q[3];
sx q[3];
rz(2.8255879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2720211) q[2];
sx q[2];
rz(-1.973899) q[2];
sx q[2];
rz(2.0167548) q[2];
rz(-0.62075067) q[3];
sx q[3];
rz(-0.97095942) q[3];
sx q[3];
rz(0.11515215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426185) q[0];
sx q[0];
rz(-1.9816575) q[0];
sx q[0];
rz(-2.9677891) q[0];
rz(-2.4558892) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(2.3033843) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099412709) q[0];
sx q[0];
rz(-2.0297219) q[0];
sx q[0];
rz(2.6290429) q[0];
x q[1];
rz(-0.72553749) q[2];
sx q[2];
rz(-0.71070403) q[2];
sx q[2];
rz(2.2827471) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61889583) q[1];
sx q[1];
rz(-2.451183) q[1];
sx q[1];
rz(0.22367475) q[1];
rz(-0.075013667) q[3];
sx q[3];
rz(-2.677644) q[3];
sx q[3];
rz(1.2926276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.54245) q[2];
sx q[2];
rz(-1.4582381) q[2];
sx q[2];
rz(0.35169265) q[2];
rz(1.9463978) q[3];
sx q[3];
rz(-1.1870793) q[3];
sx q[3];
rz(-3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.675932) q[0];
sx q[0];
rz(-0.52023482) q[0];
sx q[0];
rz(-0.35292536) q[0];
rz(2.832761) q[1];
sx q[1];
rz(-1.0363204) q[1];
sx q[1];
rz(-3.1130863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089680044) q[0];
sx q[0];
rz(-1.3962922) q[0];
sx q[0];
rz(-2.1580218) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9602922) q[2];
sx q[2];
rz(-1.2151698) q[2];
sx q[2];
rz(0.86330044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18896477) q[1];
sx q[1];
rz(-0.82023925) q[1];
sx q[1];
rz(2.8431975) q[1];
rz(1.1358374) q[3];
sx q[3];
rz(-1.3490145) q[3];
sx q[3];
rz(1.3835761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9286524) q[2];
sx q[2];
rz(-0.063455909) q[2];
sx q[2];
rz(-2.1707936) q[2];
rz(-1.0468696) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(-0.48986062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(0.090959892) q[0];
rz(-1.92314) q[1];
sx q[1];
rz(-0.33088845) q[1];
sx q[1];
rz(2.5023696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1341742) q[0];
sx q[0];
rz(-1.8960516) q[0];
sx q[0];
rz(-1.7132617) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7090061) q[2];
sx q[2];
rz(-1.4497529) q[2];
sx q[2];
rz(-2.2908186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7057695) q[1];
sx q[1];
rz(-0.87638777) q[1];
sx q[1];
rz(2.4468876) q[1];
x q[2];
rz(-3.1396041) q[3];
sx q[3];
rz(-1.8973777) q[3];
sx q[3];
rz(-1.1359362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0794534) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(2.6677168) q[2];
rz(1.8114932) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(2.7661095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629958) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(-0.22258776) q[0];
rz(2.0780156) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(-1.9482013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7442856) q[0];
sx q[0];
rz(-2.3145876) q[0];
sx q[0];
rz(1.4124407) q[0];
rz(-pi) q[1];
rz(-1.0859126) q[2];
sx q[2];
rz(-0.25561324) q[2];
sx q[2];
rz(-0.26299082) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4970253) q[1];
sx q[1];
rz(-2.679832) q[1];
sx q[1];
rz(0.30739947) q[1];
rz(-pi) q[2];
rz(-2.1126595) q[3];
sx q[3];
rz(-0.30908424) q[3];
sx q[3];
rz(-2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1412763) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(-2.5013962) q[2];
rz(3.1196502) q[3];
sx q[3];
rz(-1.7317737) q[3];
sx q[3];
rz(0.35023165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.491965) q[0];
sx q[0];
rz(-1.3977298) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(-0.87977663) q[1];
sx q[1];
rz(-1.9089411) q[1];
sx q[1];
rz(-0.89281503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8917903) q[0];
sx q[0];
rz(-3.043104) q[0];
sx q[0];
rz(-0.47827025) q[0];
rz(-2.3273349) q[2];
sx q[2];
rz(-0.99254464) q[2];
sx q[2];
rz(-1.1922497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2485321) q[1];
sx q[1];
rz(-1.2888146) q[1];
sx q[1];
rz(2.3470013) q[1];
rz(-pi) q[2];
rz(2.625301) q[3];
sx q[3];
rz(-1.1974575) q[3];
sx q[3];
rz(0.0084127154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8799379) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(-3.0547764) q[2];
rz(0.9451198) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(1.1300348) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574808) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(0.064067319) q[0];
rz(-2.0059026) q[1];
sx q[1];
rz(-1.7013197) q[1];
sx q[1];
rz(-2.6293829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50822608) q[0];
sx q[0];
rz(-0.53589833) q[0];
sx q[0];
rz(0.6750489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2733803) q[2];
sx q[2];
rz(-1.3033452) q[2];
sx q[2];
rz(2.3600459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87742245) q[1];
sx q[1];
rz(-1.8296731) q[1];
sx q[1];
rz(1.6593462) q[1];
rz(-0.40887399) q[3];
sx q[3];
rz(-1.5561473) q[3];
sx q[3];
rz(1.15575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3247437) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(2.0762439) q[2];
rz(1.4893074) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(-0.092441946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866933) q[0];
sx q[0];
rz(-2.7296992) q[0];
sx q[0];
rz(0.33083415) q[0];
rz(1.9031485) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(-0.27473658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392259) q[0];
sx q[0];
rz(-2.5792092) q[0];
sx q[0];
rz(-2.2455162) q[0];
rz(-pi) q[1];
rz(1.6034699) q[2];
sx q[2];
rz(-1.297566) q[2];
sx q[2];
rz(0.88116026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.292784) q[1];
sx q[1];
rz(-1.0674137) q[1];
sx q[1];
rz(2.9817957) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.069271) q[3];
sx q[3];
rz(-2.157271) q[3];
sx q[3];
rz(-2.9972635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9773418) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(-1.3573307) q[2];
rz(-2.6409798) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(-1.4969426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5951344) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(3.0733227) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(2.7082607) q[2];
sx q[2];
rz(-2.7637252) q[2];
sx q[2];
rz(-2.6890474) q[2];
rz(2.3292062) q[3];
sx q[3];
rz(-2.8858622) q[3];
sx q[3];
rz(1.1271911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
