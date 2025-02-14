OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(0.49459767) q[0];
rz(-5.1410723) q[1];
sx q[1];
rz(2.1358868) q[1];
sx q[1];
rz(11.436643) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0725192) q[0];
sx q[0];
rz(-1.7290218) q[0];
sx q[0];
rz(2.4618966) q[0];
rz(-pi) q[1];
rz(-2.6886875) q[2];
sx q[2];
rz(-2.0257086) q[2];
sx q[2];
rz(-2.186113) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5482481) q[1];
sx q[1];
rz(-1.3181837) q[1];
sx q[1];
rz(-2.3139364) q[1];
rz(-2.7143728) q[3];
sx q[3];
rz(-2.1764206) q[3];
sx q[3];
rz(0.96715121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8226681) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(0.85282105) q[2];
rz(1.3301814) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(-3.0947963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(-1.6774696) q[0];
rz(1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(-1.9333855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0532329) q[0];
sx q[0];
rz(-2.3908983) q[0];
sx q[0];
rz(2.4093767) q[0];
x q[1];
rz(2.0457532) q[2];
sx q[2];
rz(-1.37687) q[2];
sx q[2];
rz(-0.79481193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.959001) q[1];
sx q[1];
rz(-1.4495069) q[1];
sx q[1];
rz(-1.3022997) q[1];
rz(0.45470806) q[3];
sx q[3];
rz(-1.7268001) q[3];
sx q[3];
rz(-1.7059513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6055484) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.991796) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-0.27161828) q[0];
rz(-2.244921) q[1];
sx q[1];
rz(-2.6397557) q[1];
sx q[1];
rz(-1.2976049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51894278) q[0];
sx q[0];
rz(-2.033816) q[0];
sx q[0];
rz(0.011269022) q[0];
rz(-pi) q[1];
rz(2.6372725) q[2];
sx q[2];
rz(-1.9078622) q[2];
sx q[2];
rz(0.45318174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7097292) q[1];
sx q[1];
rz(-1.456902) q[1];
sx q[1];
rz(1.8191871) q[1];
x q[2];
rz(1.5758697) q[3];
sx q[3];
rz(-2.3327565) q[3];
sx q[3];
rz(0.047732959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4012332) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(-0.05376251) q[2];
rz(-1.2939804) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(-0.23538858) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2006705) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.4642375) q[0];
rz(1.1306521) q[1];
sx q[1];
rz(-0.73431763) q[1];
sx q[1];
rz(-0.11437036) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30442023) q[0];
sx q[0];
rz(-2.3995598) q[0];
sx q[0];
rz(-0.64226182) q[0];
x q[1];
rz(2.0691772) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(1.4630813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76941608) q[1];
sx q[1];
rz(-1.5972563) q[1];
sx q[1];
rz(-2.2464804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19802494) q[3];
sx q[3];
rz(-1.5772444) q[3];
sx q[3];
rz(-2.7355268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6681246) q[2];
sx q[2];
rz(-1.6056085) q[2];
sx q[2];
rz(-3.065897) q[2];
rz(2.7068052) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(-0.079631478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(2.721526) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(1.8992281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5899701) q[0];
sx q[0];
rz(-1.1484572) q[0];
sx q[0];
rz(0.7606272) q[0];
rz(-pi) q[1];
rz(-1.8625284) q[2];
sx q[2];
rz(-1.1559249) q[2];
sx q[2];
rz(1.3049098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1429726) q[1];
sx q[1];
rz(-3.0398453) q[1];
sx q[1];
rz(-1.7228026) q[1];
rz(-pi) q[2];
rz(0.22132921) q[3];
sx q[3];
rz(-1.5240655) q[3];
sx q[3];
rz(2.5525301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8215948) q[2];
sx q[2];
rz(-0.91511202) q[2];
sx q[2];
rz(2.7670009) q[2];
rz(1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11151611) q[0];
sx q[0];
rz(-0.498963) q[0];
sx q[0];
rz(2.7110355) q[0];
rz(0.14398362) q[1];
sx q[1];
rz(-2.0591683) q[1];
sx q[1];
rz(-0.63124257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3595561) q[0];
sx q[0];
rz(-0.61114531) q[0];
sx q[0];
rz(0.094159889) q[0];
rz(-1.8885884) q[2];
sx q[2];
rz(-0.43995198) q[2];
sx q[2];
rz(-0.41340128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5507372) q[1];
sx q[1];
rz(-0.8181347) q[1];
sx q[1];
rz(2.4356151) q[1];
x q[2];
rz(-2.1926375) q[3];
sx q[3];
rz(-2.6027711) q[3];
sx q[3];
rz(-2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1708019) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(-1.4405174) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.454527) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(0.92426306) q[0];
rz(2.848792) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(0.80088314) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0706884) q[0];
sx q[0];
rz(-1.3137806) q[0];
sx q[0];
rz(-1.1852253) q[0];
rz(-2.7960294) q[2];
sx q[2];
rz(-2.2032732) q[2];
sx q[2];
rz(1.3407624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2877601) q[1];
sx q[1];
rz(-1.4656855) q[1];
sx q[1];
rz(-1.6002602) q[1];
rz(-0.022014736) q[3];
sx q[3];
rz(-1.5436633) q[3];
sx q[3];
rz(-0.61329816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(-1.2804383) q[2];
rz(2.473623) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(0.41332301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.2385383) q[0];
sx q[0];
rz(1.1827693) q[0];
rz(-2.9044652) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(0.059344083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049669) q[0];
sx q[0];
rz(-1.4511746) q[0];
sx q[0];
rz(-1.2960394) q[0];
rz(-pi) q[1];
rz(-2.4164532) q[2];
sx q[2];
rz(-1.1584917) q[2];
sx q[2];
rz(2.7453842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0036111) q[1];
sx q[1];
rz(-1.6754158) q[1];
sx q[1];
rz(-1.305718) q[1];
rz(-pi) q[2];
rz(1.343325) q[3];
sx q[3];
rz(-2.8099647) q[3];
sx q[3];
rz(0.92078269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(1.7891368) q[2];
rz(1.3672359) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(-2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(-0.45968858) q[0];
rz(0.028060878) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(1.9557767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494649) q[0];
sx q[0];
rz(-0.26927265) q[0];
sx q[0];
rz(-2.9721391) q[0];
rz(2.2345951) q[2];
sx q[2];
rz(-0.12842783) q[2];
sx q[2];
rz(-2.3333486) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1255155) q[1];
sx q[1];
rz(-2.0488157) q[1];
sx q[1];
rz(-1.0325055) q[1];
x q[2];
rz(-2.5236058) q[3];
sx q[3];
rz(-1.3387965) q[3];
sx q[3];
rz(-2.8232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(2.9848177) q[2];
rz(2.5458941) q[3];
sx q[3];
rz(-2.7955293) q[3];
sx q[3];
rz(0.025207635) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6237685) q[0];
sx q[0];
rz(-1.0305923) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(-0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-2.877291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83885709) q[0];
sx q[0];
rz(-2.502901) q[0];
sx q[0];
rz(-2.6810718) q[0];
rz(-pi) q[1];
rz(-2.5223612) q[2];
sx q[2];
rz(-2.3546017) q[2];
sx q[2];
rz(-2.3004722) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42929998) q[1];
sx q[1];
rz(-1.4136864) q[1];
sx q[1];
rz(-0.18204851) q[1];
rz(-pi) q[2];
rz(-0.49510689) q[3];
sx q[3];
rz(-1.6551842) q[3];
sx q[3];
rz(-1.5051248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(-1.7372355) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(-0.082988113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(2.4304541) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(2.0696832) q[2];
sx q[2];
rz(-1.5317393) q[2];
sx q[2];
rz(-2.692937) q[2];
rz(2.5855999) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
