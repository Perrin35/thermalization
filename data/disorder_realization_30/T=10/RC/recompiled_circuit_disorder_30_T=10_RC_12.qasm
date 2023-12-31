OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
x q[1];
rz(2.1985487) q[2];
sx q[2];
rz(-0.53484166) q[2];
sx q[2];
rz(-0.56378555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1026417) q[1];
sx q[1];
rz(-2.0248374) q[1];
sx q[1];
rz(2.7137043) q[1];
rz(-pi) q[2];
rz(-0.37813152) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(-3.0331628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(2.7976024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.7864236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8665168) q[0];
sx q[0];
rz(-1.5059024) q[0];
sx q[0];
rz(2.1823723) q[0];
x q[1];
rz(0.16263527) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(1.4974809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.018108798) q[1];
sx q[1];
rz(-1.6347486) q[1];
sx q[1];
rz(-1.121184) q[1];
rz(-pi) q[2];
x q[2];
rz(1.387272) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(-0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.5887567) q[0];
x q[1];
rz(0.89945729) q[2];
sx q[2];
rz(-2.3231635) q[2];
sx q[2];
rz(2.1607272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.47485891) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(2.1875847) q[1];
rz(-2.2933526) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7751559) q[0];
sx q[0];
rz(-1.1661134) q[0];
sx q[0];
rz(-1.0767656) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3035994) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(-2.879564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30918446) q[1];
sx q[1];
rz(-0.60317457) q[1];
sx q[1];
rz(0.96997728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7632145) q[3];
sx q[3];
rz(-1.2580066) q[3];
sx q[3];
rz(-2.6421412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0793314) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(0.22236951) q[0];
rz(-0.58381501) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(2.3655287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27969589) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(-2.8054603) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1760686) q[3];
sx q[3];
rz(-2*pi/13) q[3];
sx q[3];
rz(2.0616639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.8959321) q[2];
rz(-1.2549531) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(0.83827034) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(-2.9340414) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-2.025827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0633495) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(-2.1372165) q[0];
rz(2.0954779) q[2];
sx q[2];
rz(-2.4002053) q[2];
sx q[2];
rz(1.807883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9876154) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(1.5495367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7498155) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(-0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(-2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-0.66837464) q[0];
sx q[0];
rz(-1.5947123) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-3.0722741) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.454969) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(-1.7186233) q[1];
x q[2];
rz(0.087248487) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(-1.6202392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103219) q[0];
sx q[0];
rz(-1.6390419) q[0];
sx q[0];
rz(-0.1971498) q[0];
rz(-pi) q[1];
rz(0.40201681) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(1.6059665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0115067) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(1.4873051) q[1];
rz(-1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(0.27063453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5198869) q[0];
sx q[0];
rz(-1.5396376) q[0];
sx q[0];
rz(3.1057538) q[0];
x q[1];
rz(2.6793849) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(0.30009899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2808025) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(-1.6153107) q[1];
rz(-pi) q[2];
rz(-2.4828033) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(-2.7487019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5737168) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(-0.042851187) q[0];
x q[1];
rz(0.020936326) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(-2.7431969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(0.93025031) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(-1.0333992) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-2.9359948) q[2];
sx q[2];
rz(-1.0369221) q[2];
sx q[2];
rz(2.592929) q[2];
rz(-2.3776688) q[3];
sx q[3];
rz(-2.7803956) q[3];
sx q[3];
rz(-1.9261388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
