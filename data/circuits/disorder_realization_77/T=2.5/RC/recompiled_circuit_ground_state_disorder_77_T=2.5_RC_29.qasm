OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9935432) q[0];
sx q[0];
rz(-2.7033959) q[0];
sx q[0];
rz(2.4105657) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(0.55404034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2386411) q[0];
sx q[0];
rz(-1.1045505) q[0];
sx q[0];
rz(-2.58642) q[0];
rz(-pi) q[1];
rz(1.0921539) q[2];
sx q[2];
rz(-1.3502099) q[2];
sx q[2];
rz(-1.5719065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82137992) q[1];
sx q[1];
rz(-0.88372916) q[1];
sx q[1];
rz(-0.043619556) q[1];
x q[2];
rz(0.4783614) q[3];
sx q[3];
rz(-2.195092) q[3];
sx q[3];
rz(-2.9251298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76250166) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(0.16538922) q[2];
rz(-0.23073828) q[3];
sx q[3];
rz(-0.24300353) q[3];
sx q[3];
rz(-2.2802584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925258) q[0];
sx q[0];
rz(-2.8133744) q[0];
sx q[0];
rz(1.7829371) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(2.3014136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019842783) q[0];
sx q[0];
rz(-2.5917834) q[0];
sx q[0];
rz(2.0998831) q[0];
rz(-2.0303867) q[2];
sx q[2];
rz(-0.75558096) q[2];
sx q[2];
rz(-2.7019175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11586861) q[1];
sx q[1];
rz(-2.3626973) q[1];
sx q[1];
rz(2.7164357) q[1];
rz(-pi) q[2];
x q[2];
rz(1.718037) q[3];
sx q[3];
rz(-0.81063945) q[3];
sx q[3];
rz(-2.7642706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3885865) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(0.93442717) q[2];
rz(1.2855444) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(-0.88750315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11338209) q[0];
sx q[0];
rz(-1.0209571) q[0];
sx q[0];
rz(-1.5519979) q[0];
rz(1.3681083) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(1.1205463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5277953) q[0];
sx q[0];
rz(-0.58301914) q[0];
sx q[0];
rz(-0.9436508) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0650313) q[2];
sx q[2];
rz(-2.6496844) q[2];
sx q[2];
rz(3.0658403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9160794) q[1];
sx q[1];
rz(-0.54359964) q[1];
sx q[1];
rz(2.9617642) q[1];
rz(-pi) q[2];
rz(0.12579671) q[3];
sx q[3];
rz(-1.9748944) q[3];
sx q[3];
rz(0.27751291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0997194) q[2];
sx q[2];
rz(-2.9283044) q[2];
sx q[2];
rz(2.8495157) q[2];
rz(-1.9000351) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860155) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(0.43169942) q[0];
rz(0.15402928) q[1];
sx q[1];
rz(-1.5715716) q[1];
sx q[1];
rz(2.0808751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6650603) q[0];
sx q[0];
rz(-1.8899584) q[0];
sx q[0];
rz(1.9299797) q[0];
x q[1];
rz(0.94580146) q[2];
sx q[2];
rz(-2.8428787) q[2];
sx q[2];
rz(0.83160366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0077483245) q[1];
sx q[1];
rz(-0.40002003) q[1];
sx q[1];
rz(-2.176214) q[1];
rz(-0.68935518) q[3];
sx q[3];
rz(-2.8971147) q[3];
sx q[3];
rz(-1.0058798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3920307) q[2];
sx q[2];
rz(-1.6101086) q[2];
sx q[2];
rz(-2.5677666) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.6618238) q[3];
sx q[3];
rz(2.3315232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5350128) q[0];
sx q[0];
rz(-0.25595328) q[0];
sx q[0];
rz(-2.8888597) q[0];
rz(2.9622954) q[1];
sx q[1];
rz(-1.3456234) q[1];
sx q[1];
rz(-1.4089233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104467) q[0];
sx q[0];
rz(-0.24872336) q[0];
sx q[0];
rz(1.876271) q[0];
rz(-pi) q[1];
rz(-0.42667146) q[2];
sx q[2];
rz(-1.2152078) q[2];
sx q[2];
rz(-1.5384962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21633575) q[1];
sx q[1];
rz(-0.79265037) q[1];
sx q[1];
rz(1.3986194) q[1];
rz(-pi) q[2];
rz(-1.9305655) q[3];
sx q[3];
rz(-1.0749987) q[3];
sx q[3];
rz(-1.825765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0130284) q[2];
sx q[2];
rz(-1.7594254) q[2];
sx q[2];
rz(1.2010835) q[2];
rz(-0.80738336) q[3];
sx q[3];
rz(-1.7816593) q[3];
sx q[3];
rz(1.132563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0627237) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(0.8412745) q[0];
rz(-1.9367283) q[1];
sx q[1];
rz(-1.8236022) q[1];
sx q[1];
rz(2.1119609) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1618274) q[0];
sx q[0];
rz(-1.0603535) q[0];
sx q[0];
rz(-3.1296041) q[0];
rz(-pi) q[1];
rz(0.44544051) q[2];
sx q[2];
rz(-2.5800319) q[2];
sx q[2];
rz(-0.5024006) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8334044) q[1];
sx q[1];
rz(-0.17886111) q[1];
sx q[1];
rz(0.043828242) q[1];
rz(-pi) q[2];
rz(-2.8383202) q[3];
sx q[3];
rz(-1.7494043) q[3];
sx q[3];
rz(-2.957119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.131375) q[2];
sx q[2];
rz(-0.75675941) q[2];
sx q[2];
rz(-2.6436464) q[2];
rz(-2.61854) q[3];
sx q[3];
rz(-1.4191351) q[3];
sx q[3];
rz(0.12652346) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.842857) q[0];
sx q[0];
rz(-0.13164483) q[0];
sx q[0];
rz(-0.81197062) q[0];
rz(-0.89093351) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-0.99937159) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3131789) q[0];
sx q[0];
rz(-1.1346319) q[0];
sx q[0];
rz(0.51491134) q[0];
rz(-pi) q[1];
rz(1.4677202) q[2];
sx q[2];
rz(-1.4749881) q[2];
sx q[2];
rz(-0.16703781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33617299) q[1];
sx q[1];
rz(-0.78485705) q[1];
sx q[1];
rz(-2.428167) q[1];
x q[2];
rz(2.8511397) q[3];
sx q[3];
rz(-0.92602611) q[3];
sx q[3];
rz(-1.6381263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.823395) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(1.8434175) q[2];
rz(3.1031109) q[3];
sx q[3];
rz(-1.1063856) q[3];
sx q[3];
rz(1.6160256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0865974) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(0.31563345) q[0];
rz(-1.2339969) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(1.2589781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4988853) q[0];
sx q[0];
rz(-0.85432893) q[0];
sx q[0];
rz(-1.1906719) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.744049) q[2];
sx q[2];
rz(-0.8129186) q[2];
sx q[2];
rz(-0.11101857) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8371305) q[1];
sx q[1];
rz(-1.9726957) q[1];
sx q[1];
rz(2.7086651) q[1];
rz(-0.057825967) q[3];
sx q[3];
rz(-0.10198051) q[3];
sx q[3];
rz(-2.0174055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9180318) q[2];
sx q[2];
rz(-0.74350244) q[2];
sx q[2];
rz(-3.1136759) q[2];
rz(-0.74808407) q[3];
sx q[3];
rz(-1.4192162) q[3];
sx q[3];
rz(2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41392031) q[0];
sx q[0];
rz(-2.1509009) q[0];
sx q[0];
rz(-2.682611) q[0];
rz(1.1363632) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(-2.9232025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439664) q[0];
sx q[0];
rz(-1.5541926) q[0];
sx q[0];
rz(0.023331716) q[0];
rz(-pi) q[1];
x q[1];
rz(1.893546) q[2];
sx q[2];
rz(-2.7392929) q[2];
sx q[2];
rz(2.2108271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11248828) q[1];
sx q[1];
rz(-2.6262754) q[1];
sx q[1];
rz(0.30179939) q[1];
x q[2];
rz(1.8828527) q[3];
sx q[3];
rz(-2.7358905) q[3];
sx q[3];
rz(-0.73075529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8454664) q[2];
sx q[2];
rz(-0.25456905) q[2];
sx q[2];
rz(2.0938342) q[2];
rz(0.74538499) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(-1.4926225) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02136852) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-2.7987203) q[0];
rz(-2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(2.6284133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1091217) q[0];
sx q[0];
rz(-2.4801804) q[0];
sx q[0];
rz(1.830807) q[0];
rz(-pi) q[1];
rz(-2.1634899) q[2];
sx q[2];
rz(-0.58916559) q[2];
sx q[2];
rz(1.1617253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3575705) q[1];
sx q[1];
rz(-2.5139132) q[1];
sx q[1];
rz(0.72369544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9916198) q[3];
sx q[3];
rz(-1.4412035) q[3];
sx q[3];
rz(-2.0146927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7036983) q[2];
sx q[2];
rz(-2.5056705) q[2];
sx q[2];
rz(-0.61100125) q[2];
rz(0.25887394) q[3];
sx q[3];
rz(-1.2114108) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543058) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(2.1517131) q[1];
sx q[1];
rz(-1.9342593) q[1];
sx q[1];
rz(0.63607279) q[1];
rz(-1.5237332) q[2];
sx q[2];
rz(-0.18885352) q[2];
sx q[2];
rz(0.39299008) q[2];
rz(-0.69901917) q[3];
sx q[3];
rz(-1.5200281) q[3];
sx q[3];
rz(2.716223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
