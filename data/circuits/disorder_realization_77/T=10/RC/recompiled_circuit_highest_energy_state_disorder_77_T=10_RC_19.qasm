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
rz(-0.51802975) q[0];
sx q[0];
rz(-2.0002444) q[0];
sx q[0];
rz(-0.040123392) q[0];
rz(-2.8934381) q[1];
sx q[1];
rz(-2.0791972) q[1];
sx q[1];
rz(0.084029347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798767) q[0];
sx q[0];
rz(-1.6296224) q[0];
sx q[0];
rz(-0.071278871) q[0];
rz(0.40496396) q[2];
sx q[2];
rz(-2.7128007) q[2];
sx q[2];
rz(2.2909209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99043118) q[1];
sx q[1];
rz(-2.6365981) q[1];
sx q[1];
rz(-2.7846365) q[1];
rz(-2.5474882) q[3];
sx q[3];
rz(-0.66451529) q[3];
sx q[3];
rz(-2.6474109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4485126) q[2];
sx q[2];
rz(-1.3658407) q[2];
sx q[2];
rz(2.8903294) q[2];
rz(1.9692028) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(-3.0917061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6867111) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(1.3739817) q[0];
rz(-0.33085597) q[1];
sx q[1];
rz(-1.4127981) q[1];
sx q[1];
rz(-0.27423283) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170256) q[0];
sx q[0];
rz(-2.8170878) q[0];
sx q[0];
rz(-2.9175678) q[0];
rz(3.129446) q[2];
sx q[2];
rz(-0.53330219) q[2];
sx q[2];
rz(-2.8664194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8195023) q[1];
sx q[1];
rz(-1.7871457) q[1];
sx q[1];
rz(-2.809193) q[1];
rz(-pi) q[2];
rz(1.7555439) q[3];
sx q[3];
rz(-2.7165301) q[3];
sx q[3];
rz(2.5109072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26068822) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(1.0880967) q[2];
rz(1.043383) q[3];
sx q[3];
rz(-0.16845307) q[3];
sx q[3];
rz(0.25922957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.008721) q[0];
sx q[0];
rz(-2.6368124) q[0];
sx q[0];
rz(-1.1016499) q[0];
rz(2.3654826) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(2.4993842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5284536) q[0];
sx q[0];
rz(-2.0106843) q[0];
sx q[0];
rz(0.47275193) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8181591) q[2];
sx q[2];
rz(-0.86630922) q[2];
sx q[2];
rz(3.041628) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8485198) q[1];
sx q[1];
rz(-0.43918931) q[1];
sx q[1];
rz(-1.3249012) q[1];
x q[2];
rz(2.4829743) q[3];
sx q[3];
rz(-2.082654) q[3];
sx q[3];
rz(-0.030125253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3099826) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(2.2115121) q[2];
rz(-0.81405226) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(-1.9109776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.880068) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(0.50450182) q[0];
rz(-0.91744676) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(-0.20828542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513826) q[0];
sx q[0];
rz(-1.8135895) q[0];
sx q[0];
rz(1.7240748) q[0];
x q[1];
rz(-1.850205) q[2];
sx q[2];
rz(-0.87051755) q[2];
sx q[2];
rz(-0.4513739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.444297) q[1];
sx q[1];
rz(-1.0047067) q[1];
sx q[1];
rz(2.5725911) q[1];
rz(-pi) q[2];
rz(-0.39842968) q[3];
sx q[3];
rz(-0.59729415) q[3];
sx q[3];
rz(0.78890991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.31671277) q[2];
sx q[2];
rz(-1.0944288) q[2];
sx q[2];
rz(-1.1701976) q[2];
rz(2.1757388) q[3];
sx q[3];
rz(-2.070277) q[3];
sx q[3];
rz(-2.9962208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.55249864) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(2.476995) q[0];
rz(2.873114) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(0.99884117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132709) q[0];
sx q[0];
rz(-1.7230464) q[0];
sx q[0];
rz(-2.1571742) q[0];
rz(1.8605609) q[2];
sx q[2];
rz(-0.43998566) q[2];
sx q[2];
rz(0.54674613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70390648) q[1];
sx q[1];
rz(-1.4515299) q[1];
sx q[1];
rz(1.4005303) q[1];
rz(-pi) q[2];
rz(-2.2625974) q[3];
sx q[3];
rz(-1.5546518) q[3];
sx q[3];
rz(-0.059034928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(3.0244136) q[2];
rz(-3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(-0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(-1.2600979) q[0];
rz(3.0267808) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(3.0462435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95923808) q[0];
sx q[0];
rz(-1.7732014) q[0];
sx q[0];
rz(1.5558045) q[0];
x q[1];
rz(-1.8666519) q[2];
sx q[2];
rz(-2.618578) q[2];
sx q[2];
rz(-0.77525866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8417175) q[1];
sx q[1];
rz(-0.63427329) q[1];
sx q[1];
rz(-2.3298765) q[1];
x q[2];
rz(1.5891777) q[3];
sx q[3];
rz(-0.52747969) q[3];
sx q[3];
rz(-0.35273409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0872588) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(0.66162649) q[2];
rz(2.3719487) q[3];
sx q[3];
rz(-1.6911643) q[3];
sx q[3];
rz(0.86696398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374283) q[0];
sx q[0];
rz(-2.7935109) q[0];
sx q[0];
rz(2.3837756) q[0];
rz(0.025253145) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(-1.1753561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379743) q[0];
sx q[0];
rz(-2.0424065) q[0];
sx q[0];
rz(1.1704117) q[0];
x q[1];
rz(-2.6525431) q[2];
sx q[2];
rz(-2.1903775) q[2];
sx q[2];
rz(1.3330158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84056329) q[1];
sx q[1];
rz(-1.1811273) q[1];
sx q[1];
rz(2.3430166) q[1];
rz(-pi) q[2];
rz(-2.712599) q[3];
sx q[3];
rz(-1.5499664) q[3];
sx q[3];
rz(0.57383895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56388277) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(1.2152524) q[2];
rz(-0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(1.5569713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4828846) q[0];
sx q[0];
rz(-1.2889129) q[0];
sx q[0];
rz(-0.38920745) q[0];
rz(0.088118531) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0813839) q[0];
sx q[0];
rz(-1.7870149) q[0];
sx q[0];
rz(0.24198089) q[0];
x q[1];
rz(-1.1168209) q[2];
sx q[2];
rz(-2.1810227) q[2];
sx q[2];
rz(2.7895841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2417422) q[1];
sx q[1];
rz(-2.012737) q[1];
sx q[1];
rz(-1.5133218) q[1];
x q[2];
rz(1.1346899) q[3];
sx q[3];
rz(-2.3824771) q[3];
sx q[3];
rz(-1.5838983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7342928) q[2];
sx q[2];
rz(-2.5278842) q[2];
sx q[2];
rz(-1.8180004) q[2];
rz(-1.343441) q[3];
sx q[3];
rz(-2.3517793) q[3];
sx q[3];
rz(0.33671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6650498) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(0.05050412) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(0.49096289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57632724) q[0];
sx q[0];
rz(-2.6748228) q[0];
sx q[0];
rz(-0.76961036) q[0];
rz(0.30740909) q[2];
sx q[2];
rz(-1.4239862) q[2];
sx q[2];
rz(1.4585782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9742187) q[1];
sx q[1];
rz(-2.189296) q[1];
sx q[1];
rz(0.87272404) q[1];
x q[2];
rz(0.36763962) q[3];
sx q[3];
rz(-2.841137) q[3];
sx q[3];
rz(-2.5641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(0.8636221) q[2];
rz(2.9577799) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(2.1553154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.42613906) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(-0.62582985) q[0];
rz(2.4848056) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(1.2084557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55126429) q[0];
sx q[0];
rz(-2.9737824) q[0];
sx q[0];
rz(-2.2237334) q[0];
x q[1];
rz(-0.65932806) q[2];
sx q[2];
rz(-1.314333) q[2];
sx q[2];
rz(1.9224482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0543218) q[1];
sx q[1];
rz(-3.0551214) q[1];
sx q[1];
rz(0.030589624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7478862) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(1.0799112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5586231) q[2];
sx q[2];
rz(-1.8331336) q[2];
sx q[2];
rz(2.6226131) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(1.3435266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2124355) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(2.9875372) q[1];
sx q[1];
rz(-0.93692056) q[1];
sx q[1];
rz(-0.9137203) q[1];
rz(-3.0689756) q[2];
sx q[2];
rz(-2.6238863) q[2];
sx q[2];
rz(-3.1095502) q[2];
rz(-1.9939353) q[3];
sx q[3];
rz(-0.71865766) q[3];
sx q[3];
rz(-1.3620993) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
