OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82792038) q[0];
sx q[0];
rz(-1.6452687) q[0];
sx q[0];
rz(2.9287455) q[0];
rz(-1.4053474) q[1];
sx q[1];
rz(-1.5264629) q[1];
sx q[1];
rz(2.5636173) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62993193) q[0];
sx q[0];
rz(-1.970147) q[0];
sx q[0];
rz(-2.6005247) q[0];
rz(0.87066381) q[2];
sx q[2];
rz(-1.7144702) q[2];
sx q[2];
rz(-2.0266909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4446775) q[1];
sx q[1];
rz(-1.5648876) q[1];
sx q[1];
rz(-1.9828933) q[1];
rz(0.36823456) q[3];
sx q[3];
rz(-1.4289749) q[3];
sx q[3];
rz(0.54556812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2026244) q[2];
sx q[2];
rz(-1.5105931) q[2];
sx q[2];
rz(2.5228339) q[2];
rz(2.7135811) q[3];
sx q[3];
rz(-1.2330387) q[3];
sx q[3];
rz(1.4324043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0564868) q[0];
sx q[0];
rz(-0.036186941) q[0];
sx q[0];
rz(0.65271839) q[0];
rz(2.808049) q[1];
sx q[1];
rz(-2.3338552) q[1];
sx q[1];
rz(0.52887598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760813) q[0];
sx q[0];
rz(-1.2032857) q[0];
sx q[0];
rz(1.0239081) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1471932) q[2];
sx q[2];
rz(-1.6644239) q[2];
sx q[2];
rz(-0.67568278) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6436504) q[1];
sx q[1];
rz(-1.4406246) q[1];
sx q[1];
rz(1.4243717) q[1];
x q[2];
rz(-2.3102898) q[3];
sx q[3];
rz(-1.5244433) q[3];
sx q[3];
rz(1.6024737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7876106) q[2];
sx q[2];
rz(-1.748961) q[2];
sx q[2];
rz(0.37484136) q[2];
rz(1.2104872) q[3];
sx q[3];
rz(-2.6755302) q[3];
sx q[3];
rz(-1.2208285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1061363) q[0];
sx q[0];
rz(-1.9566256) q[0];
sx q[0];
rz(2.5842066) q[0];
rz(-2.4222971) q[1];
sx q[1];
rz(-0.44984111) q[1];
sx q[1];
rz(2.5296899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61479622) q[0];
sx q[0];
rz(-0.4810534) q[0];
sx q[0];
rz(-0.78562268) q[0];
x q[1];
rz(1.3054661) q[2];
sx q[2];
rz(-2.4396787) q[2];
sx q[2];
rz(0.64218015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0590089) q[1];
sx q[1];
rz(-1.6477668) q[1];
sx q[1];
rz(2.1441035) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7385036) q[3];
sx q[3];
rz(-1.4018019) q[3];
sx q[3];
rz(0.72746659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3376075) q[2];
sx q[2];
rz(-1.3472202) q[2];
sx q[2];
rz(-2.1862629) q[2];
rz(-2.5721926) q[3];
sx q[3];
rz(-1.1207213) q[3];
sx q[3];
rz(1.3620522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3595381) q[0];
sx q[0];
rz(-2.8466917) q[0];
sx q[0];
rz(-0.094757946) q[0];
rz(1.5358198) q[1];
sx q[1];
rz(-1.1631807) q[1];
sx q[1];
rz(-0.43704978) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.248413) q[0];
sx q[0];
rz(-1.6706658) q[0];
sx q[0];
rz(0.084253913) q[0];
x q[1];
rz(-0.37462285) q[2];
sx q[2];
rz(-0.46430507) q[2];
sx q[2];
rz(1.2863867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.028965) q[1];
sx q[1];
rz(-1.8323168) q[1];
sx q[1];
rz(3.1413743) q[1];
x q[2];
rz(0.38744827) q[3];
sx q[3];
rz(-0.94765462) q[3];
sx q[3];
rz(0.98549313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7729418) q[2];
sx q[2];
rz(-2.3029885) q[2];
sx q[2];
rz(1.8013901) q[2];
rz(0.36315) q[3];
sx q[3];
rz(-0.64195091) q[3];
sx q[3];
rz(-2.6974881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-1.7813613) q[0];
sx q[0];
rz(-2.4960127) q[0];
sx q[0];
rz(3.0907104) q[0];
rz(0.52967349) q[1];
sx q[1];
rz(-0.61057225) q[1];
sx q[1];
rz(0.018459056) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7537024) q[0];
sx q[0];
rz(-1.2648996) q[0];
sx q[0];
rz(-2.1056771) q[0];
rz(-1.370212) q[2];
sx q[2];
rz(-2.2772352) q[2];
sx q[2];
rz(2.4141224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5010656) q[1];
sx q[1];
rz(-0.94087761) q[1];
sx q[1];
rz(0.34754158) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8419739) q[3];
sx q[3];
rz(-0.54237759) q[3];
sx q[3];
rz(-0.49132916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20794836) q[2];
sx q[2];
rz(-1.7174145) q[2];
sx q[2];
rz(-1.2733744) q[2];
rz(-0.4452855) q[3];
sx q[3];
rz(-2.4280426) q[3];
sx q[3];
rz(0.85550296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2263366) q[0];
sx q[0];
rz(-2.7610918) q[0];
sx q[0];
rz(-2.8476727) q[0];
rz(1.8982915) q[1];
sx q[1];
rz(-1.8445804) q[1];
sx q[1];
rz(1.4909202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4478233) q[0];
sx q[0];
rz(-0.64720479) q[0];
sx q[0];
rz(1.291841) q[0];
rz(1.6664198) q[2];
sx q[2];
rz(-1.8954512) q[2];
sx q[2];
rz(-2.4839475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.053664899) q[1];
sx q[1];
rz(-1.355989) q[1];
sx q[1];
rz(1.3707036) q[1];
x q[2];
rz(-1.1209247) q[3];
sx q[3];
rz(-0.48512019) q[3];
sx q[3];
rz(-2.7435722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2095498) q[2];
sx q[2];
rz(-1.2806634) q[2];
sx q[2];
rz(3.0397084) q[2];
rz(-2.3515676) q[3];
sx q[3];
rz(-0.70370379) q[3];
sx q[3];
rz(1.2557282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387966) q[0];
sx q[0];
rz(-3.0468472) q[0];
sx q[0];
rz(-2.5139659) q[0];
rz(-1.7812642) q[1];
sx q[1];
rz(-1.5584713) q[1];
sx q[1];
rz(2.9291709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0414114) q[0];
sx q[0];
rz(-0.79754058) q[0];
sx q[0];
rz(1.0486288) q[0];
rz(-pi) q[1];
rz(-1.1721344) q[2];
sx q[2];
rz(-0.63471925) q[2];
sx q[2];
rz(0.29951698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0052988) q[1];
sx q[1];
rz(-1.7583656) q[1];
sx q[1];
rz(1.0636399) q[1];
rz(-pi) q[2];
rz(0.77712975) q[3];
sx q[3];
rz(-2.5087959) q[3];
sx q[3];
rz(1.599904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75954306) q[2];
sx q[2];
rz(-0.86882773) q[2];
sx q[2];
rz(2.7810435) q[2];
rz(-0.77066317) q[3];
sx q[3];
rz(-1.2040851) q[3];
sx q[3];
rz(2.8785021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3537972) q[0];
sx q[0];
rz(-2.0833092) q[0];
sx q[0];
rz(1.298792) q[0];
rz(0.53954387) q[1];
sx q[1];
rz(-1.6863457) q[1];
sx q[1];
rz(1.65666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34299193) q[0];
sx q[0];
rz(-1.2534849) q[0];
sx q[0];
rz(0.57595944) q[0];
rz(-pi) q[1];
rz(-2.6752744) q[2];
sx q[2];
rz(-2.6728485) q[2];
sx q[2];
rz(-3.0006204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6816652) q[1];
sx q[1];
rz(-1.8742537) q[1];
sx q[1];
rz(3.1292679) q[1];
x q[2];
rz(-1.1670145) q[3];
sx q[3];
rz(-2.2444668) q[3];
sx q[3];
rz(-0.63577494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3761042) q[2];
sx q[2];
rz(-1.8643943) q[2];
sx q[2];
rz(2.6831324) q[2];
rz(1.1485398) q[3];
sx q[3];
rz(-0.12980041) q[3];
sx q[3];
rz(0.54500088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5371567) q[0];
sx q[0];
rz(-2.899709) q[0];
sx q[0];
rz(-0.81445527) q[0];
rz(-2.3711329) q[1];
sx q[1];
rz(-1.6043681) q[1];
sx q[1];
rz(1.8155712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0367397) q[0];
sx q[0];
rz(-1.827404) q[0];
sx q[0];
rz(-0.19427257) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6739266) q[2];
sx q[2];
rz(-0.86926341) q[2];
sx q[2];
rz(1.3814545) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0596366) q[1];
sx q[1];
rz(-2.9841712) q[1];
sx q[1];
rz(1.6379959) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2758851) q[3];
sx q[3];
rz(-2.9031347) q[3];
sx q[3];
rz(-1.8973647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0731547) q[2];
sx q[2];
rz(-0.96403733) q[2];
sx q[2];
rz(-2.5835719) q[2];
rz(0.26992118) q[3];
sx q[3];
rz(-2.8840265) q[3];
sx q[3];
rz(-2.423438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.54616) q[0];
sx q[0];
rz(-1.8159001) q[0];
sx q[0];
rz(-2.2823855) q[0];
rz(0.72256207) q[1];
sx q[1];
rz(-2.1911306) q[1];
sx q[1];
rz(-1.7019466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0006726) q[0];
sx q[0];
rz(-1.5790916) q[0];
sx q[0];
rz(-0.0061160864) q[0];
rz(-1.7990803) q[2];
sx q[2];
rz(-1.9365657) q[2];
sx q[2];
rz(1.4073142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7754) q[1];
sx q[1];
rz(-0.79093191) q[1];
sx q[1];
rz(-0.25021199) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4018114) q[3];
sx q[3];
rz(-1.3483182) q[3];
sx q[3];
rz(-2.1232186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.667111) q[2];
sx q[2];
rz(-1.2789395) q[2];
sx q[2];
rz(-1.4943592) q[2];
rz(-0.4161559) q[3];
sx q[3];
rz(-1.6437203) q[3];
sx q[3];
rz(3.0709628) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1304929) q[0];
sx q[0];
rz(-0.53642219) q[0];
sx q[0];
rz(-0.82905967) q[0];
rz(-2.2566055) q[1];
sx q[1];
rz(-2.5288455) q[1];
sx q[1];
rz(1.7987342) q[1];
rz(1.5985684) q[2];
sx q[2];
rz(-1.0226316) q[2];
sx q[2];
rz(1.2050261) q[2];
rz(2.8881843) q[3];
sx q[3];
rz(-0.80066917) q[3];
sx q[3];
rz(-1.0193326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
