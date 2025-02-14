OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80630535) q[0];
sx q[0];
rz(4.181554) q[0];
sx q[0];
rz(9.7747533) q[0];
rz(-1.6236053) q[1];
sx q[1];
rz(-0.94437683) q[1];
sx q[1];
rz(-1.1854393) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703263) q[0];
sx q[0];
rz(-1.2375273) q[0];
sx q[0];
rz(2.6993178) q[0];
rz(2.6529409) q[2];
sx q[2];
rz(-0.46670318) q[2];
sx q[2];
rz(1.8420417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2054641) q[1];
sx q[1];
rz(-1.3359114) q[1];
sx q[1];
rz(-1.0297736) q[1];
rz(0.18326944) q[3];
sx q[3];
rz(-0.96310593) q[3];
sx q[3];
rz(-2.303108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27652937) q[2];
sx q[2];
rz(-0.41215602) q[2];
sx q[2];
rz(-0.40974799) q[2];
rz(-2.3661738) q[3];
sx q[3];
rz(-1.6131468) q[3];
sx q[3];
rz(-2.8342136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85768098) q[0];
sx q[0];
rz(-0.49864054) q[0];
sx q[0];
rz(2.6890802) q[0];
rz(-2.9127938) q[1];
sx q[1];
rz(-2.6401873) q[1];
sx q[1];
rz(-2.9948044) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3892666) q[0];
sx q[0];
rz(-1.379836) q[0];
sx q[0];
rz(1.2021078) q[0];
rz(-pi) q[1];
rz(0.54988523) q[2];
sx q[2];
rz(-2.7061979) q[2];
sx q[2];
rz(2.3741436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9917069) q[1];
sx q[1];
rz(-3.001431) q[1];
sx q[1];
rz(-2.8542942) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5761174) q[3];
sx q[3];
rz(-0.74339044) q[3];
sx q[3];
rz(-0.88228031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9734834) q[2];
sx q[2];
rz(-0.3349458) q[2];
sx q[2];
rz(-2.2701263) q[2];
rz(-2.787309) q[3];
sx q[3];
rz(-0.93577093) q[3];
sx q[3];
rz(1.3080477) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63534921) q[0];
sx q[0];
rz(-0.62017089) q[0];
sx q[0];
rz(-1.1059906) q[0];
rz(-2.1982819) q[1];
sx q[1];
rz(-2.3425075) q[1];
sx q[1];
rz(1.4897289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2882799) q[0];
sx q[0];
rz(-1.5605456) q[0];
sx q[0];
rz(-1.6985189) q[0];
rz(-0.42647894) q[2];
sx q[2];
rz(-0.22828776) q[2];
sx q[2];
rz(-0.84289614) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6506279) q[1];
sx q[1];
rz(-0.59655439) q[1];
sx q[1];
rz(-0.067452879) q[1];
rz(-pi) q[2];
rz(-0.2977887) q[3];
sx q[3];
rz(-2.892916) q[3];
sx q[3];
rz(0.77588785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51365435) q[2];
sx q[2];
rz(-1.1955806) q[2];
sx q[2];
rz(-0.99620831) q[2];
rz(0.57322383) q[3];
sx q[3];
rz(-2.6300391) q[3];
sx q[3];
rz(-2.5073124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117677) q[0];
sx q[0];
rz(-1.2720164) q[0];
sx q[0];
rz(1.5559394) q[0];
rz(-0.32444435) q[1];
sx q[1];
rz(-2.3093846) q[1];
sx q[1];
rz(1.197804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4909013) q[0];
sx q[0];
rz(-1.317121) q[0];
sx q[0];
rz(-2.5462697) q[0];
rz(2.156331) q[2];
sx q[2];
rz(-1.4466431) q[2];
sx q[2];
rz(-0.92701605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3052605) q[1];
sx q[1];
rz(-1.700749) q[1];
sx q[1];
rz(-2.0559681) q[1];
x q[2];
rz(2.9778756) q[3];
sx q[3];
rz(-1.8293899) q[3];
sx q[3];
rz(1.05815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.045362554) q[2];
sx q[2];
rz(-0.54916334) q[2];
sx q[2];
rz(1.1693003) q[2];
rz(0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(-0.72871488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8600334) q[0];
sx q[0];
rz(-0.39880729) q[0];
sx q[0];
rz(-1.1997724) q[0];
rz(1.8343605) q[1];
sx q[1];
rz(-2.2016826) q[1];
sx q[1];
rz(1.7050381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2970151) q[0];
sx q[0];
rz(-2.7488193) q[0];
sx q[0];
rz(1.412084) q[0];
rz(-pi) q[1];
rz(-1.0588516) q[2];
sx q[2];
rz(-0.26341715) q[2];
sx q[2];
rz(2.9204766) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5561394) q[1];
sx q[1];
rz(-1.9227429) q[1];
sx q[1];
rz(-1.8896249) q[1];
x q[2];
rz(2.2186231) q[3];
sx q[3];
rz(-1.9119091) q[3];
sx q[3];
rz(0.34527147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78937078) q[2];
sx q[2];
rz(-1.7033966) q[2];
sx q[2];
rz(-2.5686725) q[2];
rz(-2.2574183) q[3];
sx q[3];
rz(-2.1638162) q[3];
sx q[3];
rz(-2.6463215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5837412) q[0];
sx q[0];
rz(-1.252625) q[0];
sx q[0];
rz(-1.0766693) q[0];
rz(-2.63511) q[1];
sx q[1];
rz(-0.61512893) q[1];
sx q[1];
rz(1.8411676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69715766) q[0];
sx q[0];
rz(-1.5570672) q[0];
sx q[0];
rz(1.5221217) q[0];
x q[1];
rz(-0.22407786) q[2];
sx q[2];
rz(-0.59106088) q[2];
sx q[2];
rz(3.0853809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2357851) q[1];
sx q[1];
rz(-1.7828568) q[1];
sx q[1];
rz(-3.0074988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9229284) q[3];
sx q[3];
rz(-1.7588968) q[3];
sx q[3];
rz(-1.355483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7252245) q[2];
sx q[2];
rz(-2.29795) q[2];
sx q[2];
rz(1.1923403) q[2];
rz(0.013817712) q[3];
sx q[3];
rz(-2.4255224) q[3];
sx q[3];
rz(0.97318399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8927807) q[0];
sx q[0];
rz(-2.2154494) q[0];
sx q[0];
rz(2.7698351) q[0];
rz(-0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(-2.5977792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0293074) q[0];
sx q[0];
rz(-0.72192955) q[0];
sx q[0];
rz(-0.63755905) q[0];
x q[1];
rz(0.45525785) q[2];
sx q[2];
rz(-0.90418679) q[2];
sx q[2];
rz(-0.30991679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2993967) q[1];
sx q[1];
rz(-2.7666515) q[1];
sx q[1];
rz(-3.123466) q[1];
x q[2];
rz(1.722808) q[3];
sx q[3];
rz(-1.4626039) q[3];
sx q[3];
rz(0.044807981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.030293327) q[2];
sx q[2];
rz(-2.0798422) q[2];
sx q[2];
rz(2.566805) q[2];
rz(-0.8542257) q[3];
sx q[3];
rz(-1.4762907) q[3];
sx q[3];
rz(2.0265719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92693555) q[0];
sx q[0];
rz(-2.4605926) q[0];
sx q[0];
rz(-0.32447234) q[0];
rz(0.60943162) q[1];
sx q[1];
rz(-2.29988) q[1];
sx q[1];
rz(0.51469523) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99583944) q[0];
sx q[0];
rz(-2.0059364) q[0];
sx q[0];
rz(0.010220411) q[0];
rz(-pi) q[1];
rz(2.7829172) q[2];
sx q[2];
rz(-2.2513362) q[2];
sx q[2];
rz(-2.2941065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.534651) q[1];
sx q[1];
rz(-1.80503) q[1];
sx q[1];
rz(2.5638084) q[1];
rz(1.8466161) q[3];
sx q[3];
rz(-1.6758741) q[3];
sx q[3];
rz(-0.17227473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9784166) q[2];
sx q[2];
rz(-1.5055483) q[2];
sx q[2];
rz(-2.8669538) q[2];
rz(-2.2091673) q[3];
sx q[3];
rz(-2.7115188) q[3];
sx q[3];
rz(-2.8453804) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68359971) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(0.13548166) q[0];
rz(0.97700351) q[1];
sx q[1];
rz(-0.75825399) q[1];
sx q[1];
rz(3.1405084) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8086804) q[0];
sx q[0];
rz(-1.7263733) q[0];
sx q[0];
rz(-2.6434487) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.094467) q[2];
sx q[2];
rz(-1.4677248) q[2];
sx q[2];
rz(-1.6873941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8813156) q[1];
sx q[1];
rz(-1.0147569) q[1];
sx q[1];
rz(-2.4350321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2204543) q[3];
sx q[3];
rz(-2.2153184) q[3];
sx q[3];
rz(-0.92512475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5486501) q[2];
sx q[2];
rz(-1.0385916) q[2];
sx q[2];
rz(-2.9366117) q[2];
rz(-0.18260469) q[3];
sx q[3];
rz(-2.9349116) q[3];
sx q[3];
rz(0.73777795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36271998) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(2.6642098) q[0];
rz(-3.0248771) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(0.83888549) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29541839) q[0];
sx q[0];
rz(-1.3059761) q[0];
sx q[0];
rz(2.2841262) q[0];
rz(-pi) q[1];
rz(2.7493189) q[2];
sx q[2];
rz(-2.3773411) q[2];
sx q[2];
rz(-1.4832254) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2005208) q[1];
sx q[1];
rz(-1.7571736) q[1];
sx q[1];
rz(-0.69520562) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3319043) q[3];
sx q[3];
rz(-1.5996859) q[3];
sx q[3];
rz(2.1521371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3597151) q[2];
sx q[2];
rz(-2.8350267) q[2];
sx q[2];
rz(-0.26620418) q[2];
rz(-2.6340458) q[3];
sx q[3];
rz(-2.4188953) q[3];
sx q[3];
rz(2.7005196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626641) q[0];
sx q[0];
rz(-1.6829818) q[0];
sx q[0];
rz(-0.80243954) q[0];
rz(0.91213999) q[1];
sx q[1];
rz(-1.9504539) q[1];
sx q[1];
rz(0.84074195) q[1];
rz(-1.4100762) q[2];
sx q[2];
rz(-1.8158603) q[2];
sx q[2];
rz(2.496904) q[2];
rz(-2.38337) q[3];
sx q[3];
rz(-0.52715404) q[3];
sx q[3];
rz(-0.84670443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
