OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3759484) q[0];
sx q[0];
rz(3.5557616) q[0];
sx q[0];
rz(10.226305) q[0];
rz(3.1385359) q[1];
sx q[1];
rz(-2.2771775) q[1];
sx q[1];
rz(-3.0473696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202336) q[0];
sx q[0];
rz(-0.84478837) q[0];
sx q[0];
rz(2.0840706) q[0];
rz(-pi) q[1];
rz(3.1108833) q[2];
sx q[2];
rz(-1.2982663) q[2];
sx q[2];
rz(2.3365793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-2.5818425) q[1];
sx q[1];
rz(-1.00882) q[1];
rz(-pi) q[2];
rz(-0.15259907) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(-1.847037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6039383) q[2];
sx q[2];
rz(-0.64017576) q[2];
sx q[2];
rz(1.5113277) q[2];
rz(2.9553735) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(-2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6179825) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(-0.13667662) q[0];
rz(-3.0467721) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(0.76900855) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.757763) q[0];
sx q[0];
rz(-1.7065598) q[0];
sx q[0];
rz(2.7635937) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4479577) q[2];
sx q[2];
rz(-1.8497548) q[2];
sx q[2];
rz(0.50194937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8254346) q[1];
sx q[1];
rz(-0.76091754) q[1];
sx q[1];
rz(-0.23657152) q[1];
rz(-1.8212998) q[3];
sx q[3];
rz(-3.0022096) q[3];
sx q[3];
rz(1.2047307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39917055) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(1.7102309) q[2];
rz(2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2630149) q[0];
sx q[0];
rz(-0.37913015) q[0];
sx q[0];
rz(1.0947134) q[0];
rz(-1.3021944) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(3.1375695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510295) q[0];
sx q[0];
rz(-0.73712611) q[0];
sx q[0];
rz(-1.766979) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3962069) q[2];
sx q[2];
rz(-2.2507387) q[2];
sx q[2];
rz(-0.58733515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.04847542) q[1];
sx q[1];
rz(-0.25893908) q[1];
sx q[1];
rz(-3.0228651) q[1];
x q[2];
rz(-0.28272776) q[3];
sx q[3];
rz(-0.73388956) q[3];
sx q[3];
rz(-0.80724387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5502988) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(0.44814056) q[2];
rz(2.6552933) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(-1.1264616) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9124741) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(1.9122596) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(-0.545048) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6908643) q[0];
sx q[0];
rz(-2.0929027) q[0];
sx q[0];
rz(-2.1073116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2322234) q[2];
sx q[2];
rz(-2.0880359) q[2];
sx q[2];
rz(-1.478372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5601756) q[1];
sx q[1];
rz(-2.9184074) q[1];
sx q[1];
rz(-1.8471922) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6757984) q[3];
sx q[3];
rz(-2.287103) q[3];
sx q[3];
rz(-2.4972625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58686078) q[2];
sx q[2];
rz(-1.502159) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(-1.933291) q[3];
sx q[3];
rz(-2.8390563) q[3];
sx q[3];
rz(-1.6944073) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(-0.30715352) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-0.53877962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8545495) q[0];
sx q[0];
rz(-1.7998621) q[0];
sx q[0];
rz(1.4255217) q[0];
x q[1];
rz(-0.42947763) q[2];
sx q[2];
rz(-1.6971849) q[2];
sx q[2];
rz(2.1548558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6314124) q[1];
sx q[1];
rz(-2.4459527) q[1];
sx q[1];
rz(2.5744777) q[1];
x q[2];
rz(2.6461824) q[3];
sx q[3];
rz(-2.3461079) q[3];
sx q[3];
rz(-1.0496248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4469845) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(-1.851932) q[2];
rz(1.5223711) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(1.6600018) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30037844) q[0];
sx q[0];
rz(-0.9391681) q[0];
sx q[0];
rz(3.0389431) q[0];
rz(2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(-2.5644459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6711376) q[0];
sx q[0];
rz(-1.5744493) q[0];
sx q[0];
rz(1.5994344) q[0];
rz(2.8961298) q[2];
sx q[2];
rz(-2.5405014) q[2];
sx q[2];
rz(-0.18537755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53458339) q[1];
sx q[1];
rz(-1.4614236) q[1];
sx q[1];
rz(-1.6499728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8982557) q[3];
sx q[3];
rz(-0.84210448) q[3];
sx q[3];
rz(-3.1407243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.131549) q[2];
sx q[2];
rz(-0.71655822) q[2];
sx q[2];
rz(-1.6183759) q[2];
rz(3.0736382) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(-0.84097356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1286569) q[0];
sx q[0];
rz(-2.1067297) q[0];
sx q[0];
rz(2.3701684) q[0];
rz(-2.6029288) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(2.0753986) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2832269) q[0];
sx q[0];
rz(-1.3182501) q[0];
sx q[0];
rz(-3.0622803) q[0];
rz(-pi) q[1];
rz(0.23287878) q[2];
sx q[2];
rz(-1.0693197) q[2];
sx q[2];
rz(1.5053144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8437517) q[1];
sx q[1];
rz(-0.21756309) q[1];
sx q[1];
rz(3.080653) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4769884) q[3];
sx q[3];
rz(-2.7258337) q[3];
sx q[3];
rz(-2.1615054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6923339) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(0.44245693) q[2];
rz(-2.9509406) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(2.382522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2451179) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(-1.3099439) q[0];
rz(-2.7438121) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(-1.7479053) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1398449) q[0];
sx q[0];
rz(-2.4187517) q[0];
sx q[0];
rz(0.87840898) q[0];
x q[1];
rz(2.2289235) q[2];
sx q[2];
rz(-2.7670112) q[2];
sx q[2];
rz(-0.41824579) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19119769) q[1];
sx q[1];
rz(-1.1069595) q[1];
sx q[1];
rz(1.9056803) q[1];
x q[2];
rz(2.05768) q[3];
sx q[3];
rz(-0.22398914) q[3];
sx q[3];
rz(1.6897702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(-0.69463378) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6756814) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(-0.95941108) q[0];
rz(-2.361946) q[1];
sx q[1];
rz(-1.0028853) q[1];
sx q[1];
rz(2.8693105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355541) q[0];
sx q[0];
rz(-1.5840496) q[0];
sx q[0];
rz(-0.0051965836) q[0];
rz(-pi) q[1];
rz(-1.9302888) q[2];
sx q[2];
rz(-2.4036416) q[2];
sx q[2];
rz(2.8518845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4466557) q[1];
sx q[1];
rz(-0.65498252) q[1];
sx q[1];
rz(2.2769364) q[1];
rz(-pi) q[2];
rz(-2.168448) q[3];
sx q[3];
rz(-2.2242507) q[3];
sx q[3];
rz(-2.0367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4043364) q[2];
sx q[2];
rz(-1.6360838) q[2];
sx q[2];
rz(0.20757248) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-2.6360377) q[3];
sx q[3];
rz(1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614457) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(-2.0848059) q[0];
rz(-2.5959004) q[1];
sx q[1];
rz(-2.1634407) q[1];
sx q[1];
rz(-0.32870764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0595039) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(-2.2057981) q[0];
rz(0.64735554) q[2];
sx q[2];
rz(-3.050905) q[2];
sx q[2];
rz(-0.91435223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.199177) q[1];
sx q[1];
rz(-2.7991382) q[1];
sx q[1];
rz(1.6204932) q[1];
rz(-pi) q[2];
rz(2.3700033) q[3];
sx q[3];
rz(-1.7696915) q[3];
sx q[3];
rz(1.0071013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9823965) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(2.0170434) q[2];
rz(-0.076796181) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(2.2176149) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9818253) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(2.3962085) q[2];
sx q[2];
rz(-1.1911662) q[2];
sx q[2];
rz(0.62698812) q[2];
rz(-2.7362551) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
