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
rz(-3.0012335) q[0];
sx q[0];
rz(-1.4599414) q[0];
sx q[0];
rz(-0.85293823) q[0];
rz(1.9536904) q[1];
sx q[1];
rz(3.3878769) q[1];
sx q[1];
rz(9.090957) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66957314) q[0];
sx q[0];
rz(-1.3206351) q[0];
sx q[0];
rz(-1.0521558) q[0];
x q[1];
rz(2.0707971) q[2];
sx q[2];
rz(-2.1437533) q[2];
sx q[2];
rz(0.9129325) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2961408) q[1];
sx q[1];
rz(-0.23794623) q[1];
sx q[1];
rz(-2.4541357) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0404011) q[3];
sx q[3];
rz(-1.6552166) q[3];
sx q[3];
rz(-2.3197378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9217801) q[2];
sx q[2];
rz(-2.2647936) q[2];
sx q[2];
rz(-0.63966695) q[2];
rz(2.2031247) q[3];
sx q[3];
rz(-0.42249051) q[3];
sx q[3];
rz(1.2708906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6476145) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(1.6139503) q[0];
rz(-0.24213067) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(-2.4193144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5764389) q[0];
sx q[0];
rz(-1.1761888) q[0];
sx q[0];
rz(2.7662686) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26210384) q[2];
sx q[2];
rz(-2.3293709) q[2];
sx q[2];
rz(-1.1103528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4237883) q[1];
sx q[1];
rz(-0.77149888) q[1];
sx q[1];
rz(-0.81591925) q[1];
rz(-pi) q[2];
rz(-0.94808156) q[3];
sx q[3];
rz(-1.4332606) q[3];
sx q[3];
rz(-2.8869473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0904514) q[2];
sx q[2];
rz(-2.800056) q[2];
sx q[2];
rz(-0.086645834) q[2];
rz(0.58049774) q[3];
sx q[3];
rz(-1.9167506) q[3];
sx q[3];
rz(-0.95180029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8101623) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(-1.4587559) q[0];
rz(-0.1156062) q[1];
sx q[1];
rz(-1.9435725) q[1];
sx q[1];
rz(-2.9748532) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519823) q[0];
sx q[0];
rz(-1.7777183) q[0];
sx q[0];
rz(-0.99159411) q[0];
x q[1];
rz(1.2010135) q[2];
sx q[2];
rz(-1.175011) q[2];
sx q[2];
rz(-1.41768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84126842) q[1];
sx q[1];
rz(-2.1144697) q[1];
sx q[1];
rz(0.44692301) q[1];
rz(-pi) q[2];
rz(2.9882713) q[3];
sx q[3];
rz(-2.3905131) q[3];
sx q[3];
rz(-2.4052503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.705767) q[2];
sx q[2];
rz(-2.9220118) q[2];
sx q[2];
rz(1.1337918) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.8688801) q[3];
sx q[3];
rz(-0.72499609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13919203) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(-2.8644526) q[0];
rz(-2.3587522) q[1];
sx q[1];
rz(-1.6354086) q[1];
sx q[1];
rz(0.35071075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6445054) q[0];
sx q[0];
rz(-0.86928029) q[0];
sx q[0];
rz(3.0374466) q[0];
rz(-pi) q[1];
rz(2.8160118) q[2];
sx q[2];
rz(-1.5991889) q[2];
sx q[2];
rz(-0.41159002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8102032) q[1];
sx q[1];
rz(-0.94600979) q[1];
sx q[1];
rz(-2.5532755) q[1];
rz(-pi) q[2];
rz(3.0983885) q[3];
sx q[3];
rz(-1.4625878) q[3];
sx q[3];
rz(-1.4747335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0413401) q[2];
sx q[2];
rz(-1.8625872) q[2];
sx q[2];
rz(1.9913541) q[2];
rz(-2.9669115) q[3];
sx q[3];
rz(-2.8476069) q[3];
sx q[3];
rz(3.0512419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.028285099) q[0];
sx q[0];
rz(-2.568013) q[0];
sx q[0];
rz(-1.1962525) q[0];
rz(0.47736564) q[1];
sx q[1];
rz(-2.3148675) q[1];
sx q[1];
rz(-2.5921519) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7350205) q[0];
sx q[0];
rz(-2.1303249) q[0];
sx q[0];
rz(-0.69243777) q[0];
x q[1];
rz(2.7844593) q[2];
sx q[2];
rz(-0.56740848) q[2];
sx q[2];
rz(-2.0878938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82936278) q[1];
sx q[1];
rz(-2.7155502) q[1];
sx q[1];
rz(-1.8660353) q[1];
x q[2];
rz(-1.1138238) q[3];
sx q[3];
rz(-1.1871871) q[3];
sx q[3];
rz(-2.3758604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4598733) q[2];
sx q[2];
rz(-2.7698066) q[2];
sx q[2];
rz(0.2447153) q[2];
rz(1.6230445) q[3];
sx q[3];
rz(-2.0270429) q[3];
sx q[3];
rz(0.092930704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2731584) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(0.42700818) q[0];
rz(-1.931841) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(3.1337813) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88290993) q[0];
sx q[0];
rz(-2.3091303) q[0];
sx q[0];
rz(-0.29405221) q[0];
rz(-1.5862238) q[2];
sx q[2];
rz(-2.5944105) q[2];
sx q[2];
rz(1.5146966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58578929) q[1];
sx q[1];
rz(-2.2206056) q[1];
sx q[1];
rz(1.6342249) q[1];
rz(0.2711556) q[3];
sx q[3];
rz(-2.249345) q[3];
sx q[3];
rz(-0.0038879768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4012287) q[2];
sx q[2];
rz(-2.2682891) q[2];
sx q[2];
rz(2.002423) q[2];
rz(0.27979699) q[3];
sx q[3];
rz(-1.3236902) q[3];
sx q[3];
rz(1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4918268) q[0];
sx q[0];
rz(-2.7597646) q[0];
sx q[0];
rz(-2.5873798) q[0];
rz(-3.0795433) q[1];
sx q[1];
rz(-0.65826145) q[1];
sx q[1];
rz(-0.87472349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59959182) q[0];
sx q[0];
rz(-2.3440147) q[0];
sx q[0];
rz(-0.16982689) q[0];
x q[1];
rz(-1.563233) q[2];
sx q[2];
rz(-1.5763792) q[2];
sx q[2];
rz(-1.9703608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3372015) q[1];
sx q[1];
rz(-2.534165) q[1];
sx q[1];
rz(1.1056221) q[1];
rz(2.234972) q[3];
sx q[3];
rz(-1.8795089) q[3];
sx q[3];
rz(-1.8083835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41500652) q[2];
sx q[2];
rz(-2.1241302) q[2];
sx q[2];
rz(-2.2028108) q[2];
rz(-0.69432652) q[3];
sx q[3];
rz(-0.66803473) q[3];
sx q[3];
rz(-0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0701377) q[0];
sx q[0];
rz(-0.55971611) q[0];
sx q[0];
rz(2.8330084) q[0];
rz(1.4831108) q[1];
sx q[1];
rz(-1.7959271) q[1];
sx q[1];
rz(2.9187091) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6312941) q[0];
sx q[0];
rz(-1.1985072) q[0];
sx q[0];
rz(-0.75081236) q[0];
rz(-pi) q[1];
rz(0.42474681) q[2];
sx q[2];
rz(-2.8966689) q[2];
sx q[2];
rz(-0.67451678) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8376779) q[1];
sx q[1];
rz(-0.24824809) q[1];
sx q[1];
rz(1.7294238) q[1];
rz(-0.65309955) q[3];
sx q[3];
rz(-0.79721236) q[3];
sx q[3];
rz(-1.2739158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3234723) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(-2.4774614) q[2];
rz(-1.1431665) q[3];
sx q[3];
rz(-1.4799456) q[3];
sx q[3];
rz(-2.0460879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42300647) q[0];
sx q[0];
rz(-1.9976595) q[0];
sx q[0];
rz(-0.018420694) q[0];
rz(0.1704692) q[1];
sx q[1];
rz(-1.8959931) q[1];
sx q[1];
rz(0.19759321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9182943) q[0];
sx q[0];
rz(-0.16897783) q[0];
sx q[0];
rz(-1.6264982) q[0];
x q[1];
rz(-0.05240106) q[2];
sx q[2];
rz(-1.2049071) q[2];
sx q[2];
rz(3.0405557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8379544) q[1];
sx q[1];
rz(-1.8686864) q[1];
sx q[1];
rz(2.9998051) q[1];
rz(0.13473265) q[3];
sx q[3];
rz(-2.7425457) q[3];
sx q[3];
rz(-1.4302554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0240872) q[2];
sx q[2];
rz(-2.0427637) q[2];
sx q[2];
rz(-0.39829028) q[2];
rz(0.5106709) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(0.85810703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9061822) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(2.9421575) q[0];
rz(-2.8681352) q[1];
sx q[1];
rz(-2.7077935) q[1];
sx q[1];
rz(1.010703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9977545) q[0];
sx q[0];
rz(-1.8324277) q[0];
sx q[0];
rz(2.4701443) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9062073) q[2];
sx q[2];
rz(-1.5512677) q[2];
sx q[2];
rz(0.69695401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6124552) q[1];
sx q[1];
rz(-2.3061064) q[1];
sx q[1];
rz(1.6409954) q[1];
x q[2];
rz(1.2730678) q[3];
sx q[3];
rz(-1.862251) q[3];
sx q[3];
rz(0.028226559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.345674) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(-2.3460713) q[2];
rz(2.7703721) q[3];
sx q[3];
rz(-1.755654) q[3];
sx q[3];
rz(0.25446874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.026074792) q[0];
sx q[0];
rz(-1.7255029) q[0];
sx q[0];
rz(1.2988476) q[0];
rz(2.5904291) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(1.7953095) q[2];
sx q[2];
rz(-2.2835352) q[2];
sx q[2];
rz(2.78751) q[2];
rz(-2.9007111) q[3];
sx q[3];
rz(-1.1987562) q[3];
sx q[3];
rz(-0.22507122) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
