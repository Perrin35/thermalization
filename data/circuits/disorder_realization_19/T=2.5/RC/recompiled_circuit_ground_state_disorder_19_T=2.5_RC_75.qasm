OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(1.0220802) q[1];
sx q[1];
rz(-0.66749579) q[1];
sx q[1];
rz(0.8134841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047377) q[0];
sx q[0];
rz(-1.2053524) q[0];
sx q[0];
rz(-2.9904537) q[0];
x q[1];
rz(2.6935319) q[2];
sx q[2];
rz(-1.8984814) q[2];
sx q[2];
rz(2.0304012) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.95619394) q[1];
sx q[1];
rz(-0.46727249) q[1];
sx q[1];
rz(2.8059792) q[1];
x q[2];
rz(-1.7343821) q[3];
sx q[3];
rz(-1.6601052) q[3];
sx q[3];
rz(-0.79998868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6525314) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(-0.92864621) q[2];
rz(-1.5422025) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(-1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20392513) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(3.0145338) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-2.3703221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64957261) q[0];
sx q[0];
rz(-2.3207773) q[0];
sx q[0];
rz(2.0628906) q[0];
x q[1];
rz(-1.2331687) q[2];
sx q[2];
rz(-1.6281096) q[2];
sx q[2];
rz(-1.7261417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88395547) q[1];
sx q[1];
rz(-1.2246545) q[1];
sx q[1];
rz(-0.37948541) q[1];
rz(-pi) q[2];
rz(0.57621376) q[3];
sx q[3];
rz(-2.0348843) q[3];
sx q[3];
rz(-2.4903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1473006) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(-3.1331983) q[2];
rz(-2.4781135) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405149) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(-0.98835522) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-3.006014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5268742) q[0];
sx q[0];
rz(-1.5784987) q[0];
sx q[0];
rz(-0.016250261) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4246321) q[2];
sx q[2];
rz(-1.8646984) q[2];
sx q[2];
rz(-2.4999121) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6171268) q[1];
sx q[1];
rz(-1.9248665) q[1];
sx q[1];
rz(2.0984142) q[1];
x q[2];
rz(2.9445678) q[3];
sx q[3];
rz(-0.84614119) q[3];
sx q[3];
rz(-2.0565095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93379891) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(-0.03820339) q[2];
rz(-0.52538747) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(0.6853404) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4811089) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(-1.3350217) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(-0.23922051) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804028) q[0];
sx q[0];
rz(-1.0188658) q[0];
sx q[0];
rz(-0.30776382) q[0];
rz(1.1785422) q[2];
sx q[2];
rz(-1.129727) q[2];
sx q[2];
rz(-0.46229306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.071453302) q[1];
sx q[1];
rz(-1.2483368) q[1];
sx q[1];
rz(0.47688213) q[1];
rz(-1.3705105) q[3];
sx q[3];
rz(-1.1173751) q[3];
sx q[3];
rz(-3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(0.0017496721) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.9521451) q[3];
sx q[3];
rz(-0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2413498) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-2.8040366) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(-1.6036124) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46066901) q[0];
sx q[0];
rz(-0.78993714) q[0];
sx q[0];
rz(1.8415175) q[0];
rz(0.29422167) q[2];
sx q[2];
rz(-0.18202848) q[2];
sx q[2];
rz(2.1862669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91572054) q[1];
sx q[1];
rz(-0.32494007) q[1];
sx q[1];
rz(2.48808) q[1];
x q[2];
rz(2.7790623) q[3];
sx q[3];
rz(-0.22386079) q[3];
sx q[3];
rz(0.26765841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43188492) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(-2.0951927) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-2.5939202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980314) q[0];
sx q[0];
rz(-1.2579608) q[0];
sx q[0];
rz(-2.4936254) q[0];
rz(1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(1.3290149) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2752987) q[0];
sx q[0];
rz(-1.4614551) q[0];
sx q[0];
rz(1.7924395) q[0];
rz(-pi) q[1];
rz(-2.8515896) q[2];
sx q[2];
rz(-1.776374) q[2];
sx q[2];
rz(-2.9221591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.700579) q[1];
sx q[1];
rz(-2.3668336) q[1];
sx q[1];
rz(2.0225594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59698128) q[3];
sx q[3];
rz(-2.5425445) q[3];
sx q[3];
rz(-1.057665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7602188) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(2.2729592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182564) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(0.94183952) q[0];
rz(2.9815004) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036445905) q[0];
sx q[0];
rz(-0.22686401) q[0];
sx q[0];
rz(-1.4099717) q[0];
x q[1];
rz(2.0925557) q[2];
sx q[2];
rz(-1.6055487) q[2];
sx q[2];
rz(-0.032973789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6814179) q[1];
sx q[1];
rz(-2.3263086) q[1];
sx q[1];
rz(-0.47487201) q[1];
rz(-1.6113847) q[3];
sx q[3];
rz(-0.49792624) q[3];
sx q[3];
rz(0.17526173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3840702) q[2];
sx q[2];
rz(-1.2382058) q[2];
sx q[2];
rz(-0.43357098) q[2];
rz(1.5844257) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.3442511) q[0];
rz(-1.2035707) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(-1.7787836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31263154) q[0];
sx q[0];
rz(-1.5488708) q[0];
sx q[0];
rz(-1.6005105) q[0];
rz(-0.8261189) q[2];
sx q[2];
rz(-1.4419793) q[2];
sx q[2];
rz(-1.6507738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92719936) q[1];
sx q[1];
rz(-1.3424557) q[1];
sx q[1];
rz(-1.4424999) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8239488) q[3];
sx q[3];
rz(-2.2322725) q[3];
sx q[3];
rz(0.7484439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56889304) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(0.9838689) q[2];
rz(2.900506) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.6702061) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(-1.341691) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-1.0460269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90993308) q[0];
sx q[0];
rz(-2.446736) q[0];
sx q[0];
rz(-0.52082638) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7685031) q[2];
sx q[2];
rz(-1.4914163) q[2];
sx q[2];
rz(0.52620906) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99524263) q[1];
sx q[1];
rz(-2.1533794) q[1];
sx q[1];
rz(1.4928774) q[1];
x q[2];
rz(-0.041954354) q[3];
sx q[3];
rz(-2.0734412) q[3];
sx q[3];
rz(-1.6142538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.86299738) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(-2.7421303) q[2];
rz(-2.5755889) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.852916) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(-0.5823108) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(-0.80642548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18412735) q[0];
sx q[0];
rz(-3.0054207) q[0];
sx q[0];
rz(0.86958142) q[0];
rz(-pi) q[1];
rz(0.66566531) q[2];
sx q[2];
rz(-2.7241926) q[2];
sx q[2];
rz(-1.1102499) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2405673) q[1];
sx q[1];
rz(-1.7738924) q[1];
sx q[1];
rz(2.3939449) q[1];
rz(2.7943198) q[3];
sx q[3];
rz(-1.4662305) q[3];
sx q[3];
rz(-0.73062752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23239423) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(0.8052899) q[2];
rz(-2.9946839) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2522226) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(1.5421142) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(-1.7003822) q[2];
sx q[2];
rz(-2.414188) q[2];
sx q[2];
rz(1.3368171) q[2];
rz(0.44381683) q[3];
sx q[3];
rz(-0.6946105) q[3];
sx q[3];
rz(1.9048922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
