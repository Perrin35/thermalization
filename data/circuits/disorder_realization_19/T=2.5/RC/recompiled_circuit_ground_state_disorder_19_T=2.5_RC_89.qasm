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
rz(1.5751155) q[0];
sx q[0];
rz(10.549904) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4018702) q[0];
sx q[0];
rz(-0.39416322) q[0];
sx q[0];
rz(1.1959082) q[0];
rz(-2.4762597) q[2];
sx q[2];
rz(-2.5931381) q[2];
sx q[2];
rz(-1.0496333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95619394) q[1];
sx q[1];
rz(-0.46727249) q[1];
sx q[1];
rz(-2.8059792) q[1];
rz(-0.090510719) q[3];
sx q[3];
rz(-1.4078684) q[3];
sx q[3];
rz(0.75608692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6525314) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(-2.2129464) q[2];
rz(-1.5422025) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20392513) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(0.12705886) q[0];
rz(-0.98310414) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(0.7712706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5707023) q[0];
sx q[0];
rz(-1.2178019) q[0];
sx q[0];
rz(2.3285026) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7423058) q[2];
sx q[2];
rz(-2.7993188) q[2];
sx q[2];
rz(2.8245087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88395547) q[1];
sx q[1];
rz(-1.2246545) q[1];
sx q[1];
rz(2.7621072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5653789) q[3];
sx q[3];
rz(-2.0348843) q[3];
sx q[3];
rz(2.4903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9942921) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10107772) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(0.98835522) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-0.13557869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39865935) q[0];
sx q[0];
rz(-3.1236095) q[0];
sx q[0];
rz(0.44264408) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71696059) q[2];
sx q[2];
rz(-1.2768942) q[2];
sx q[2];
rz(-2.4999121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49040321) q[1];
sx q[1];
rz(-0.62593834) q[1];
sx q[1];
rz(2.2041026) q[1];
x q[2];
rz(-0.83643416) q[3];
sx q[3];
rz(-1.7178917) q[3];
sx q[3];
rz(0.35415748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2077937) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(-0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66048375) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(2.5307181) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(2.9023721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67475457) q[0];
sx q[0];
rz(-1.8316937) q[0];
sx q[0];
rz(-0.99715085) q[0];
rz(1.1785422) q[2];
sx q[2];
rz(-1.129727) q[2];
sx q[2];
rz(-0.46229306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4799616) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(-1.2110787) q[1];
x q[2];
rz(2.753965) q[3];
sx q[3];
rz(-0.49285965) q[3];
sx q[3];
rz(-0.29407497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7188344) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(-3.139843) q[2];
rz(0.29378978) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(-2.8747115) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2413498) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(0.84306651) q[0];
rz(-2.8040366) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(1.5379803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8385411) q[0];
sx q[0];
rz(-1.3796796) q[0];
sx q[0];
rz(0.79940807) q[0];
rz(-pi) q[1];
rz(2.9672253) q[2];
sx q[2];
rz(-1.5182759) q[2];
sx q[2];
rz(-0.32584056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5953491) q[1];
sx q[1];
rz(-1.3145295) q[1];
sx q[1];
rz(1.3687737) q[1];
x q[2];
rz(-0.36253039) q[3];
sx q[3];
rz(-2.9177319) q[3];
sx q[3];
rz(-0.26765841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43188492) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(2.0951927) q[2];
rz(-2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(2.5939202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(-2.4936254) q[0];
rz(-1.7421534) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(1.3290149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15547046) q[0];
sx q[0];
rz(-0.24674812) q[0];
sx q[0];
rz(2.0339436) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3565265) q[2];
sx q[2];
rz(-1.2870711) q[2];
sx q[2];
rz(-1.8510712) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4635515) q[1];
sx q[1];
rz(-1.8811418) q[1];
sx q[1];
rz(2.2928659) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2043318) q[3];
sx q[3];
rz(-1.0856589) q[3];
sx q[3];
rz(-1.7462891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7602188) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(3.0717487) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182564) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(-0.94183952) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401682) q[0];
sx q[0];
rz(-1.7946825) q[0];
sx q[0];
rz(3.1046449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6404387) q[2];
sx q[2];
rz(-0.52280871) q[2];
sx q[2];
rz(-1.4774587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3698764) q[1];
sx q[1];
rz(-1.9100921) q[1];
sx q[1];
rz(-2.3849065) q[1];
rz(-pi) q[2];
rz(1.6113847) q[3];
sx q[3];
rz(-0.49792624) q[3];
sx q[3];
rz(-0.17526173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3840702) q[2];
sx q[2];
rz(-1.2382058) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(1.5571669) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.7973416) q[0];
rz(1.2035707) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(1.7787836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2479073) q[0];
sx q[0];
rz(-3.1046668) q[0];
sx q[0];
rz(-0.93494995) q[0];
x q[1];
rz(0.8261189) q[2];
sx q[2];
rz(-1.4419793) q[2];
sx q[2];
rz(1.6507738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4688022) q[1];
sx q[1];
rz(-1.6957449) q[1];
sx q[1];
rz(2.9114257) q[1];
rz(-2.8239488) q[3];
sx q[3];
rz(-0.90932019) q[3];
sx q[3];
rz(2.3931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56889304) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(-0.9838689) q[2];
rz(-2.900506) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702061) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(-0.27467003) q[0];
rz(-1.7999016) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(2.0955657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5902025) q[0];
sx q[0];
rz(-2.1596163) q[0];
sx q[0];
rz(-1.9639652) q[0];
rz(-pi) q[1];
rz(1.1860023) q[2];
sx q[2];
rz(-0.21285393) q[2];
sx q[2];
rz(1.720088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6089692) q[1];
sx q[1];
rz(-1.6358422) q[1];
sx q[1];
rz(2.5576127) q[1];
rz(-pi) q[2];
rz(1.4946497) q[3];
sx q[3];
rz(-0.5042432) q[3];
sx q[3];
rz(-1.7011736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(-2.7421303) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852916) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-2.5592819) q[0];
rz(2.9583926) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(2.3351672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2516625) q[0];
sx q[0];
rz(-1.6747054) q[0];
sx q[0];
rz(-0.088168747) q[0];
x q[1];
rz(1.3034794) q[2];
sx q[2];
rz(-1.8952994) q[2];
sx q[2];
rz(-0.40057202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88397056) q[1];
sx q[1];
rz(-2.3720354) q[1];
sx q[1];
rz(0.29410024) q[1];
rz(2.8424758) q[3];
sx q[3];
rz(-0.36206216) q[3];
sx q[3];
rz(-2.5821834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(-2.3363028) q[2];
rz(0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.88937) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(1.5421142) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(2.2940214) q[2];
sx q[2];
rz(-1.4847652) q[2];
sx q[2];
rz(2.8105856) q[2];
rz(-0.64503786) q[3];
sx q[3];
rz(-1.8492263) q[3];
sx q[3];
rz(3.1254569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
