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
rz(-1.2008774) q[0];
sx q[0];
rz(5.3826217) q[0];
sx q[0];
rz(9.1917888) q[0];
rz(1.6147344) q[1];
sx q[1];
rz(-1.072071) q[1];
sx q[1];
rz(2.022068) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699864) q[0];
sx q[0];
rz(-1.3363839) q[0];
sx q[0];
rz(3.1283911) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5801611) q[2];
sx q[2];
rz(-1.560964) q[2];
sx q[2];
rz(-1.9468284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0488551) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(-1.2176179) q[1];
rz(-2.5324838) q[3];
sx q[3];
rz(-1.2909596) q[3];
sx q[3];
rz(-2.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69922525) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(0.11288682) q[2];
rz(-2.8804307) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(-1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(-0.054340266) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.4726787) q[1];
sx q[1];
rz(0.36453077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88663) q[0];
sx q[0];
rz(-2.1658362) q[0];
sx q[0];
rz(0.9592077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0436552) q[2];
sx q[2];
rz(-2.4424565) q[2];
sx q[2];
rz(-1.3329891) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6682525) q[1];
sx q[1];
rz(-2.1917731) q[1];
sx q[1];
rz(-1.6644068) q[1];
x q[2];
rz(1.0899312) q[3];
sx q[3];
rz(-0.37853795) q[3];
sx q[3];
rz(2.7906281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(2.0330632) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(1.6746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7473258) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(1.0362097) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(-0.55589693) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5591878) q[0];
sx q[0];
rz(-2.3055861) q[0];
sx q[0];
rz(-0.57484267) q[0];
rz(-pi) q[1];
rz(2.0228407) q[2];
sx q[2];
rz(-0.085918203) q[2];
sx q[2];
rz(-0.8762067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4208918) q[1];
sx q[1];
rz(-1.7016422) q[1];
sx q[1];
rz(0.86549211) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44938748) q[3];
sx q[3];
rz(-1.3419749) q[3];
sx q[3];
rz(-0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7613775) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(1.9492487) q[2];
rz(-3.0377667) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035606774) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(0.22044388) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0649234) q[0];
sx q[0];
rz(-0.60416302) q[0];
sx q[0];
rz(-2.2518065) q[0];
rz(-pi) q[1];
rz(2.1313558) q[2];
sx q[2];
rz(-1.8228442) q[2];
sx q[2];
rz(2.857936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.087763) q[1];
sx q[1];
rz(-1.9265895) q[1];
sx q[1];
rz(-3.0213814) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7679522) q[3];
sx q[3];
rz(-1.9585525) q[3];
sx q[3];
rz(-1.0755634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38349884) q[2];
sx q[2];
rz(-1.3145964) q[2];
sx q[2];
rz(-2.626075) q[2];
rz(0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(-2.3084124) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715217) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(0.66437379) q[0];
rz(-2.6145256) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(-1.1767496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19163469) q[0];
sx q[0];
rz(-1.598888) q[0];
sx q[0];
rz(-1.4856443) q[0];
rz(-pi) q[1];
rz(-1.0254622) q[2];
sx q[2];
rz(-0.76280138) q[2];
sx q[2];
rz(-2.5632586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6709395) q[1];
sx q[1];
rz(-1.3325508) q[1];
sx q[1];
rz(1.5907423) q[1];
x q[2];
rz(2.9620911) q[3];
sx q[3];
rz(-1.0065026) q[3];
sx q[3];
rz(1.5998942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1145757) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(-2.2034755) q[2];
rz(-2.0959057) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.806458) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(-1.8916116) q[0];
rz(2.9373923) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(2.2883889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1565022) q[0];
sx q[0];
rz(-1.9130052) q[0];
sx q[0];
rz(-1.2669105) q[0];
rz(-pi) q[1];
rz(0.77682747) q[2];
sx q[2];
rz(-0.78571253) q[2];
sx q[2];
rz(1.43249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0895929) q[1];
sx q[1];
rz(-1.9885364) q[1];
sx q[1];
rz(-2.3485648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9167658) q[3];
sx q[3];
rz(-2.7802573) q[3];
sx q[3];
rz(1.122235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(-1.9419144) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.6658039) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(0.98261181) q[0];
rz(1.8334552) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(-1.7024202) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5403598) q[0];
sx q[0];
rz(-1.2315016) q[0];
sx q[0];
rz(-1.9215867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7050883) q[2];
sx q[2];
rz(-2.6156313) q[2];
sx q[2];
rz(-2.6924722) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8285259) q[1];
sx q[1];
rz(-1.8983574) q[1];
sx q[1];
rz(2.2799973) q[1];
rz(-pi) q[2];
rz(2.6643326) q[3];
sx q[3];
rz(-1.6972194) q[3];
sx q[3];
rz(1.8074769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-2.2480201) q[2];
sx q[2];
rz(-2.5992375) q[2];
rz(-2.1584623) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.5777609) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(2.3610709) q[0];
rz(1.1753987) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(1.0343879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3607517) q[0];
sx q[0];
rz(-0.409161) q[0];
sx q[0];
rz(-0.43171127) q[0];
rz(1.4277568) q[2];
sx q[2];
rz(-1.489822) q[2];
sx q[2];
rz(1.9097999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0757622) q[1];
sx q[1];
rz(-2.7236863) q[1];
sx q[1];
rz(-0.28553812) q[1];
rz(-pi) q[2];
rz(-2.7195549) q[3];
sx q[3];
rz(-2.5521818) q[3];
sx q[3];
rz(-0.60822076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9913651) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(-3.0312209) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.3519663) q[3];
sx q[3];
rz(0.29236326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.89123911) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(-1.891834) q[0];
rz(-0.091014422) q[1];
sx q[1];
rz(-0.72151557) q[1];
sx q[1];
rz(-2.4299842) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2219951) q[0];
sx q[0];
rz(-1.5775497) q[0];
sx q[0];
rz(0.77829495) q[0];
rz(2.3468412) q[2];
sx q[2];
rz(-2.134179) q[2];
sx q[2];
rz(1.3842441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6775637) q[1];
sx q[1];
rz(-1.11598) q[1];
sx q[1];
rz(-1.9536665) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36104155) q[3];
sx q[3];
rz(-0.78892498) q[3];
sx q[3];
rz(1.6411244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(1.1727772) q[2];
rz(-1.5277537) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(-3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461819) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(1.5661731) q[0];
rz(-0.35789403) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(-0.90248743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76644635) q[0];
sx q[0];
rz(-1.3129108) q[0];
sx q[0];
rz(0.93574406) q[0];
rz(-pi) q[1];
rz(3.112215) q[2];
sx q[2];
rz(-2.4973719) q[2];
sx q[2];
rz(0.51960301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8621626) q[1];
sx q[1];
rz(-1.9050652) q[1];
sx q[1];
rz(2.3390655) q[1];
rz(-pi) q[2];
rz(-1.020458) q[3];
sx q[3];
rz(-1.6452879) q[3];
sx q[3];
rz(2.1410774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1050528) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(1.499929) q[2];
rz(0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70984107) q[0];
sx q[0];
rz(-0.7702282) q[0];
sx q[0];
rz(-1.7256398) q[0];
rz(3.1406828) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(-1.4191237) q[2];
sx q[2];
rz(-0.57874741) q[2];
sx q[2];
rz(-0.49606965) q[2];
rz(2.3933181) q[3];
sx q[3];
rz(-2.1088441) q[3];
sx q[3];
rz(0.90499457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
