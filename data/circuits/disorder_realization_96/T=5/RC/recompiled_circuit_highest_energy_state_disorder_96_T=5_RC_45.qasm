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
rz(-2.8366635) q[0];
sx q[0];
rz(-3.0753758) q[0];
sx q[0];
rz(2.9856292) q[0];
rz(-1.9744385) q[1];
sx q[1];
rz(-0.65216291) q[1];
sx q[1];
rz(0.48643938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6682566) q[0];
sx q[0];
rz(-1.0205246) q[0];
sx q[0];
rz(2.2815198) q[0];
x q[1];
rz(1.0503642) q[2];
sx q[2];
rz(-1.2005521) q[2];
sx q[2];
rz(-0.56132078) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52594324) q[1];
sx q[1];
rz(-1.2541703) q[1];
sx q[1];
rz(-1.8773882) q[1];
x q[2];
rz(2.3144929) q[3];
sx q[3];
rz(-1.7286282) q[3];
sx q[3];
rz(-0.99770297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8638986) q[2];
sx q[2];
rz(-2.4319477) q[2];
sx q[2];
rz(-0.42897439) q[2];
rz(-3.0947558) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2752537) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(0.66542768) q[0];
rz(2.0003419) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(2.4990987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25099558) q[0];
sx q[0];
rz(-2.3038376) q[0];
sx q[0];
rz(-1.5023853) q[0];
x q[1];
rz(-2.7879509) q[2];
sx q[2];
rz(-2.414915) q[2];
sx q[2];
rz(-2.7893989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4027109) q[1];
sx q[1];
rz(-0.69873525) q[1];
sx q[1];
rz(0.043234392) q[1];
rz(-pi) q[2];
rz(-1.0297433) q[3];
sx q[3];
rz(-1.4192038) q[3];
sx q[3];
rz(-2.9261287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82681727) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(-2.5891506) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(2.4655931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33573547) q[0];
sx q[0];
rz(-0.72000802) q[0];
sx q[0];
rz(-0.82255256) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(1.6792345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8203634) q[0];
sx q[0];
rz(-1.9495533) q[0];
sx q[0];
rz(2.6418153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.64531) q[2];
sx q[2];
rz(-1.9828116) q[2];
sx q[2];
rz(0.81531422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84356704) q[1];
sx q[1];
rz(-1.1192392) q[1];
sx q[1];
rz(0.3275015) q[1];
rz(-pi) q[2];
rz(-0.30841804) q[3];
sx q[3];
rz(-1.2756516) q[3];
sx q[3];
rz(2.4663188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86639577) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(0.59494507) q[2];
rz(2.7368937) q[3];
sx q[3];
rz(-1.1070822) q[3];
sx q[3];
rz(-1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054656) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(-2.8986616) q[0];
rz(-0.88492197) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(3.1387709) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045319917) q[0];
sx q[0];
rz(-2.3191873) q[0];
sx q[0];
rz(-2.1069645) q[0];
rz(-pi) q[1];
rz(2.2523802) q[2];
sx q[2];
rz(-0.74478645) q[2];
sx q[2];
rz(1.7303263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0012008) q[1];
sx q[1];
rz(-2.2691548) q[1];
sx q[1];
rz(-3.0496518) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97798621) q[3];
sx q[3];
rz(-2.1007256) q[3];
sx q[3];
rz(-2.0265614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5367077) q[2];
sx q[2];
rz(-1.9900091) q[2];
sx q[2];
rz(-0.73337698) q[2];
rz(-0.65507656) q[3];
sx q[3];
rz(-0.30787444) q[3];
sx q[3];
rz(-2.5797599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4578399) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(0.14232464) q[0];
rz(-1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(-0.49837643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55492586) q[0];
sx q[0];
rz(-2.2827671) q[0];
sx q[0];
rz(-3.0601383) q[0];
rz(-pi) q[1];
rz(3.0875838) q[2];
sx q[2];
rz(-2.3819792) q[2];
sx q[2];
rz(2.1646966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.315334) q[1];
sx q[1];
rz(-0.79763258) q[1];
sx q[1];
rz(-2.7263097) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33739319) q[3];
sx q[3];
rz(-2.3327347) q[3];
sx q[3];
rz(2.0615426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0087697) q[2];
sx q[2];
rz(-1.2558197) q[2];
sx q[2];
rz(-2.447017) q[2];
rz(-0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(0.25025234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4271127) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(0.66202778) q[0];
rz(1.9561249) q[1];
sx q[1];
rz(-2.1123501) q[1];
sx q[1];
rz(1.9756636) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53138779) q[0];
sx q[0];
rz(-0.90739668) q[0];
sx q[0];
rz(1.8276455) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2189617) q[2];
sx q[2];
rz(-0.6195407) q[2];
sx q[2];
rz(2.9094686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5422029) q[1];
sx q[1];
rz(-1.6205317) q[1];
sx q[1];
rz(1.4982144) q[1];
rz(-pi) q[2];
rz(2.1816129) q[3];
sx q[3];
rz(-2.0997542) q[3];
sx q[3];
rz(0.66416603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7621496) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(-2.4318648) q[2];
rz(-2.659667) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(-0.019439241) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703281) q[0];
sx q[0];
rz(-2.1605991) q[0];
sx q[0];
rz(-1.5671267) q[0];
rz(1.6471242) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(0.87237298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0013638) q[0];
sx q[0];
rz(-1.2622381) q[0];
sx q[0];
rz(0.17805992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.383423) q[2];
sx q[2];
rz(-0.85722605) q[2];
sx q[2];
rz(3.0564412) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9530198) q[1];
sx q[1];
rz(-1.5404048) q[1];
sx q[1];
rz(1.6769483) q[1];
rz(-2.2254785) q[3];
sx q[3];
rz(-1.5903842) q[3];
sx q[3];
rz(0.71924984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81277728) q[2];
sx q[2];
rz(-2.2810563) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(0.40337107) q[3];
sx q[3];
rz(-1.6962467) q[3];
sx q[3];
rz(-0.71322125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23680747) q[0];
sx q[0];
rz(-1.9984364) q[0];
sx q[0];
rz(1.2971725) q[0];
rz(0.5212658) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(3.0788132) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15069467) q[0];
sx q[0];
rz(-2.6932062) q[0];
sx q[0];
rz(1.7240811) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36022236) q[2];
sx q[2];
rz(-1.7884322) q[2];
sx q[2];
rz(-1.9955903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52555841) q[1];
sx q[1];
rz(-1.3241321) q[1];
sx q[1];
rz(-0.14825578) q[1];
rz(-pi) q[2];
rz(-2.6046329) q[3];
sx q[3];
rz(-0.99565047) q[3];
sx q[3];
rz(-1.5123658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89821833) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(0.96837366) q[2];
rz(-2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6043855) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(2.3576417) q[0];
rz(1.2046658) q[1];
sx q[1];
rz(-2.7684863) q[1];
sx q[1];
rz(-1.7663667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7989144) q[0];
sx q[0];
rz(-1.7851789) q[0];
sx q[0];
rz(-0.27002295) q[0];
rz(1.398574) q[2];
sx q[2];
rz(-0.15744124) q[2];
sx q[2];
rz(-2.1724608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3578525) q[1];
sx q[1];
rz(-2.2688818) q[1];
sx q[1];
rz(-0.74905209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0779984) q[3];
sx q[3];
rz(-2.7873899) q[3];
sx q[3];
rz(2.2627047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(-0.44573477) q[2];
rz(-2.5388057) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(-0.84356892) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8730901) q[0];
sx q[0];
rz(-2.6188445) q[0];
sx q[0];
rz(2.6173746) q[0];
rz(1.2007319) q[1];
sx q[1];
rz(-1.2501161) q[1];
sx q[1];
rz(1.9573617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4864593) q[0];
sx q[0];
rz(-2.7441027) q[0];
sx q[0];
rz(-0.50661202) q[0];
rz(3.0734407) q[2];
sx q[2];
rz(-1.8859204) q[2];
sx q[2];
rz(-1.610705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0457234) q[1];
sx q[1];
rz(-1.5830401) q[1];
sx q[1];
rz(0.13282383) q[1];
rz(-1.6053725) q[3];
sx q[3];
rz(-0.94217052) q[3];
sx q[3];
rz(2.2104267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21620096) q[2];
sx q[2];
rz(-2.6915458) q[2];
sx q[2];
rz(1.0518543) q[2];
rz(0.22107407) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(2.2033447) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939519) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(0.37638695) q[1];
sx q[1];
rz(-0.86224894) q[1];
sx q[1];
rz(-1.1663762) q[1];
rz(3.0920269) q[2];
sx q[2];
rz(-2.5709346) q[2];
sx q[2];
rz(1.554629) q[2];
rz(-2.4276499) q[3];
sx q[3];
rz(-1.1957914) q[3];
sx q[3];
rz(1.316432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
