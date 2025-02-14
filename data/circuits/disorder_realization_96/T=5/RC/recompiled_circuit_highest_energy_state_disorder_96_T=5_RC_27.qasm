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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(0.15596341) q[0];
rz(1.1671542) q[1];
sx q[1];
rz(3.7937556) q[1];
sx q[1];
rz(8.9383386) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64336813) q[0];
sx q[0];
rz(-0.8684477) q[0];
sx q[0];
rz(-0.81612103) q[0];
rz(-pi) q[1];
rz(2.2335792) q[2];
sx q[2];
rz(-2.5129299) q[2];
sx q[2];
rz(-2.6952621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52594324) q[1];
sx q[1];
rz(-1.8874223) q[1];
sx q[1];
rz(-1.2642045) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3144929) q[3];
sx q[3];
rz(-1.7286282) q[3];
sx q[3];
rz(0.99770297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8638986) q[2];
sx q[2];
rz(-2.4319477) q[2];
sx q[2];
rz(0.42897439) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2752537) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(-2.476165) q[0];
rz(-1.1412507) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(-0.64249396) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25099558) q[0];
sx q[0];
rz(-0.83775508) q[0];
sx q[0];
rz(1.5023853) q[0];
x q[1];
rz(1.8694473) q[2];
sx q[2];
rz(-0.89787102) q[2];
sx q[2];
rz(-2.330614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4027109) q[1];
sx q[1];
rz(-0.69873525) q[1];
sx q[1];
rz(-3.0983583) q[1];
rz(-pi) q[2];
rz(1.2824545) q[3];
sx q[3];
rz(-0.55984646) q[3];
sx q[3];
rz(1.1091055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3147754) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(2.5891506) q[2];
rz(-1.4779444) q[3];
sx q[3];
rz(-1.2976357) q[3];
sx q[3];
rz(-0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8058572) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(2.3190401) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(1.4623581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0912446) q[0];
sx q[0];
rz(-2.0322662) q[0];
sx q[0];
rz(-1.9965003) q[0];
rz(-pi) q[1];
rz(-0.41303582) q[2];
sx q[2];
rz(-1.6390642) q[2];
sx q[2];
rz(-2.3562246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2980256) q[1];
sx q[1];
rz(-2.0223534) q[1];
sx q[1];
rz(0.3275015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8796647) q[3];
sx q[3];
rz(-1.2761242) q[3];
sx q[3];
rz(-2.1536649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86639577) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(2.5466476) q[2];
rz(2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(-1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054656) q[0];
sx q[0];
rz(-1.2768224) q[0];
sx q[0];
rz(2.8986616) q[0];
rz(2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(-3.1387709) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3784175) q[0];
sx q[0];
rz(-2.2522914) q[0];
sx q[0];
rz(0.50294089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88921241) q[2];
sx q[2];
rz(-0.74478645) q[2];
sx q[2];
rz(-1.7303263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8587715) q[1];
sx q[1];
rz(-2.438218) q[1];
sx q[1];
rz(1.6797296) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97798621) q[3];
sx q[3];
rz(-2.1007256) q[3];
sx q[3];
rz(1.1150313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(-0.73337698) q[2];
rz(2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68375278) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(2.999268) q[0];
rz(-1.7798452) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(-2.6432162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866668) q[0];
sx q[0];
rz(-0.85882551) q[0];
sx q[0];
rz(-0.081454309) q[0];
rz(-pi) q[1];
rz(-1.6220197) q[2];
sx q[2];
rz(-0.81256676) q[2];
sx q[2];
rz(-1.0513154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.315334) q[1];
sx q[1];
rz(-2.3439601) q[1];
sx q[1];
rz(-2.7263097) q[1];
rz(-2.8041995) q[3];
sx q[3];
rz(-0.80885799) q[3];
sx q[3];
rz(1.0800501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.132823) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(-0.69457561) q[2];
rz(0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-0.25025234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4271127) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(2.4795649) q[0];
rz(1.9561249) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(-1.9756636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2625354) q[0];
sx q[0];
rz(-1.7723119) q[0];
sx q[0];
rz(0.67964566) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4170961) q[2];
sx q[2];
rz(-0.96818334) q[2];
sx q[2];
rz(-3.1069047) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97502134) q[1];
sx q[1];
rz(-1.4983043) q[1];
sx q[1];
rz(-3.0917262) q[1];
x q[2];
rz(-2.1816129) q[3];
sx q[3];
rz(-2.0997542) q[3];
sx q[3];
rz(2.4774266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.379443) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(0.70972788) q[2];
rz(-0.48192561) q[3];
sx q[3];
rz(-2.66633) q[3];
sx q[3];
rz(3.1221534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67126453) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(-1.5671267) q[0];
rz(-1.6471242) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(2.2692197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013638) q[0];
sx q[0];
rz(-1.2622381) q[0];
sx q[0];
rz(-0.17805992) q[0];
rz(0.21193223) q[2];
sx q[2];
rz(-2.4080347) q[2];
sx q[2];
rz(2.9447945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9530198) q[1];
sx q[1];
rz(-1.6011878) q[1];
sx q[1];
rz(1.6769483) q[1];
rz(-pi) q[2];
rz(-0.91611417) q[3];
sx q[3];
rz(-1.5512084) q[3];
sx q[3];
rz(-2.4223428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81277728) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(-2.7382216) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-2.4283714) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23680747) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(-1.2971725) q[0];
rz(2.6203268) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(-3.0788132) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2817665) q[0];
sx q[0];
rz(-1.504557) q[0];
sx q[0];
rz(-2.0145922) q[0];
x q[1];
rz(0.36022236) q[2];
sx q[2];
rz(-1.3531605) q[2];
sx q[2];
rz(1.1460024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.023397327) q[1];
sx q[1];
rz(-0.28701008) q[1];
sx q[1];
rz(-2.1013409) q[1];
x q[2];
rz(0.90274685) q[3];
sx q[3];
rz(-2.3760736) q[3];
sx q[3];
rz(-0.79878858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2433743) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(-2.173219) q[2];
rz(-0.51356703) q[3];
sx q[3];
rz(-1.7889203) q[3];
sx q[3];
rz(2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6043855) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(0.78395098) q[0];
rz(1.2046658) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(1.7663667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3426783) q[0];
sx q[0];
rz(-1.3564137) q[0];
sx q[0];
rz(-2.8715697) q[0];
x q[1];
rz(1.4156467) q[2];
sx q[2];
rz(-1.5439234) q[2];
sx q[2];
rz(2.7100615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32574496) q[1];
sx q[1];
rz(-1.0220075) q[1];
sx q[1];
rz(-2.4239848) q[1];
rz(1.0635942) q[3];
sx q[3];
rz(-2.7873899) q[3];
sx q[3];
rz(0.87888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2969926) q[2];
sx q[2];
rz(-2.832909) q[2];
sx q[2];
rz(0.44573477) q[2];
rz(-0.60278696) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(0.84356892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8730901) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(2.6173746) q[0];
rz(1.2007319) q[1];
sx q[1];
rz(-1.2501161) q[1];
sx q[1];
rz(-1.1842309) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94468695) q[0];
sx q[0];
rz(-1.9160998) q[0];
sx q[0];
rz(-1.3698335) q[0];
rz(0.068151926) q[2];
sx q[2];
rz(-1.8859204) q[2];
sx q[2];
rz(-1.5308876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47329119) q[1];
sx q[1];
rz(-1.4379825) q[1];
sx q[1];
rz(-1.5831489) q[1];
x q[2];
rz(-3.0940787) q[3];
sx q[3];
rz(-0.62944747) q[3];
sx q[3];
rz(2.269182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9253917) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(2.0897384) q[2];
rz(-2.9205186) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(-0.93824798) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939519) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(-0.37638695) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(2.5714931) q[2];
sx q[2];
rz(-1.5975633) q[2];
sx q[2];
rz(-0.057889197) q[2];
rz(1.0906278) q[3];
sx q[3];
rz(-0.91560293) q[3];
sx q[3];
rz(-0.56165725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
