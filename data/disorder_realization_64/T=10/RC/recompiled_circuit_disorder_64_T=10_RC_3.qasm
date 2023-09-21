OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(1.8145476) q[1];
sx q[1];
rz(10.756455) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(-2.1190686) q[0];
rz(-pi) q[1];
rz(-1.9900436) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-0.56564769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3803346) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(-0.18591979) q[1];
x q[2];
rz(-2.0256151) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(-1.7636553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(2.0781793) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20365276) q[0];
sx q[0];
rz(-1.5601888) q[0];
sx q[0];
rz(-1.525735) q[0];
rz(-pi) q[1];
rz(-0.64237853) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(0.013052879) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9284294) q[1];
sx q[1];
rz(-0.82248102) q[1];
sx q[1];
rz(-0.33555062) q[1];
x q[2];
rz(-0.1699154) q[3];
sx q[3];
rz(-1.3193469) q[3];
sx q[3];
rz(-0.91472317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(1.089383) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3141146) q[0];
sx q[0];
rz(-2.4347322) q[0];
sx q[0];
rz(2.2543738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9956279) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(-0.59031634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3315862) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(-1.7942795) q[1];
rz(-pi) q[2];
rz(-1.009922) q[3];
sx q[3];
rz(-0.86285931) q[3];
sx q[3];
rz(-0.55299711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4041876) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(-0.23097158) q[0];
rz(2.9927822) q[2];
sx q[2];
rz(-2.4850922) q[2];
sx q[2];
rz(-1.3615001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6105146) q[1];
sx q[1];
rz(-2.5602166) q[1];
sx q[1];
rz(0.73552144) q[1];
rz(-pi) q[2];
rz(2.8344645) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(0.083837282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-1.0353154) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028776289) q[0];
sx q[0];
rz(-1.1338286) q[0];
sx q[0];
rz(2.2321738) q[0];
rz(-pi) q[1];
rz(0.073280235) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(-2.5626593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(-0.095468949) q[1];
x q[2];
rz(0.27541311) q[3];
sx q[3];
rz(-1.7731035) q[3];
sx q[3];
rz(-1.3261258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.996421) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3375895) q[0];
sx q[0];
rz(-2.1875256) q[0];
sx q[0];
rz(-1.7417275) q[0];
rz(-pi) q[1];
rz(-1.1057165) q[2];
sx q[2];
rz(-2.1207223) q[2];
sx q[2];
rz(-1.20649) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0761557) q[1];
sx q[1];
rz(-1.3216615) q[1];
sx q[1];
rz(2.8912828) q[1];
rz(-pi) q[2];
rz(-0.16222555) q[3];
sx q[3];
rz(-0.82419206) q[3];
sx q[3];
rz(1.970286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4371944) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(2.722446) q[0];
rz(-1.8998434) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(-0.74908756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-2.8357361) q[1];
rz(-pi) q[2];
rz(-0.52280207) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(0.83731246) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(2.1829139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4184957) q[0];
sx q[0];
rz(-0.76508689) q[0];
sx q[0];
rz(-2.3575248) q[0];
rz(-pi) q[1];
rz(1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(1.4657071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9194591) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(-1.3608576) q[1];
rz(2.557425) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(-1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-2.6760496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(2.4337016) q[0];
x q[1];
rz(1.2836254) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(-0.62921333) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1381582) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(2.483063) q[1];
rz(-pi) q[2];
rz(0.97165473) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(-3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(2.5591992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84401417) q[0];
sx q[0];
rz(-1.3676924) q[0];
sx q[0];
rz(2.2556979) q[0];
x q[1];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.7047802) q[2];
sx q[2];
rz(-1.6756563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6813587) q[1];
sx q[1];
rz(-1.453062) q[1];
sx q[1];
rz(-1.7540936) q[1];
x q[2];
rz(1.598761) q[3];
sx q[3];
rz(-1.7879322) q[3];
sx q[3];
rz(-1.0167227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.3394042) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(2.3274158) q[2];
sx q[2];
rz(-1.7780751) q[2];
sx q[2];
rz(1.8572394) q[2];
rz(-2.3253757) q[3];
sx q[3];
rz(-0.64707884) q[3];
sx q[3];
rz(-0.94221471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];