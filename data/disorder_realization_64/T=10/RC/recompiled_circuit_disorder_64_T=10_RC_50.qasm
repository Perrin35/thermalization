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
rz(-3.0769899) q[0];
sx q[0];
rz(-0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0639122) q[0];
sx q[0];
rz(-0.67718107) q[0];
sx q[0];
rz(2.1190686) q[0];
x q[1];
rz(-1.7982499) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(1.4129461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0845619) q[1];
sx q[1];
rz(-1.44094) q[1];
sx q[1];
rz(-0.76743857) q[1];
rz(-1.1159775) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(-1.3779373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3246831) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(-0.88513199) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-0.0016454776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5434108) q[0];
sx q[0];
rz(-0.046292154) q[0];
sx q[0];
rz(1.339519) q[0];
rz(-0.64237853) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(-3.1285398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8811223) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(1.2299728) q[1];
x q[2];
rz(-2.9716773) q[3];
sx q[3];
rz(-1.8222457) q[3];
sx q[3];
rz(-0.91472317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(-2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(0.99350199) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(1.089383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3141146) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(-0.88721888) q[0];
x q[1];
rz(-2.543534) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(-1.2146815) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8100064) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(1.7942795) q[1];
rz(-0.55604071) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(1.8204821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.115288) q[0];
sx q[0];
rz(-1.3867154) q[0];
sx q[0];
rz(0.91362761) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6845409) q[2];
sx q[2];
rz(-0.92278381) q[2];
sx q[2];
rz(1.5485473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.847536) q[1];
sx q[1];
rz(-1.1514074) q[1];
sx q[1];
rz(1.9860752) q[1];
rz(-pi) q[2];
rz(1.2218277) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(1.7581913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.416052) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(2.1062772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-2.007764) q[0];
sx q[0];
rz(2.2321738) q[0];
x q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(2.5626593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6742179) q[1];
sx q[1];
rz(-1.6661223) q[1];
sx q[1];
rz(-1.5159038) q[1];
rz(-pi) q[2];
rz(-1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(-1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6672872) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(-2.5179203) q[0];
rz(0.63164288) q[2];
sx q[2];
rz(-0.70438671) q[2];
sx q[2];
rz(2.7001675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44240272) q[1];
sx q[1];
rz(-1.3283722) q[1];
sx q[1];
rz(1.3139903) q[1];
x q[2];
rz(0.16222555) q[3];
sx q[3];
rz(-0.82419206) q[3];
sx q[3];
rz(1.1713067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(0.73202837) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7043982) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-2.722446) q[0];
rz(-1.8998434) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(2.3925051) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-2.8357361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52280207) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(-0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(0.83731246) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(-3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(2.1829139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(-2.3575248) q[0];
x q[1];
rz(-2.0197228) q[2];
sx q[2];
rz(-0.58510963) q[2];
sx q[2];
rz(-1.4657071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2221335) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(-1.7807351) q[1];
x q[2];
rz(-2.9577191) q[3];
sx q[3];
rz(-2.549578) q[3];
sx q[3];
rz(-1.5541935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(-1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(2.7046955) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(-0.46554309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(-2.4337016) q[0];
rz(-1.8579673) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(-0.62921333) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8852639) q[1];
sx q[1];
rz(-1.1210821) q[1];
sx q[1];
rz(-1.9535669) q[1];
x q[2];
rz(-0.97165473) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(-0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(0.75138584) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-0.5823935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56349194) q[0];
sx q[0];
rz(-2.2390215) q[0];
sx q[0];
rz(2.8816954) q[0];
x q[1];
rz(1.8495314) q[2];
sx q[2];
rz(-2.6346452) q[2];
sx q[2];
rz(2.7915733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.460234) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.7540936) q[1];
x q[2];
rz(-2.9243745) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(-0.56009968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-2.8600678) q[2];
sx q[2];
rz(-0.83419656) q[2];
sx q[2];
rz(0.4783334) q[2];
rz(2.3253757) q[3];
sx q[3];
rz(-2.4945138) q[3];
sx q[3];
rz(2.1993779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];