OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74237139) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-0.39682927) q[0];
x q[1];
rz(1.7982499) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(1.7286466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6380438) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(1.391295) q[1];
rz(-0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(2.0781793) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(3.1399472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7749274) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(-3.1309743) q[0];
rz(-pi) q[1];
rz(1.8284945) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(1.7110273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.590608) q[1];
sx q[1];
rz(-1.3270757) q[1];
sx q[1];
rz(-2.3477712) q[1];
rz(2.9716773) q[3];
sx q[3];
rz(-1.8222457) q[3];
sx q[3];
rz(-2.2268695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-0.70811159) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8436962) q[0];
sx q[0];
rz(-1.9934405) q[0];
sx q[0];
rz(0.98590132) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.543534) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(-1.2146815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5685801) q[1];
sx q[1];
rz(-1.1862566) q[1];
sx q[1];
rz(0.09210715) q[1];
rz(0.55604071) q[3];
sx q[3];
rz(-2.2696113) q[3];
sx q[3];
rz(-1.3211105) q[3];
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
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7774178) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(1.2749519) q[0];
rz(-pi) q[1];
rz(-1.6845409) q[2];
sx q[2];
rz(-2.2188088) q[2];
sx q[2];
rz(-1.5930454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4543912) q[1];
sx q[1];
rz(-1.193421) q[1];
sx q[1];
rz(-2.6881933) q[1];
rz(-pi) q[2];
rz(-0.30712819) q[3];
sx q[3];
rz(-1.9046475) q[3];
sx q[3];
rz(-0.083837282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.416052) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(2.1062772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0400378) q[0];
sx q[0];
rz(-2.3674175) q[0];
sx q[0];
rz(2.220962) q[0];
x q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(-0.57893334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(-0.095468949) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4953793) q[3];
sx q[3];
rz(-2.8013902) q[3];
sx q[3];
rz(-0.37341213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(0.76888293) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(0.2072269) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6672872) q[0];
sx q[0];
rz(-1.4315839) q[0];
sx q[0];
rz(-2.5179203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63164288) q[2];
sx q[2];
rz(-2.4372059) q[2];
sx q[2];
rz(-0.44142516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8689649) q[1];
sx q[1];
rz(-0.3513063) q[1];
sx q[1];
rz(-2.3428194) q[1];
rz(-0.81760041) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1563819) q[0];
sx q[0];
rz(-1.9893601) q[0];
sx q[0];
rz(1.5147989) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0351699) q[2];
sx q[2];
rz(-1.2692045) q[2];
sx q[2];
rz(2.0472722) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2597255) q[1];
sx q[1];
rz(-1.423466) q[1];
sx q[1];
rz(-2.8357361) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2613867) q[3];
sx q[3];
rz(-1.0687807) q[3];
sx q[3];
rz(1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(-0.95867872) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6129235) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(-0.59707609) q[0];
x q[1];
rz(-0.27997048) q[2];
sx q[2];
rz(-1.0500056) q[2];
sx q[2];
rz(-2.1998646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6764076) q[1];
sx q[1];
rz(-1.0470111) q[1];
sx q[1];
rz(0.1233867) q[1];
rz(1.6931375) q[3];
sx q[3];
rz(-0.99007505) q[3];
sx q[3];
rz(1.8078705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(-1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-0.46554309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8252194) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(-2.4337016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8579673) q[2];
sx q[2];
rz(-0.42771491) q[2];
sx q[2];
rz(0.62921333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0034345) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(2.483063) q[1];
rz(-pi) q[2];
rz(-2.1699379) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(-3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(-0.75138584) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(0.5823935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1727027) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(-1.2560647) q[0];
rz(-pi) q[1];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.7047802) q[2];
sx q[2];
rz(-1.6756563) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.460234) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.7540936) q[1];
rz(-pi) q[2];
rz(0.12606975) q[3];
sx q[3];
rz(-2.9226916) q[3];
sx q[3];
rz(-2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(2.3790512) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(1.273524) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(2.0740261) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];