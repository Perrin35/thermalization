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
rz(1.8099161) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74237139) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-0.39682927) q[0];
rz(-pi) q[1];
rz(1.1515491) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-0.56564769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.761258) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(0.18591979) q[1];
rz(1.1159775) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(1.3779373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(3.1399472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20365276) q[0];
sx q[0];
rz(-1.5601888) q[0];
sx q[0];
rz(1.525735) q[0];
rz(-pi) q[1];
rz(0.33866377) q[2];
sx q[2];
rz(-0.67064697) q[2];
sx q[2];
rz(-1.2884969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.590608) q[1];
sx q[1];
rz(-1.814517) q[1];
sx q[1];
rz(0.79382146) q[1];
x q[2];
rz(-0.1699154) q[3];
sx q[3];
rz(-1.8222457) q[3];
sx q[3];
rz(0.91472317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(2.0522096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2978965) q[0];
sx q[0];
rz(-1.1481522) q[0];
sx q[0];
rz(-0.98590132) q[0];
rz(-pi) q[1];
rz(2.5554352) q[2];
sx q[2];
rz(-0.68550368) q[2];
sx q[2];
rz(-0.11867487) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.032420302) q[1];
sx q[1];
rz(-1.4854327) q[1];
sx q[1];
rz(-1.1847772) q[1];
rz(-2.3508164) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(-0.63000597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.737405) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(2.9106211) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6511351) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(-0.091094253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6105146) q[1];
sx q[1];
rz(-2.5602166) q[1];
sx q[1];
rz(0.73552144) q[1];
x q[2];
rz(-2.8344645) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(-0.083837282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.416052) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-1.1338286) q[0];
sx q[0];
rz(-2.2321738) q[0];
x q[1];
rz(-3.0683124) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(-0.57893334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94381911) q[1];
sx q[1];
rz(-3.0316331) q[1];
sx q[1];
rz(2.6206559) q[1];
rz(2.4953793) q[3];
sx q[3];
rz(-2.8013902) q[3];
sx q[3];
rz(0.37341213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(0.76888293) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(-1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795905) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80400318) q[0];
sx q[0];
rz(-0.9540671) q[0];
sx q[0];
rz(-1.7417275) q[0];
rz(-pi) q[1];
rz(2.5099498) q[2];
sx q[2];
rz(-0.70438671) q[2];
sx q[2];
rz(-2.7001675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8689649) q[1];
sx q[1];
rz(-0.3513063) q[1];
sx q[1];
rz(-2.3428194) q[1];
x q[2];
rz(-2.3239922) q[3];
sx q[3];
rz(-1.6896276) q[3];
sx q[3];
rz(-2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8032288) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-2.8815564) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
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
rz(-1.0193298) q[0];
sx q[0];
rz(-0.42207345) q[0];
sx q[0];
rz(0.12515573) q[0];
x q[1];
rz(3.0351699) q[2];
sx q[2];
rz(-1.8723882) q[2];
sx q[2];
rz(1.0943204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26475036) q[1];
sx q[1];
rz(-1.8732338) q[1];
sx q[1];
rz(1.4164063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52280207) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(3.0772046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(-2.1829139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129235) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(2.5445166) q[0];
x q[1];
rz(-2.8616222) q[2];
sx q[2];
rz(-1.0500056) q[2];
sx q[2];
rz(2.1998646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9194591) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(-1.3608576) q[1];
rz(-2.557425) q[3];
sx q[3];
rz(-1.468588) q[3];
sx q[3];
rz(2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(0.46554309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(0.70789106) q[0];
rz(1.2836254) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(-0.62921333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0034345) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(0.65852965) q[1];
rz(-1.2113308) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(-2.204493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(-1.643606) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1727027) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(-1.885528) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2920612) q[2];
sx q[2];
rz(-0.50694743) q[2];
sx q[2];
rz(-0.35001937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.13233391) q[1];
sx q[1];
rz(-1.3887822) q[1];
sx q[1];
rz(-0.11972129) q[1];
x q[2];
rz(-0.12606975) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(0.88760469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.3394042) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(1.273524) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
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