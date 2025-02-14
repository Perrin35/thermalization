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
rz(2.9003484) q[0];
sx q[0];
rz(-2.9304507) q[0];
sx q[0];
rz(-0.40786064) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(-0.064755138) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57015837) q[0];
sx q[0];
rz(-1.4553833) q[0];
sx q[0];
rz(-3.136271) q[0];
rz(-pi) q[1];
rz(0.98528905) q[2];
sx q[2];
rz(-1.3006214) q[2];
sx q[2];
rz(-2.5463242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9038856) q[1];
sx q[1];
rz(-0.14169417) q[1];
sx q[1];
rz(3.0631439) q[1];
rz(-pi) q[2];
rz(1.466732) q[3];
sx q[3];
rz(-0.95115653) q[3];
sx q[3];
rz(2.5980169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9593418) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(1.135929) q[2];
rz(0.7302537) q[3];
sx q[3];
rz(-1.3948995) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.092904329) q[0];
sx q[0];
rz(-2.1850259) q[0];
sx q[0];
rz(-3.0862869) q[0];
rz(-2.7703908) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(-0.41195437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5434951) q[0];
sx q[0];
rz(-2.6109222) q[0];
sx q[0];
rz(3.0125488) q[0];
rz(-pi) q[1];
rz(0.29779224) q[2];
sx q[2];
rz(-1.8911368) q[2];
sx q[2];
rz(-3.0646119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69393278) q[1];
sx q[1];
rz(-2.3946163) q[1];
sx q[1];
rz(-0.20836094) q[1];
x q[2];
rz(-0.92276908) q[3];
sx q[3];
rz(-1.8776181) q[3];
sx q[3];
rz(2.0056779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61341316) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(-0.8826274) q[2];
rz(1.708185) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9751137) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(-0.56667462) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.5387225) q[1];
sx q[1];
rz(1.3272939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9629563) q[0];
sx q[0];
rz(-1.4861075) q[0];
sx q[0];
rz(0.87015193) q[0];
rz(-1.7812875) q[2];
sx q[2];
rz(-2.5062525) q[2];
sx q[2];
rz(-0.47266211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3396521) q[1];
sx q[1];
rz(-1.2627456) q[1];
sx q[1];
rz(3.1402474) q[1];
rz(-2.0762844) q[3];
sx q[3];
rz(-0.79366854) q[3];
sx q[3];
rz(2.6738559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7495482) q[2];
sx q[2];
rz(-2.0508524) q[2];
sx q[2];
rz(0.90847477) q[2];
rz(-2.3467017) q[3];
sx q[3];
rz(-2.0386233) q[3];
sx q[3];
rz(-0.046549646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6598776) q[0];
sx q[0];
rz(-0.91232038) q[0];
sx q[0];
rz(0.6657486) q[0];
rz(0.84716973) q[1];
sx q[1];
rz(-0.24996346) q[1];
sx q[1];
rz(0.4096823) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9315791) q[0];
sx q[0];
rz(-1.0780452) q[0];
sx q[0];
rz(-1.3973622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3329725) q[2];
sx q[2];
rz(-1.8372155) q[2];
sx q[2];
rz(3.1297822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.078881513) q[1];
sx q[1];
rz(-2.1797175) q[1];
sx q[1];
rz(-0.2746065) q[1];
x q[2];
rz(1.631103) q[3];
sx q[3];
rz(-0.44476393) q[3];
sx q[3];
rz(-2.841265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4117671) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(1.3789619) q[2];
rz(-1.1369368) q[3];
sx q[3];
rz(-1.3552908) q[3];
sx q[3];
rz(-0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1788504) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(0.41819292) q[0];
rz(1.7182982) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(-1.6421912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2432775) q[0];
sx q[0];
rz(-1.1278129) q[0];
sx q[0];
rz(-1.8421296) q[0];
rz(-1.8756494) q[2];
sx q[2];
rz(-1.8902167) q[2];
sx q[2];
rz(1.6736506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5852927) q[1];
sx q[1];
rz(-1.5130116) q[1];
sx q[1];
rz(1.3520266) q[1];
rz(0.065369936) q[3];
sx q[3];
rz(-1.4413068) q[3];
sx q[3];
rz(-1.1060028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8666009) q[2];
sx q[2];
rz(-1.1555221) q[2];
sx q[2];
rz(0.17821136) q[2];
rz(-2.7675659) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(-1.1427574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571534) q[0];
sx q[0];
rz(-1.3244119) q[0];
sx q[0];
rz(-1.4763747) q[0];
rz(-0.58552512) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(2.4077328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39291362) q[0];
sx q[0];
rz(-2.3322912) q[0];
sx q[0];
rz(2.5038263) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30690212) q[2];
sx q[2];
rz(-1.8820805) q[2];
sx q[2];
rz(2.5800623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3702995) q[1];
sx q[1];
rz(-1.7669189) q[1];
sx q[1];
rz(-2.0972154) q[1];
rz(-2.9328314) q[3];
sx q[3];
rz(-2.0144297) q[3];
sx q[3];
rz(2.7334653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.057498589) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(-1.499184) q[2];
rz(1.5205787) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(2.9852941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3868189) q[0];
sx q[0];
rz(-2.2235625) q[0];
sx q[0];
rz(-1.7479489) q[0];
rz(0.58760324) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(-2.844152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5597683) q[0];
sx q[0];
rz(-0.80162573) q[0];
sx q[0];
rz(2.0770492) q[0];
rz(-pi) q[1];
rz(-1.3002505) q[2];
sx q[2];
rz(-1.5890117) q[2];
sx q[2];
rz(-1.2684938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6828281) q[1];
sx q[1];
rz(-1.5804844) q[1];
sx q[1];
rz(2.2504456) q[1];
x q[2];
rz(-2.2138811) q[3];
sx q[3];
rz(-1.7684019) q[3];
sx q[3];
rz(1.2090982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8224767) q[2];
sx q[2];
rz(-1.7812704) q[2];
sx q[2];
rz(2.759867) q[2];
rz(-0.60025674) q[3];
sx q[3];
rz(-0.67265066) q[3];
sx q[3];
rz(-2.8762347) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(-1.8040682) q[0];
rz(-0.0094825347) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(1.7705852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776896) q[0];
sx q[0];
rz(-0.36720095) q[0];
sx q[0];
rz(-0.3911256) q[0];
rz(2.2393537) q[2];
sx q[2];
rz(-1.1663365) q[2];
sx q[2];
rz(2.8079833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43804719) q[1];
sx q[1];
rz(-0.79709541) q[1];
sx q[1];
rz(-0.17532562) q[1];
x q[2];
rz(-2.9389006) q[3];
sx q[3];
rz(-0.42169844) q[3];
sx q[3];
rz(1.7360753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81801549) q[2];
sx q[2];
rz(-2.24021) q[2];
sx q[2];
rz(-0.14249194) q[2];
rz(0.88322181) q[3];
sx q[3];
rz(-1.764069) q[3];
sx q[3];
rz(1.6787329) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0505117) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(1.2485414) q[0];
rz(0.40731353) q[1];
sx q[1];
rz(-1.7946323) q[1];
sx q[1];
rz(-1.9479082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7418967) q[0];
sx q[0];
rz(-0.4238766) q[0];
sx q[0];
rz(2.2370643) q[0];
rz(-pi) q[1];
rz(-0.87393729) q[2];
sx q[2];
rz(-0.94183445) q[2];
sx q[2];
rz(-0.12128092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5496918) q[1];
sx q[1];
rz(-0.61335269) q[1];
sx q[1];
rz(-2.5332622) q[1];
rz(-pi) q[2];
x q[2];
rz(1.570669) q[3];
sx q[3];
rz(-2.8817892) q[3];
sx q[3];
rz(-2.3854431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(-2.0286782) q[2];
rz(3.0985966) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-2.903741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064903108) q[0];
sx q[0];
rz(-1.4142798) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(-1.9211357) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(0.97741309) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61924926) q[0];
sx q[0];
rz(-0.27740208) q[0];
sx q[0];
rz(-2.7707556) q[0];
x q[1];
rz(-1.270411) q[2];
sx q[2];
rz(-1.4255878) q[2];
sx q[2];
rz(1.5086439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4949242) q[1];
sx q[1];
rz(-0.99993247) q[1];
sx q[1];
rz(-0.60653703) q[1];
rz(-pi) q[2];
rz(0.39267003) q[3];
sx q[3];
rz(-1.304198) q[3];
sx q[3];
rz(-2.6368898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0105373) q[2];
sx q[2];
rz(-2.8752893) q[2];
sx q[2];
rz(-1.4492501) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-2.0802458) q[3];
sx q[3];
rz(1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5129678) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.3799008) q[1];
sx q[1];
rz(-1.6086144) q[1];
sx q[1];
rz(-0.10133941) q[1];
rz(-2.9607282) q[2];
sx q[2];
rz(-1.7500063) q[2];
sx q[2];
rz(1.0371006) q[2];
rz(0.64639277) q[3];
sx q[3];
rz(-1.7714995) q[3];
sx q[3];
rz(-2.5378791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
