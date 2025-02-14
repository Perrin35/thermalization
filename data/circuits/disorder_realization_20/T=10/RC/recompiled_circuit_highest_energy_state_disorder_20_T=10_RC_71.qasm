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
rz(-2.1751997) q[0];
sx q[0];
rz(4.8973358) q[0];
sx q[0];
rz(10.020221) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(-0.59102494) q[1];
sx q[1];
rz(-3.0182867) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6791413) q[0];
sx q[0];
rz(-1.9872905) q[0];
sx q[0];
rz(0.97521675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79084556) q[2];
sx q[2];
rz(-2.186785) q[2];
sx q[2];
rz(-0.14033422) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5414303) q[1];
sx q[1];
rz(-1.3568391) q[1];
sx q[1];
rz(-3.0453277) q[1];
x q[2];
rz(2.7601065) q[3];
sx q[3];
rz(-2.7499928) q[3];
sx q[3];
rz(2.7158383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3143602) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(-3.0296791) q[2];
rz(-1.3917475) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157442) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-2.9916905) q[0];
rz(0.28597486) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(2.012595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99230951) q[0];
sx q[0];
rz(-0.2899) q[0];
sx q[0];
rz(1.4647746) q[0];
rz(-pi) q[1];
rz(1.8104042) q[2];
sx q[2];
rz(-3.0802266) q[2];
sx q[2];
rz(-0.97739109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6817878) q[1];
sx q[1];
rz(-0.39783172) q[1];
sx q[1];
rz(2.7017606) q[1];
x q[2];
rz(-0.82650696) q[3];
sx q[3];
rz(-0.77946589) q[3];
sx q[3];
rz(-2.731088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67843208) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(-1.9909667) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(2.7599957) q[0];
rz(1.2941788) q[1];
sx q[1];
rz(-2.0353863) q[1];
sx q[1];
rz(2.9433184) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564045) q[0];
sx q[0];
rz(-1.745178) q[0];
sx q[0];
rz(1.8013493) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4691822) q[2];
sx q[2];
rz(-0.67199113) q[2];
sx q[2];
rz(0.72173126) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1041996) q[1];
sx q[1];
rz(-1.3379119) q[1];
sx q[1];
rz(1.7892455) q[1];
x q[2];
rz(0.020129344) q[3];
sx q[3];
rz(-1.3960037) q[3];
sx q[3];
rz(2.54769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76564378) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(0.51885968) q[2];
rz(1.2505924) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64312235) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(-0.44922391) q[0];
rz(-1.6953702) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(2.5208688) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45263824) q[0];
sx q[0];
rz(-1.4264548) q[0];
sx q[0];
rz(2.4792433) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6942107) q[2];
sx q[2];
rz(-2.2006249) q[2];
sx q[2];
rz(0.067867756) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6606969) q[1];
sx q[1];
rz(-0.45240739) q[1];
sx q[1];
rz(0.47641944) q[1];
x q[2];
rz(-3.0884864) q[3];
sx q[3];
rz(-2.1087007) q[3];
sx q[3];
rz(-2.2369179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(-0.16109666) q[2];
rz(0.39786878) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(-2.9919992) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75905269) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(-2.8845442) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.2535198) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2855447) q[0];
sx q[0];
rz(-3.1012707) q[0];
sx q[0];
rz(-2.0364385) q[0];
x q[1];
rz(-2.9078616) q[2];
sx q[2];
rz(-2.3178551) q[2];
sx q[2];
rz(-1.9875789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4684126) q[1];
sx q[1];
rz(-0.77408964) q[1];
sx q[1];
rz(1.5062529) q[1];
rz(3.1367132) q[3];
sx q[3];
rz(-2.7944428) q[3];
sx q[3];
rz(0.43970575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72838655) q[2];
sx q[2];
rz(-1.9325958) q[2];
sx q[2];
rz(1.8750635) q[2];
rz(-2.5201216) q[3];
sx q[3];
rz(-2.5340243) q[3];
sx q[3];
rz(2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(3.1165282) q[0];
rz(1.2114245) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-0.2549583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716041) q[0];
sx q[0];
rz(-1.5758262) q[0];
sx q[0];
rz(0.049335376) q[0];
rz(-pi) q[1];
rz(2.5549501) q[2];
sx q[2];
rz(-1.987118) q[2];
sx q[2];
rz(-3.0468009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7656169) q[1];
sx q[1];
rz(-1.0436397) q[1];
sx q[1];
rz(-0.65816452) q[1];
rz(-2.1433349) q[3];
sx q[3];
rz(-1.7504331) q[3];
sx q[3];
rz(-1.6999753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24594626) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(2.6825405) q[2];
rz(-1.6445271) q[3];
sx q[3];
rz(-3.0247757) q[3];
sx q[3];
rz(-3.0204401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(-0.91880265) q[1];
sx q[1];
rz(-2.0252392) q[1];
sx q[1];
rz(1.2109717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.092441) q[0];
sx q[0];
rz(-1.7730129) q[0];
sx q[0];
rz(-0.33009712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2077256) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(-3.0702555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1001491) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(-0.99975296) q[1];
rz(-pi) q[2];
rz(1.2065229) q[3];
sx q[3];
rz(-1.2281872) q[3];
sx q[3];
rz(-0.87876696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81898895) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(0.71713478) q[2];
rz(-1.0273733) q[3];
sx q[3];
rz(-1.4077978) q[3];
sx q[3];
rz(-2.7023442) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913919) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(0.014658654) q[0];
rz(-0.75153366) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(-1.470648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178407) q[0];
sx q[0];
rz(-1.7589594) q[0];
sx q[0];
rz(-1.4298196) q[0];
rz(-pi) q[1];
rz(-0.29701091) q[2];
sx q[2];
rz(-2.8737465) q[2];
sx q[2];
rz(1.0849407) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4336509) q[1];
sx q[1];
rz(-1.6835064) q[1];
sx q[1];
rz(-1.9757084) q[1];
rz(-3.1378391) q[3];
sx q[3];
rz(-1.225718) q[3];
sx q[3];
rz(1.0964637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8810001) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(-1.1547487) q[3];
sx q[3];
rz(-0.65360779) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(1.0754841) q[0];
rz(3.0134046) q[1];
sx q[1];
rz(-2.2905541) q[1];
sx q[1];
rz(0.34559616) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3227279) q[0];
sx q[0];
rz(-0.68536192) q[0];
sx q[0];
rz(-1.3166053) q[0];
rz(2.1394452) q[2];
sx q[2];
rz(-1.8446184) q[2];
sx q[2];
rz(-1.3783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3010878) q[1];
sx q[1];
rz(-1.6337758) q[1];
sx q[1];
rz(-2.5999864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7365428) q[3];
sx q[3];
rz(-2.1624544) q[3];
sx q[3];
rz(1.4480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(1.6205988) q[2];
rz(1.3711551) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(-2.7267406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82776752) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-0.23571043) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.3689573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018250759) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-0.98066179) q[0];
rz(-2.4659285) q[2];
sx q[2];
rz(-1.079529) q[2];
sx q[2];
rz(1.5671687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9006834) q[1];
sx q[1];
rz(-0.59421173) q[1];
sx q[1];
rz(-1.2628984) q[1];
x q[2];
rz(0.36054109) q[3];
sx q[3];
rz(-1.4829794) q[3];
sx q[3];
rz(2.5537864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(-3.0885922) q[2];
rz(-2.3897589) q[3];
sx q[3];
rz(-2.1078608) q[3];
sx q[3];
rz(-0.92541614) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5175405) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(1.3507631) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(2.4443632) q[2];
sx q[2];
rz(-1.4705428) q[2];
sx q[2];
rz(-1.85208) q[2];
rz(2.3334669) q[3];
sx q[3];
rz(-1.1302409) q[3];
sx q[3];
rz(-1.7779868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
