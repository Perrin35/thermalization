OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4176183) q[0];
sx q[0];
rz(-1.4899878) q[0];
sx q[0];
rz(2.2111501) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(-2.0342483) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005172) q[0];
sx q[0];
rz(-2.374875) q[0];
sx q[0];
rz(2.6630152) q[0];
rz(2.1030175) q[2];
sx q[2];
rz(-2.3468809) q[2];
sx q[2];
rz(0.46193916) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90802898) q[1];
sx q[1];
rz(-1.6828013) q[1];
sx q[1];
rz(1.2909375) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34875617) q[3];
sx q[3];
rz(-0.96491279) q[3];
sx q[3];
rz(-0.44266686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-2.1556373) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-2.3243288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.9190448) q[0];
rz(-pi) q[1];
rz(0.21821071) q[2];
sx q[2];
rz(-1.0901703) q[2];
sx q[2];
rz(2.6618119) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4546928) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(-1.8886186) q[1];
rz(-pi) q[2];
rz(2.9528584) q[3];
sx q[3];
rz(-1.1047603) q[3];
sx q[3];
rz(2.1621494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(0.36402738) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36689511) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(0.96631518) q[0];
rz(-0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.4470709) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9896302) q[0];
sx q[0];
rz(-1.785779) q[0];
sx q[0];
rz(1.5528029) q[0];
rz(-pi) q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(1.9384055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5406394) q[1];
sx q[1];
rz(-2.1689231) q[1];
sx q[1];
rz(-1.1126493) q[1];
rz(-1.0560016) q[3];
sx q[3];
rz(-2.6714973) q[3];
sx q[3];
rz(-1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7029999) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(-1.1052216) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-0.38898653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3272414) q[0];
sx q[0];
rz(-1.8530493) q[0];
sx q[0];
rz(1.3319356) q[0];
rz(1.9289891) q[2];
sx q[2];
rz(-0.51414031) q[2];
sx q[2];
rz(-2.6280623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7547238) q[1];
sx q[1];
rz(-1.1479953) q[1];
sx q[1];
rz(-0.60476426) q[1];
rz(1.7763406) q[3];
sx q[3];
rz(-2.3081429) q[3];
sx q[3];
rz(-1.8635686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7230364) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(2.8592498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-0.49450758) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617103) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(0.42219992) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.102166) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(-0.53265041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8066986) q[1];
sx q[1];
rz(-2.1999199) q[1];
sx q[1];
rz(-1.4966399) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4440998) q[3];
sx q[3];
rz(-1.5478304) q[3];
sx q[3];
rz(2.3820153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(-0.29850706) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.235199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081731) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(-0.74605201) q[0];
rz(-pi) q[1];
rz(1.9094798) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(-1.0519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0732121) q[1];
sx q[1];
rz(-0.74406032) q[1];
sx q[1];
rz(-2.5030604) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6545717) q[3];
sx q[3];
rz(-1.8063074) q[3];
sx q[3];
rz(-1.868353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(0.43169272) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(-0.11211638) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33681413) q[0];
sx q[0];
rz(-1.3114197) q[0];
sx q[0];
rz(-0.39774261) q[0];
rz(-2.8261975) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(2.18404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0472764) q[1];
sx q[1];
rz(-2.2502406) q[1];
sx q[1];
rz(1.9654771) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1464305) q[3];
sx q[3];
rz(-2.6570435) q[3];
sx q[3];
rz(0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-3.0155638) q[2];
rz(2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.9326899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-0.27420843) q[0];
sx q[0];
rz(-1.8673531) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65931321) q[2];
sx q[2];
rz(-0.94617832) q[2];
sx q[2];
rz(0.63559947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3371256) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(-2.0130403) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8672089) q[3];
sx q[3];
rz(-2.489438) q[3];
sx q[3];
rz(-1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.3195999) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(0.035141703) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(-1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-2.7499054) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93145934) q[0];
sx q[0];
rz(-2.0053232) q[0];
sx q[0];
rz(1.5772485) q[0];
rz(2.7462256) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(1.5026827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(0.4391567) q[1];
rz(0.12556062) q[3];
sx q[3];
rz(-0.52913044) q[3];
sx q[3];
rz(2.373113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(0.76198602) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(-1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-2.9796519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83308342) q[0];
sx q[0];
rz(-2.1662795) q[0];
sx q[0];
rz(-2.7916629) q[0];
rz(-pi) q[1];
rz(-1.6185206) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(-1.7392841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3871799) q[1];
sx q[1];
rz(-1.3993235) q[1];
sx q[1];
rz(-3.0110554) q[1];
rz(-3.083769) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(3.0081089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1841715) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.4282248) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(-1.7715001) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(1.7182405) q[2];
sx q[2];
rz(-1.4530008) q[2];
sx q[2];
rz(3.0086649) q[2];
rz(-2.6876642) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
