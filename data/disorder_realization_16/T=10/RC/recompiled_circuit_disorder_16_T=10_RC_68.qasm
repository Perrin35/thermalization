OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96903893) q[0];
sx q[0];
rz(-1.2533029) q[0];
sx q[0];
rz(2.88455) q[0];
x q[1];
rz(2.7726735) q[2];
sx q[2];
rz(-2.2152165) q[2];
sx q[2];
rz(1.3117787) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6602064) q[1];
sx q[1];
rz(-1.2848789) q[1];
sx q[1];
rz(1.5155161) q[1];
x q[2];
rz(-1.5987414) q[3];
sx q[3];
rz(-2.6290647) q[3];
sx q[3];
rz(-0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(-0.075803444) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(0.064963438) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(-1.2423135) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-2.8613052) q[0];
sx q[0];
rz(-2.6602402) q[0];
x q[1];
rz(1.785948) q[2];
sx q[2];
rz(-2.4977376) q[2];
sx q[2];
rz(-1.7807963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8426659) q[1];
sx q[1];
rz(-2.0999523) q[1];
sx q[1];
rz(-1.3859205) q[1];
rz(-0.55450704) q[3];
sx q[3];
rz(-0.79075659) q[3];
sx q[3];
rz(-2.8046372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80766455) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(-3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(2.1247991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.5350071) q[0];
sx q[0];
rz(-3.0599942) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52428464) q[2];
sx q[2];
rz(-2.0031843) q[2];
sx q[2];
rz(-1.0334894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9079202) q[1];
sx q[1];
rz(-2.3111812) q[1];
sx q[1];
rz(2.9552712) q[1];
rz(-pi) q[2];
rz(-1.73909) q[3];
sx q[3];
rz(-0.28545359) q[3];
sx q[3];
rz(-1.1807549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(2.8947815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.719559) q[0];
sx q[0];
rz(2.2674198) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7505789) q[2];
sx q[2];
rz(-0.51119971) q[2];
sx q[2];
rz(2.9875987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3226763) q[1];
sx q[1];
rz(-0.37441844) q[1];
sx q[1];
rz(-0.14426343) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.85535) q[3];
sx q[3];
rz(-0.97385588) q[3];
sx q[3];
rz(-2.3846574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.56303) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(0.24838233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800956) q[0];
sx q[0];
rz(-1.0307612) q[0];
sx q[0];
rz(1.6105152) q[0];
x q[1];
rz(2.9551198) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(-1.1764256) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(-2.8413248) q[1];
rz(-pi) q[2];
rz(1.8099269) q[3];
sx q[3];
rz(-2.2056747) q[3];
sx q[3];
rz(-0.90415108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1335063) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.7957934) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4437618) q[0];
sx q[0];
rz(-2.9450581) q[0];
sx q[0];
rz(2.6854808) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7331946) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(-2.1085395) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9166959) q[1];
sx q[1];
rz(-1.6273013) q[1];
sx q[1];
rz(2.4042261) q[1];
x q[2];
rz(0.23645225) q[3];
sx q[3];
rz(-1.7835622) q[3];
sx q[3];
rz(0.80602431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(2.52264) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(2.0708864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8585513) q[0];
sx q[0];
rz(-2.0548477) q[0];
sx q[0];
rz(-1.5714684) q[0];
x q[1];
rz(-2.7878739) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(-3.0596717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13094014) q[1];
sx q[1];
rz(-2.3619235) q[1];
sx q[1];
rz(-0.14203771) q[1];
rz(-pi) q[2];
rz(2.8093852) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(-1.2063684) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(2.1941197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652126) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(-3.0239848) q[0];
rz(-pi) q[1];
rz(0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(1.8159602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84711134) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(2.860445) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9318337) q[3];
sx q[3];
rz(-2.3449538) q[3];
sx q[3];
rz(2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(2.8619213) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(2.9387617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664353) q[0];
sx q[0];
rz(-0.046910722) q[0];
sx q[0];
rz(-2.5527918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3896396) q[2];
sx q[2];
rz(-0.30138902) q[2];
sx q[2];
rz(-1.4434659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1111787) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(1.4700252) q[1];
x q[2];
rz(-2.8392302) q[3];
sx q[3];
rz(-0.48067579) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-0.28082401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.330864) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(2.3848563) q[0];
rz(-1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(0.10373058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22051375) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(1.6284579) q[1];
x q[2];
rz(1.0481846) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(2.2446333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(-2.1394219) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(0.87396809) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(0.80550823) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(-2.279083) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
