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
rz(1.3951294) q[0];
sx q[0];
rz(-2.3869618) q[0];
sx q[0];
rz(-1.0838497) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9636757) q[0];
sx q[0];
rz(-0.48030329) q[0];
sx q[0];
rz(2.1558574) q[0];
rz(2.0332912) q[2];
sx q[2];
rz(-1.9841737) q[2];
sx q[2];
rz(-0.55962901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1281761) q[1];
sx q[1];
rz(-2.8833296) q[1];
sx q[1];
rz(-1.5034666) q[1];
rz(-pi) q[2];
rz(-0.45248078) q[3];
sx q[3];
rz(-0.47305952) q[3];
sx q[3];
rz(-0.61247952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1276663) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(-2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056331228) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(-2.7714609) q[0];
rz(-0.47259694) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(-0.95742375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2070302) q[0];
sx q[0];
rz(-2.0683204) q[0];
sx q[0];
rz(-0.091188641) q[0];
rz(-pi) q[1];
rz(2.967011) q[2];
sx q[2];
rz(-2.0160297) q[2];
sx q[2];
rz(-3.0246459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70499252) q[1];
sx q[1];
rz(-0.7448405) q[1];
sx q[1];
rz(3.0158494) q[1];
rz(-pi) q[2];
rz(2.7368746) q[3];
sx q[3];
rz(-2.4232695) q[3];
sx q[3];
rz(-2.5248034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1364253) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(2.1180507) q[2];
rz(1.7510022) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(-1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(2.5787831) q[0];
rz(1.5142745) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-0.079158457) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25185967) q[0];
sx q[0];
rz(-1.2556228) q[0];
sx q[0];
rz(0.51879518) q[0];
rz(-pi) q[1];
rz(-0.49071006) q[2];
sx q[2];
rz(-1.6067925) q[2];
sx q[2];
rz(-1.7902059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1677163) q[1];
sx q[1];
rz(-1.9922246) q[1];
sx q[1];
rz(-0.80444471) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3615846) q[3];
sx q[3];
rz(-0.87015188) q[3];
sx q[3];
rz(-2.8740945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9074273) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(0.35587707) q[3];
sx q[3];
rz(-1.1938813) q[3];
sx q[3];
rz(-2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45977649) q[0];
sx q[0];
rz(-1.9720607) q[0];
sx q[0];
rz(1.1112777) q[0];
rz(0.92012826) q[1];
sx q[1];
rz(-1.5891113) q[1];
sx q[1];
rz(-0.2535325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7074575) q[0];
sx q[0];
rz(-0.61096707) q[0];
sx q[0];
rz(2.5649628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6854671) q[2];
sx q[2];
rz(-1.7301699) q[2];
sx q[2];
rz(-1.8346952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.580198) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(-2.8785588) q[1];
rz(-1.0658137) q[3];
sx q[3];
rz(-1.8441213) q[3];
sx q[3];
rz(-1.2697288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3229052) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(-1.4595855) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.9394268) q[3];
sx q[3];
rz(-0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(-1.7942418) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(-1.7151054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7280557) q[0];
sx q[0];
rz(-1.9621092) q[0];
sx q[0];
rz(-1.9584459) q[0];
x q[1];
rz(2.5002062) q[2];
sx q[2];
rz(-1.25433) q[2];
sx q[2];
rz(-0.14009012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0795143) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(0.18194992) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35512202) q[3];
sx q[3];
rz(-1.5105845) q[3];
sx q[3];
rz(-2.9400241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9023989) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(2.9854767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657144) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(-2.2789047) q[0];
rz(-2.8942096) q[1];
sx q[1];
rz(-1.8377973) q[1];
sx q[1];
rz(-2.6695796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.056045) q[0];
sx q[0];
rz(-1.2331095) q[0];
sx q[0];
rz(-2.4187536) q[0];
rz(-2.9905867) q[2];
sx q[2];
rz(-1.0698294) q[2];
sx q[2];
rz(0.46268625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1122163) q[1];
sx q[1];
rz(-2.7885155) q[1];
sx q[1];
rz(-0.53332163) q[1];
rz(-1.2788762) q[3];
sx q[3];
rz(-2.4987767) q[3];
sx q[3];
rz(-2.0187279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(-1.3989353) q[2];
rz(-1.6861964) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(2.5808891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48651925) q[0];
sx q[0];
rz(-1.9955248) q[0];
sx q[0];
rz(2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(-1.8124883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603921) q[0];
sx q[0];
rz(-2.5787163) q[0];
sx q[0];
rz(2.2928574) q[0];
rz(1.7460004) q[2];
sx q[2];
rz(-2.8602798) q[2];
sx q[2];
rz(0.11872053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74187255) q[1];
sx q[1];
rz(-2.406027) q[1];
sx q[1];
rz(1.5217785) q[1];
x q[2];
rz(0.01984486) q[3];
sx q[3];
rz(-2.1571419) q[3];
sx q[3];
rz(-1.8137518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3507877) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(-0.038185509) q[2];
rz(-0.83032483) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(-3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4627948) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(-1.7713254) q[0];
rz(-2.7554152) q[1];
sx q[1];
rz(-1.736707) q[1];
sx q[1];
rz(-0.6699627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46028194) q[0];
sx q[0];
rz(-0.71397266) q[0];
sx q[0];
rz(-2.332469) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0490446) q[2];
sx q[2];
rz(-2.4227648) q[2];
sx q[2];
rz(-2.3859442) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8893056) q[1];
sx q[1];
rz(-0.60599594) q[1];
sx q[1];
rz(-1.7701946) q[1];
x q[2];
rz(1.3588328) q[3];
sx q[3];
rz(-1.5342714) q[3];
sx q[3];
rz(-1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.92114821) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(1.05668) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073762745) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(-2.1347866) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(-0.83008343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842889) q[0];
sx q[0];
rz(-2.8606133) q[0];
sx q[0];
rz(0.12019867) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3557498) q[2];
sx q[2];
rz(-0.98282645) q[2];
sx q[2];
rz(2.7933592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0574903) q[1];
sx q[1];
rz(-1.7508645) q[1];
sx q[1];
rz(1.4157285) q[1];
rz(-pi) q[2];
rz(-3.1287635) q[3];
sx q[3];
rz(-2.0137824) q[3];
sx q[3];
rz(-1.6010546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-1.0969561) q[2];
rz(1.8638301) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(-0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(0.37503234) q[0];
rz(3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(-0.92432252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40709201) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(1.2248216) q[0];
rz(-pi) q[1];
rz(1.0607988) q[2];
sx q[2];
rz(-0.96335232) q[2];
sx q[2];
rz(0.80112544) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.754318) q[1];
sx q[1];
rz(-1.5421093) q[1];
sx q[1];
rz(-2.7806675) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64226867) q[3];
sx q[3];
rz(-2.6005201) q[3];
sx q[3];
rz(-0.16995811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(-2.2484153) q[2];
rz(1.1183974) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669582) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(-0.18192667) q[1];
sx q[1];
rz(-2.5839099) q[1];
sx q[1];
rz(-0.90669496) q[1];
rz(1.2557658) q[2];
sx q[2];
rz(-1.6581104) q[2];
sx q[2];
rz(2.8078512) q[2];
rz(0.652486) q[3];
sx q[3];
rz(-2.5204044) q[3];
sx q[3];
rz(1.6922097) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
