OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(5.2390715) q[1];
sx q[1];
rz(10.704438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5582433) q[0];
sx q[0];
rz(-1.8177176) q[0];
sx q[0];
rz(2.9980744) q[0];
rz(1.7945292) q[2];
sx q[2];
rz(-1.9828084) q[2];
sx q[2];
rz(1.7681233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0919839) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(-2.0642573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23366191) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(-1.83028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90855234) q[0];
sx q[0];
rz(-1.6516343) q[0];
sx q[0];
rz(0.064706133) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1335282) q[2];
sx q[2];
rz(-1.9568866) q[2];
sx q[2];
rz(-1.776945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0812379) q[1];
sx q[1];
rz(-1.758467) q[1];
sx q[1];
rz(-2.9802122) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4653141) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(-2.7200384) q[0];
x q[1];
rz(-2.9163755) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-2.1622216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.682444) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(0.1174121) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95007105) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(-2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-2.2213675) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5266787) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(-2.148669) q[0];
rz(-pi) q[1];
rz(-1.6572861) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(-2.6038225) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8060382) q[1];
sx q[1];
rz(-2.028855) q[1];
sx q[1];
rz(-0.6946509) q[1];
rz(-0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(-2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.2801923) q[2];
rz(0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950726) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(-1.1016269) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3327417) q[2];
sx q[2];
rz(-0.14698725) q[2];
sx q[2];
rz(0.57839314) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0799775) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(1.3498989) q[1];
rz(-2.1987678) q[3];
sx q[3];
rz(-0.70126611) q[3];
sx q[3];
rz(-3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.5135117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-0.61074722) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.020535) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(-0.72292098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7102393) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(0.61183521) q[1];
rz(-0.87168872) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-0.17091621) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8531187) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-2.8209646) q[0];
x q[1];
rz(-1.600012) q[2];
sx q[2];
rz(-1.6299106) q[2];
sx q[2];
rz(2.7616449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4759051) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(1.3728632) q[1];
rz(0.62992923) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(1.5671974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-0.9986977) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(2.1112679) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1925826) q[0];
sx q[0];
rz(-2.2826676) q[0];
sx q[0];
rz(-2.683995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0439992) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(-0.85925697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1134539) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(1.0949275) q[1];
rz(2.2780667) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(-1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410256) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(2.7450949) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1124434) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(-1.0708035) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2444297) q[1];
sx q[1];
rz(-2.0532554) q[1];
sx q[1];
rz(3.0753115) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9477356) q[3];
sx q[3];
rz(-1.590168) q[3];
sx q[3];
rz(-1.7957578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(-2.9454254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575532) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-0.89417017) q[0];
rz(-1.7061383) q[2];
sx q[2];
rz(-1.6207098) q[2];
sx q[2];
rz(-0.70336715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13388053) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(-1.5478565) q[1];
x q[2];
rz(0.65532834) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(0.67656803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3907884) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-1.0749764) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.3148057) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(-0.20028533) q[2];
sx q[2];
rz(-0.23858783) q[2];
sx q[2];
rz(-1.5993652) q[2];
rz(0.046394596) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];