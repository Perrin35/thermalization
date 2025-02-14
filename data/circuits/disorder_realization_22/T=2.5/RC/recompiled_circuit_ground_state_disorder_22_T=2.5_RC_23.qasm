OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.206267) q[0];
sx q[0];
rz(-1.6214108) q[0];
sx q[0];
rz(-1.7229463) q[0];
rz(-3.0980134) q[1];
sx q[1];
rz(-1.8585304) q[1];
sx q[1];
rz(-0.15951523) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1698604) q[0];
sx q[0];
rz(-1.5992478) q[0];
sx q[0];
rz(-1.4048496) q[0];
rz(-pi) q[1];
rz(-0.61181991) q[2];
sx q[2];
rz(-1.7828336) q[2];
sx q[2];
rz(3.0481047) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39594138) q[1];
sx q[1];
rz(-0.49761727) q[1];
sx q[1];
rz(1.2785589) q[1];
x q[2];
rz(0.36002145) q[3];
sx q[3];
rz(-2.2593479) q[3];
sx q[3];
rz(3.036397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45304766) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(3.0513406) q[2];
rz(0.074617535) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(2.727865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2523044) q[0];
sx q[0];
rz(-1.1587208) q[0];
sx q[0];
rz(-1.4537551) q[0];
rz(2.3989035) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(-0.76505605) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1677645) q[0];
sx q[0];
rz(-2.544738) q[0];
sx q[0];
rz(3.0227227) q[0];
x q[1];
rz(2.3397331) q[2];
sx q[2];
rz(-0.86881402) q[2];
sx q[2];
rz(-3.1131256) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7530427) q[1];
sx q[1];
rz(-1.8990128) q[1];
sx q[1];
rz(-0.80372827) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54799796) q[3];
sx q[3];
rz(-1.5685076) q[3];
sx q[3];
rz(1.8016004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9146156) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(-1.4545308) q[2];
rz(-2.4948273) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(-1.2796992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.0525368) q[0];
sx q[0];
rz(-1.6707957) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(1.2308944) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(2.0848354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0522359) q[0];
sx q[0];
rz(-1.0228906) q[0];
sx q[0];
rz(0.38378451) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0750404) q[2];
sx q[2];
rz(-2.1556052) q[2];
sx q[2];
rz(0.57384795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18348098) q[1];
sx q[1];
rz(-0.90881244) q[1];
sx q[1];
rz(2.8173911) q[1];
x q[2];
rz(2.908024) q[3];
sx q[3];
rz(-2.2456777) q[3];
sx q[3];
rz(-0.44665953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4061654) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(0.65262922) q[2];
rz(-1.2930219) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3339313) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(-1.3941143) q[0];
rz(-1.8380503) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(-2.0533144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1980285) q[0];
sx q[0];
rz(-1.2965512) q[0];
sx q[0];
rz(-2.345577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2919664) q[2];
sx q[2];
rz(-1.0996981) q[2];
sx q[2];
rz(1.8249408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7409121) q[1];
sx q[1];
rz(-2.0517382) q[1];
sx q[1];
rz(-0.6948284) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6234947) q[3];
sx q[3];
rz(-2.0695588) q[3];
sx q[3];
rz(-2.3757405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3950222) q[2];
sx q[2];
rz(-1.7638548) q[2];
sx q[2];
rz(-0.48040381) q[2];
rz(-2.8754442) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0323623) q[0];
sx q[0];
rz(-0.6210331) q[0];
sx q[0];
rz(-1.9367628) q[0];
rz(-0.042639848) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(-2.1991594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54357994) q[0];
sx q[0];
rz(-2.6930598) q[0];
sx q[0];
rz(1.8878114) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93330002) q[2];
sx q[2];
rz(-2.0752873) q[2];
sx q[2];
rz(-1.8613033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3417899) q[1];
sx q[1];
rz(-1.1936578) q[1];
sx q[1];
rz(-2.1076635) q[1];
rz(-pi) q[2];
rz(2.5025401) q[3];
sx q[3];
rz(-1.6585232) q[3];
sx q[3];
rz(1.7114342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11613906) q[2];
sx q[2];
rz(-0.68486539) q[2];
sx q[2];
rz(-1.9174513) q[2];
rz(-0.72743607) q[3];
sx q[3];
rz(-1.6614611) q[3];
sx q[3];
rz(3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392451) q[0];
sx q[0];
rz(-2.6363063) q[0];
sx q[0];
rz(-2.8832866) q[0];
rz(-0.91336617) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(-2.1736653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0570404) q[0];
sx q[0];
rz(-1.0121317) q[0];
sx q[0];
rz(-1.1112122) q[0];
rz(-2.2931678) q[2];
sx q[2];
rz(-1.2944739) q[2];
sx q[2];
rz(1.487731) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0022626) q[1];
sx q[1];
rz(-1.5353629) q[1];
sx q[1];
rz(-2.6639498) q[1];
rz(-pi) q[2];
rz(2.4960732) q[3];
sx q[3];
rz(-1.1214167) q[3];
sx q[3];
rz(1.1741032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6496868) q[2];
sx q[2];
rz(-2.457149) q[2];
sx q[2];
rz(0.41927949) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(1.7694582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.075031) q[0];
sx q[0];
rz(-0.12396585) q[0];
sx q[0];
rz(-0.01874622) q[0];
rz(-0.092451267) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(3.0466373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43272797) q[0];
sx q[0];
rz(-1.0148541) q[0];
sx q[0];
rz(-1.1880607) q[0];
rz(-pi) q[1];
rz(0.59573458) q[2];
sx q[2];
rz(-1.2168665) q[2];
sx q[2];
rz(0.46094337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5508037) q[1];
sx q[1];
rz(-1.1491286) q[1];
sx q[1];
rz(-1.8487537) q[1];
x q[2];
rz(-1.2025039) q[3];
sx q[3];
rz(-2.7164258) q[3];
sx q[3];
rz(1.7452546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2350601) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(0.85901421) q[2];
rz(-2.5281995) q[3];
sx q[3];
rz(-0.20039138) q[3];
sx q[3];
rz(-1.8092417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61619806) q[0];
sx q[0];
rz(-1.8680251) q[0];
sx q[0];
rz(1.4816351) q[0];
rz(-2.0741277) q[1];
sx q[1];
rz(-2.4165202) q[1];
sx q[1];
rz(1.3190837) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6292012) q[0];
sx q[0];
rz(-0.90921697) q[0];
sx q[0];
rz(0.11920905) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5379792) q[2];
sx q[2];
rz(-2.5381662) q[2];
sx q[2];
rz(-2.6923384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47506491) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(-0.98824309) q[1];
rz(1.4433925) q[3];
sx q[3];
rz(-1.9249328) q[3];
sx q[3];
rz(-2.4417801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8760425) q[2];
sx q[2];
rz(-2.7713113) q[2];
sx q[2];
rz(3.1004356) q[2];
rz(2.0187078) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48162833) q[0];
sx q[0];
rz(-0.97844231) q[0];
sx q[0];
rz(-1.4725257) q[0];
rz(-2.4166079) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(0.57458007) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1102421) q[0];
sx q[0];
rz(-1.6185805) q[0];
sx q[0];
rz(0.33804531) q[0];
x q[1];
rz(-0.2821183) q[2];
sx q[2];
rz(-1.1638255) q[2];
sx q[2];
rz(-0.59772432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2255937) q[1];
sx q[1];
rz(-0.8527506) q[1];
sx q[1];
rz(-1.6101933) q[1];
x q[2];
rz(-0.88317849) q[3];
sx q[3];
rz(-1.6856582) q[3];
sx q[3];
rz(2.3695947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64728111) q[2];
sx q[2];
rz(-0.24792555) q[2];
sx q[2];
rz(0.9606804) q[2];
rz(-0.48323768) q[3];
sx q[3];
rz(-1.5471231) q[3];
sx q[3];
rz(-1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305337) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(-2.6997987) q[0];
rz(-2.4590625) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(-2.0880879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2164648) q[0];
sx q[0];
rz(-1.3138714) q[0];
sx q[0];
rz(0.83576699) q[0];
rz(-pi) q[1];
rz(1.1978713) q[2];
sx q[2];
rz(-1.4879543) q[2];
sx q[2];
rz(-1.5035455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5610841) q[1];
sx q[1];
rz(-1.2496767) q[1];
sx q[1];
rz(-2.2301484) q[1];
x q[2];
rz(-2.9656762) q[3];
sx q[3];
rz(-2.7421167) q[3];
sx q[3];
rz(-0.3829869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(-2.919) q[2];
rz(1.3709204) q[3];
sx q[3];
rz(-1.7226847) q[3];
sx q[3];
rz(2.5194949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.901578) q[0];
sx q[0];
rz(-1.0116901) q[0];
sx q[0];
rz(-2.1160175) q[0];
rz(1.9758132) q[1];
sx q[1];
rz(-2.5318601) q[1];
sx q[1];
rz(-0.21808521) q[1];
rz(-1.0819306) q[2];
sx q[2];
rz(-2.0578544) q[2];
sx q[2];
rz(-0.73357972) q[2];
rz(-1.2491741) q[3];
sx q[3];
rz(-1.5703441) q[3];
sx q[3];
rz(2.3342568) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
