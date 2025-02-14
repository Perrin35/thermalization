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
rz(2.1386327) q[0];
sx q[0];
rz(-0.47337368) q[0];
sx q[0];
rz(2.0546497) q[0];
rz(-0.76397693) q[1];
sx q[1];
rz(-0.44372258) q[1];
sx q[1];
rz(-0.8313764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2646541) q[0];
sx q[0];
rz(-1.042141) q[0];
sx q[0];
rz(-2.1318046) q[0];
x q[1];
rz(2.0152807) q[2];
sx q[2];
rz(-0.51156564) q[2];
sx q[2];
rz(-1.4064521) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2519166) q[1];
sx q[1];
rz(-1.6908501) q[1];
sx q[1];
rz(-3.037583) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0696297) q[3];
sx q[3];
rz(-1.5441893) q[3];
sx q[3];
rz(-0.047544971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5272556) q[2];
sx q[2];
rz(-1.273512) q[2];
sx q[2];
rz(0.41198507) q[2];
rz(-0.24474239) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.7016542) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(-0.42020759) q[0];
rz(-0.48928753) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(1.18321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91285465) q[0];
sx q[0];
rz(-1.1381554) q[0];
sx q[0];
rz(-0.88254024) q[0];
rz(-pi) q[1];
rz(2.3406896) q[2];
sx q[2];
rz(-2.704853) q[2];
sx q[2];
rz(0.85608053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1258613) q[1];
sx q[1];
rz(-0.9974879) q[1];
sx q[1];
rz(-1.1124951) q[1];
x q[2];
rz(-1.7443827) q[3];
sx q[3];
rz(-1.4549482) q[3];
sx q[3];
rz(2.4811452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45541397) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(0.14032826) q[2];
rz(2.8651107) q[3];
sx q[3];
rz(-2.3639207) q[3];
sx q[3];
rz(-0.28237835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.33271933) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(-1.7669539) q[0];
rz(-2.2287492) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(-1.0309781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46488097) q[0];
sx q[0];
rz(-1.8959989) q[0];
sx q[0];
rz(0.19285481) q[0];
rz(-pi) q[1];
rz(-2.7430277) q[2];
sx q[2];
rz(-1.0462648) q[2];
sx q[2];
rz(2.2664859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1304969) q[1];
sx q[1];
rz(-1.5807932) q[1];
sx q[1];
rz(1.1586055) q[1];
rz(1.3783021) q[3];
sx q[3];
rz(-1.2519426) q[3];
sx q[3];
rz(1.1553193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-0.63808179) q[2];
rz(3.0408527) q[3];
sx q[3];
rz(-2.1295857) q[3];
sx q[3];
rz(0.49873275) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8345555) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(-1.3822973) q[0];
rz(0.31940976) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(-2.3841948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46947458) q[0];
sx q[0];
rz(-1.645477) q[0];
sx q[0];
rz(1.5622219) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4057233) q[2];
sx q[2];
rz(-2.9506548) q[2];
sx q[2];
rz(1.5519993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0075052) q[1];
sx q[1];
rz(-1.1267091) q[1];
sx q[1];
rz(0.92625463) q[1];
rz(-pi) q[2];
rz(1.0469466) q[3];
sx q[3];
rz(-1.0312087) q[3];
sx q[3];
rz(1.3803409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3697529) q[2];
sx q[2];
rz(-0.87551337) q[2];
sx q[2];
rz(-2.3337951) q[2];
rz(1.7317023) q[3];
sx q[3];
rz(-1.8818703) q[3];
sx q[3];
rz(-0.65214777) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-2.763651) q[0];
sx q[0];
rz(-2.156303) q[0];
rz(-1.6458884) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(-0.92934242) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791361) q[0];
sx q[0];
rz(-2.1497288) q[0];
sx q[0];
rz(1.2499362) q[0];
rz(-0.55589755) q[2];
sx q[2];
rz(-0.68624845) q[2];
sx q[2];
rz(0.85418073) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29581901) q[1];
sx q[1];
rz(-0.56987485) q[1];
sx q[1];
rz(-1.1954346) q[1];
rz(-3.1272917) q[3];
sx q[3];
rz(-1.0096483) q[3];
sx q[3];
rz(0.98706571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0890961) q[2];
sx q[2];
rz(-1.3546983) q[2];
sx q[2];
rz(-2.055114) q[2];
rz(0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(1.121678) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97359598) q[0];
sx q[0];
rz(-0.43742988) q[0];
sx q[0];
rz(-2.9214389) q[0];
rz(-1.6379697) q[1];
sx q[1];
rz(-1.3238246) q[1];
sx q[1];
rz(0.55353177) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65992113) q[0];
sx q[0];
rz(-1.685647) q[0];
sx q[0];
rz(2.9433795) q[0];
rz(-1.7479182) q[2];
sx q[2];
rz(-1.2125041) q[2];
sx q[2];
rz(-0.19233335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0281535) q[1];
sx q[1];
rz(-2.4467797) q[1];
sx q[1];
rz(-3.0656612) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8643478) q[3];
sx q[3];
rz(-2.436815) q[3];
sx q[3];
rz(2.1516906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0452051) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(-0.20827797) q[2];
rz(-1.1194718) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.646362) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(-2.0679423) q[0];
rz(0.13058361) q[1];
sx q[1];
rz(-1.5740296) q[1];
sx q[1];
rz(2.1317587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45859858) q[0];
sx q[0];
rz(-1.7560335) q[0];
sx q[0];
rz(-2.7390929) q[0];
rz(-0.45744894) q[2];
sx q[2];
rz(-0.1651131) q[2];
sx q[2];
rz(0.091006309) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9562137) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(0.57422178) q[1];
rz(-2.6920986) q[3];
sx q[3];
rz(-2.426034) q[3];
sx q[3];
rz(0.64330949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69663557) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(-2.8182287) q[2];
rz(2.4671386) q[3];
sx q[3];
rz(-2.2489397) q[3];
sx q[3];
rz(0.023580624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0167639) q[0];
sx q[0];
rz(-0.49377307) q[0];
sx q[0];
rz(2.0523409) q[0];
rz(-2.543154) q[1];
sx q[1];
rz(-1.8611703) q[1];
sx q[1];
rz(0.79900297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625811) q[0];
sx q[0];
rz(-2.7220753) q[0];
sx q[0];
rz(-0.74504344) q[0];
rz(3.0472894) q[2];
sx q[2];
rz(-0.906773) q[2];
sx q[2];
rz(1.9349328) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0388698) q[1];
sx q[1];
rz(-0.33867044) q[1];
sx q[1];
rz(2.6790957) q[1];
x q[2];
rz(-1.5653506) q[3];
sx q[3];
rz(-1.3289691) q[3];
sx q[3];
rz(-2.8018045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66990718) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(-1.8890107) q[2];
rz(-1.0505098) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(2.2091776) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85225409) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(-2.6259212) q[0];
rz(-1.6853261) q[1];
sx q[1];
rz(-1.8003502) q[1];
sx q[1];
rz(-0.32404831) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1103678) q[0];
sx q[0];
rz(-0.47981167) q[0];
sx q[0];
rz(-0.57439248) q[0];
x q[1];
rz(-2.134974) q[2];
sx q[2];
rz(-0.7340275) q[2];
sx q[2];
rz(3.1285518) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9726213) q[1];
sx q[1];
rz(-1.7105127) q[1];
sx q[1];
rz(-2.0129544) q[1];
x q[2];
rz(1.0131272) q[3];
sx q[3];
rz(-1.5032167) q[3];
sx q[3];
rz(0.24012071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9825762) q[2];
sx q[2];
rz(-2.0197208) q[2];
sx q[2];
rz(0.70875657) q[2];
rz(-2.3580264) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.6243352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5358955) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(0.15644431) q[0];
rz(2.5018196) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(0.99008647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.472725) q[0];
sx q[0];
rz(-0.91356431) q[0];
sx q[0];
rz(-0.79909493) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.089870139) q[2];
sx q[2];
rz(-1.9192358) q[2];
sx q[2];
rz(2.9065098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53452021) q[1];
sx q[1];
rz(-1.4379518) q[1];
sx q[1];
rz(-2.584949) q[1];
rz(-pi) q[2];
rz(-2.5699432) q[3];
sx q[3];
rz(-0.23861966) q[3];
sx q[3];
rz(-2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(-0.61335316) q[2];
rz(-0.26398811) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(3.0494674) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(-1.4641948) q[1];
sx q[1];
rz(-0.72036998) q[1];
sx q[1];
rz(1.5710685) q[1];
rz(-1.175066) q[2];
sx q[2];
rz(-1.1968975) q[2];
sx q[2];
rz(0.71002985) q[2];
rz(-1.6102552) q[3];
sx q[3];
rz(-0.95245556) q[3];
sx q[3];
rz(0.34921619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
