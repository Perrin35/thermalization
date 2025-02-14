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
rz(0.60628211) q[0];
sx q[0];
rz(-2.0908794) q[0];
sx q[0];
rz(2.3442955) q[0];
rz(-2.1939313) q[1];
sx q[1];
rz(5.0343577) q[1];
sx q[1];
rz(9.0696637) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712747) q[0];
sx q[0];
rz(-1.6375033) q[0];
sx q[0];
rz(2.8913777) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0933286) q[2];
sx q[2];
rz(-1.588378) q[2];
sx q[2];
rz(-0.64811743) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15820643) q[1];
sx q[1];
rz(-2.4987028) q[1];
sx q[1];
rz(-0.41916533) q[1];
rz(0.83824165) q[3];
sx q[3];
rz(-0.62617338) q[3];
sx q[3];
rz(-0.29919344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.015410272) q[2];
sx q[2];
rz(-2.5842032) q[2];
sx q[2];
rz(-1.9554479) q[2];
rz(1.0960389) q[3];
sx q[3];
rz(-2.3604184) q[3];
sx q[3];
rz(-1.5196777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1315411) q[0];
sx q[0];
rz(-3.0942656) q[0];
sx q[0];
rz(1.9897687) q[0];
rz(0.25543073) q[1];
sx q[1];
rz(-1.3574182) q[1];
sx q[1];
rz(-0.40886042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9329703) q[0];
sx q[0];
rz(-1.044036) q[0];
sx q[0];
rz(2.7072858) q[0];
x q[1];
rz(-1.4299655) q[2];
sx q[2];
rz(-0.69770506) q[2];
sx q[2];
rz(-1.0753466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3656816) q[1];
sx q[1];
rz(-1.3789007) q[1];
sx q[1];
rz(0.86071162) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71315546) q[3];
sx q[3];
rz(-0.9760614) q[3];
sx q[3];
rz(0.99280533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6795808) q[2];
sx q[2];
rz(-2.2363594) q[2];
sx q[2];
rz(-2.8751664) q[2];
rz(-2.5816495) q[3];
sx q[3];
rz(-2.4661049) q[3];
sx q[3];
rz(2.717836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4567081) q[0];
sx q[0];
rz(-2.279778) q[0];
sx q[0];
rz(-1.7488712) q[0];
rz(-2.4640153) q[1];
sx q[1];
rz(-1.1303439) q[1];
sx q[1];
rz(-0.63634253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0746431) q[0];
sx q[0];
rz(-2.136629) q[0];
sx q[0];
rz(-1.1003347) q[0];
x q[1];
rz(-0.26480953) q[2];
sx q[2];
rz(-2.2356459) q[2];
sx q[2];
rz(-0.88685551) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96752466) q[1];
sx q[1];
rz(-1.2651769) q[1];
sx q[1];
rz(0.6966865) q[1];
rz(-pi) q[2];
rz(2.3597765) q[3];
sx q[3];
rz(-1.86275) q[3];
sx q[3];
rz(1.7192732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59844083) q[2];
sx q[2];
rz(-1.8676912) q[2];
sx q[2];
rz(0.72371036) q[2];
rz(1.9321457) q[3];
sx q[3];
rz(-0.8330605) q[3];
sx q[3];
rz(-2.3256653) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4010312) q[0];
sx q[0];
rz(-2.7771692) q[0];
sx q[0];
rz(0.60582274) q[0];
rz(-1.1653028) q[1];
sx q[1];
rz(-2.0348246) q[1];
sx q[1];
rz(-2.5071526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0710356) q[0];
sx q[0];
rz(-1.0939176) q[0];
sx q[0];
rz(1.5232851) q[0];
x q[1];
rz(2.6003378) q[2];
sx q[2];
rz(-1.4536023) q[2];
sx q[2];
rz(1.8126496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2812579) q[1];
sx q[1];
rz(-2.0678646) q[1];
sx q[1];
rz(-0.71038891) q[1];
x q[2];
rz(-0.090417763) q[3];
sx q[3];
rz(-1.8274698) q[3];
sx q[3];
rz(0.21714563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4623462) q[2];
sx q[2];
rz(-0.99455589) q[2];
sx q[2];
rz(-1.4086949) q[2];
rz(-2.5016221) q[3];
sx q[3];
rz(-0.43158117) q[3];
sx q[3];
rz(1.4225167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.53936476) q[0];
sx q[0];
rz(-0.27502763) q[0];
sx q[0];
rz(-1.6590903) q[0];
rz(-1.9954584) q[1];
sx q[1];
rz(-2.4303552) q[1];
sx q[1];
rz(-1.8123288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70721165) q[0];
sx q[0];
rz(-2.2298621) q[0];
sx q[0];
rz(0.7313318) q[0];
rz(0.48945697) q[2];
sx q[2];
rz(-2.202919) q[2];
sx q[2];
rz(2.7165627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0983713) q[1];
sx q[1];
rz(-0.30128208) q[1];
sx q[1];
rz(-2.2194018) q[1];
x q[2];
rz(3.077636) q[3];
sx q[3];
rz(-1.8956479) q[3];
sx q[3];
rz(-0.16983804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9596404) q[2];
sx q[2];
rz(-0.96936172) q[2];
sx q[2];
rz(2.1592965) q[2];
rz(2.9082409) q[3];
sx q[3];
rz(-1.8459277) q[3];
sx q[3];
rz(2.6006234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(3.0503466) q[0];
sx q[0];
rz(-2.8627658) q[0];
sx q[0];
rz(0.74063754) q[0];
rz(-0.29847538) q[1];
sx q[1];
rz(-1.1497755) q[1];
sx q[1];
rz(0.99075913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753526) q[0];
sx q[0];
rz(-1.9431433) q[0];
sx q[0];
rz(1.5840992) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58133881) q[2];
sx q[2];
rz(-0.55989385) q[2];
sx q[2];
rz(-2.8357504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3988855) q[1];
sx q[1];
rz(-1.1749335) q[1];
sx q[1];
rz(2.2363801) q[1];
rz(1.1634109) q[3];
sx q[3];
rz(-2.2797027) q[3];
sx q[3];
rz(1.8590301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52567452) q[2];
sx q[2];
rz(-2.2423223) q[2];
sx q[2];
rz(0.57596469) q[2];
rz(1.2032262) q[3];
sx q[3];
rz(-1.7267745) q[3];
sx q[3];
rz(-1.6171713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.7791157) q[0];
sx q[0];
rz(-0.285382) q[0];
sx q[0];
rz(-1.2364291) q[0];
rz(1.7615039) q[1];
sx q[1];
rz(-0.80685902) q[1];
sx q[1];
rz(-2.7469452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98251909) q[0];
sx q[0];
rz(-0.22501105) q[0];
sx q[0];
rz(-1.7270442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.790843) q[2];
sx q[2];
rz(-1.4119822) q[2];
sx q[2];
rz(-1.0823298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9731991) q[1];
sx q[1];
rz(-1.7500119) q[1];
sx q[1];
rz(0.94931305) q[1];
rz(-pi) q[2];
rz(-2.8498185) q[3];
sx q[3];
rz(-1.7543766) q[3];
sx q[3];
rz(-2.1734299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3820162) q[2];
sx q[2];
rz(-2.3754109) q[2];
sx q[2];
rz(1.3845328) q[2];
rz(-1.1197439) q[3];
sx q[3];
rz(-1.1625959) q[3];
sx q[3];
rz(-2.3118481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2914332) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(-0.39655381) q[0];
rz(1.5810168) q[1];
sx q[1];
rz(-2.719559) q[1];
sx q[1];
rz(0.54403725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85910857) q[0];
sx q[0];
rz(-0.73900925) q[0];
sx q[0];
rz(-3.0942731) q[0];
x q[1];
rz(-1.3048473) q[2];
sx q[2];
rz(-1.6839993) q[2];
sx q[2];
rz(1.2987657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8072784) q[1];
sx q[1];
rz(-2.6014666) q[1];
sx q[1];
rz(-0.29088104) q[1];
rz(-pi) q[2];
rz(1.7782147) q[3];
sx q[3];
rz(-0.44379297) q[3];
sx q[3];
rz(1.7621413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80514002) q[2];
sx q[2];
rz(-0.84638941) q[2];
sx q[2];
rz(-1.7479755) q[2];
rz(-2.9584068) q[3];
sx q[3];
rz(-1.2510679) q[3];
sx q[3];
rz(-0.92060554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377604) q[0];
sx q[0];
rz(-2.5354009) q[0];
sx q[0];
rz(0.071618557) q[0];
rz(2.3764745) q[1];
sx q[1];
rz(-1.3970929) q[1];
sx q[1];
rz(-2.5535233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4912754) q[0];
sx q[0];
rz(-1.2450019) q[0];
sx q[0];
rz(-0.97401818) q[0];
rz(0.88895615) q[2];
sx q[2];
rz(-1.5247022) q[2];
sx q[2];
rz(2.7323013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32826281) q[1];
sx q[1];
rz(-0.67280233) q[1];
sx q[1];
rz(-0.15209122) q[1];
x q[2];
rz(-1.16251) q[3];
sx q[3];
rz(-2.4823349) q[3];
sx q[3];
rz(-1.6867278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4582943) q[2];
sx q[2];
rz(-2.8647162) q[2];
sx q[2];
rz(1.5923306) q[2];
rz(-1.8706627) q[3];
sx q[3];
rz(-1.9586321) q[3];
sx q[3];
rz(-2.8196238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1147605) q[0];
sx q[0];
rz(-2.8883567) q[0];
sx q[0];
rz(2.7993171) q[0];
rz(0.5510785) q[1];
sx q[1];
rz(-1.2197878) q[1];
sx q[1];
rz(-0.50318998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7320815) q[0];
sx q[0];
rz(-2.4284673) q[0];
sx q[0];
rz(-0.54422191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.25496) q[2];
sx q[2];
rz(-1.6079512) q[2];
sx q[2];
rz(-2.6583918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7552265) q[1];
sx q[1];
rz(-1.0468353) q[1];
sx q[1];
rz(-0.95295277) q[1];
x q[2];
rz(0.14519174) q[3];
sx q[3];
rz(-2.593793) q[3];
sx q[3];
rz(-2.5992268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8503348) q[2];
sx q[2];
rz(-1.0185654) q[2];
sx q[2];
rz(2.2955718) q[2];
rz(-3.1076943) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(0.94023824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.7593183) q[0];
sx q[0];
rz(-1.2713534) q[0];
sx q[0];
rz(-0.8937339) q[0];
rz(-1.8779124) q[1];
sx q[1];
rz(-2.3448941) q[1];
sx q[1];
rz(5/(14*pi)) q[1];
rz(2.0533753) q[2];
sx q[2];
rz(-0.84090424) q[2];
sx q[2];
rz(3.1089208) q[2];
rz(-0.55584021) q[3];
sx q[3];
rz(-1.8636788) q[3];
sx q[3];
rz(2.004971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
