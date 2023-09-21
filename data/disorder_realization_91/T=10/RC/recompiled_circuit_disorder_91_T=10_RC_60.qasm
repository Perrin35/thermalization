OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(5.7167238) q[0];
sx q[0];
rz(9.5958435) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(-1.3309825) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(1.5050423) q[0];
rz(-pi) q[1];
rz(2.6257329) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(-3.0992103) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-0.25169262) q[1];
sx q[1];
rz(-1.6394872) q[1];
rz(2.4986476) q[3];
sx q[3];
rz(-1.7503947) q[3];
sx q[3];
rz(0.56550607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(-1.7261516) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.4527028) q[0];
rz(-0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6712924) q[0];
sx q[0];
rz(-1.6912582) q[0];
sx q[0];
rz(2.5380773) q[0];
x q[1];
rz(0.47313182) q[2];
sx q[2];
rz(-2.5169249) q[2];
sx q[2];
rz(1.9249137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86094942) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(1.8722948) q[1];
rz(-pi) q[2];
rz(-1.3817915) q[3];
sx q[3];
rz(-1.731589) q[3];
sx q[3];
rz(0.27634987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(-2.753479) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(-0.31633502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(1.578376) q[0];
rz(-pi) q[1];
rz(1.146831) q[2];
sx q[2];
rz(-2.1018873) q[2];
sx q[2];
rz(2.6742427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2410779) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(0.48836744) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9100788) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(-0.85698444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(-1.0236615) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(1.3446993) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(-1.6279189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154402) q[0];
sx q[0];
rz(-2.2427796) q[0];
sx q[0];
rz(1.5120904) q[0];
rz(1.4430181) q[2];
sx q[2];
rz(-2.5737692) q[2];
sx q[2];
rz(-2.36433) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0107683) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(-0.018647714) q[1];
rz(-pi) q[2];
rz(1.9862595) q[3];
sx q[3];
rz(-0.88298015) q[3];
sx q[3];
rz(-2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(1.9173737) q[2];
rz(0.41695693) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(-1.303724) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(3.0900893) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52841016) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(-0.86092237) q[0];
rz(-0.15437834) q[2];
sx q[2];
rz(-1.8533684) q[2];
sx q[2];
rz(-2.4992361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4645849) q[1];
sx q[1];
rz(-1.1255463) q[1];
sx q[1];
rz(0.17226179) q[1];
rz(-2.4078835) q[3];
sx q[3];
rz(-1.4791616) q[3];
sx q[3];
rz(1.1250145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14295386) q[0];
sx q[0];
rz(-2.8790701) q[0];
sx q[0];
rz(0.72857626) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7716878) q[2];
sx q[2];
rz(-2.7668014) q[2];
sx q[2];
rz(-2.7704266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-0.38227907) q[1];
sx q[1];
rz(-1.2804968) q[1];
x q[2];
rz(1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.3827682) q[2];
rz(-2.8921195) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(2.2447341) q[0];
rz(0.25009051) q[1];
sx q[1];
rz(-2.9607594) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8539124) q[0];
sx q[0];
rz(-2.2932862) q[0];
sx q[0];
rz(1.7476837) q[0];
rz(2.2824077) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(-0.32260103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84278216) q[1];
sx q[1];
rz(-2.543445) q[1];
sx q[1];
rz(0.36069718) q[1];
x q[2];
rz(1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(-0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(0.21952595) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-0.21729939) q[0];
rz(0.030933881) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(-2.4826179) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156408) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(0.036235972) q[0];
x q[1];
rz(2.4023513) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(3.0657363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3498889) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(1.2757343) q[1];
rz(2.213845) q[3];
sx q[3];
rz(-1.6258996) q[3];
sx q[3];
rz(2.7411214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.9314996) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3037198) q[0];
sx q[0];
rz(-1.923773) q[0];
sx q[0];
rz(2.9205802) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25173431) q[2];
sx q[2];
rz(-2.1605957) q[2];
sx q[2];
rz(0.41433197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71483892) q[1];
sx q[1];
rz(-0.73298448) q[1];
sx q[1];
rz(2.0183802) q[1];
x q[2];
rz(-2.5217767) q[3];
sx q[3];
rz(-1.2634988) q[3];
sx q[3];
rz(-2.4710771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.575763) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(2.2120655) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(-1.7738336) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(2.3419103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4096217) q[0];
sx q[0];
rz(-0.35847607) q[0];
sx q[0];
rz(1.326807) q[0];
rz(-pi) q[1];
rz(0.16074796) q[2];
sx q[2];
rz(-1.4904009) q[2];
sx q[2];
rz(-0.71842566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3696246) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(-0.97009138) q[1];
x q[2];
rz(-2.9959216) q[3];
sx q[3];
rz(-1.4958069) q[3];
sx q[3];
rz(-0.86736995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-0.69520673) q[2];
rz(0.55784145) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.862539) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(0.002481133) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(0.1602608) q[3];
sx q[3];
rz(-0.80453034) q[3];
sx q[3];
rz(2.9321032) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
