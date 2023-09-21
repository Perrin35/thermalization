OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53579522) q[0];
sx q[0];
rz(-1.3942766) q[0];
sx q[0];
rz(-1.4730886) q[0];
x q[1];
rz(-2.3660907) q[2];
sx q[2];
rz(-1.3368133) q[2];
sx q[2];
rz(0.50228679) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8625496) q[1];
sx q[1];
rz(-1.057813) q[1];
sx q[1];
rz(1.1898477) q[1];
x q[2];
rz(1.8331068) q[3];
sx q[3];
rz(-0.17671083) q[3];
sx q[3];
rz(1.7929329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9464232) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(2.3852824) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3023892) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(-0.50239262) q[1];
sx q[1];
rz(-0.97351176) q[1];
sx q[1];
rz(-1.5997255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8097336) q[0];
sx q[0];
rz(-2.3947358) q[0];
sx q[0];
rz(-2.5144308) q[0];
rz(1.1142271) q[2];
sx q[2];
rz(-1.4230939) q[2];
sx q[2];
rz(-2.4545124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7911885) q[1];
sx q[1];
rz(-1.5608619) q[1];
sx q[1];
rz(-0.03573461) q[1];
rz(-0.76962556) q[3];
sx q[3];
rz(-2.5105021) q[3];
sx q[3];
rz(1.0069932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(-1.4146457) q[2];
rz(-0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(-1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8787815) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(0.99197018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85342825) q[0];
sx q[0];
rz(-0.58840226) q[0];
sx q[0];
rz(1.2342274) q[0];
x q[1];
rz(-1.7730373) q[2];
sx q[2];
rz(-0.91245203) q[2];
sx q[2];
rz(-0.973268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5132644) q[1];
sx q[1];
rz(-2.6918415) q[1];
sx q[1];
rz(-0.20723923) q[1];
x q[2];
rz(0.46996689) q[3];
sx q[3];
rz(-2.6372058) q[3];
sx q[3];
rz(-0.72987635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38301864) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-0.26829159) q[2];
rz(-2.7475083) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.85369337) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(0.40507856) q[0];
rz(0.69008094) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(-0.69782034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21701248) q[0];
sx q[0];
rz(-0.092886535) q[0];
sx q[0];
rz(-1.3148944) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5562416) q[2];
sx q[2];
rz(-2.5132837) q[2];
sx q[2];
rz(1.5141687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3541672) q[1];
sx q[1];
rz(-2.0560871) q[1];
sx q[1];
rz(0.20809681) q[1];
rz(-pi) q[2];
rz(2.1711369) q[3];
sx q[3];
rz(-2.6451023) q[3];
sx q[3];
rz(-2.0127726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(1.6323803) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(1.0255381) q[0];
rz(2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(0.62932032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51879085) q[0];
sx q[0];
rz(-0.51827058) q[0];
sx q[0];
rz(-0.63664125) q[0];
rz(2.2467381) q[2];
sx q[2];
rz(-1.370315) q[2];
sx q[2];
rz(1.5948053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3931261) q[1];
sx q[1];
rz(-1.4561635) q[1];
sx q[1];
rz(2.6266891) q[1];
rz(-pi) q[2];
rz(-2.1719645) q[3];
sx q[3];
rz(-1.3519577) q[3];
sx q[3];
rz(2.9588483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(1.7957548) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(-0.19601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(-1.7350896) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.7746183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.470984) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(-1.8355808) q[0];
x q[1];
rz(-2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(2.218354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1320912) q[1];
sx q[1];
rz(-0.17772929) q[1];
sx q[1];
rz(1.0264261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9983415) q[3];
sx q[3];
rz(-1.0435836) q[3];
sx q[3];
rz(-2.9851819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24484816) q[2];
sx q[2];
rz(-2.5881793) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.7287438) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(-0.81378716) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.3758804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8517075) q[0];
sx q[0];
rz(-0.87806784) q[0];
sx q[0];
rz(-0.38498199) q[0];
rz(1.4543578) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(1.6778698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.059767698) q[1];
sx q[1];
rz(-1.6762814) q[1];
sx q[1];
rz(-3.0256773) q[1];
rz(-pi) q[2];
rz(1.4010299) q[3];
sx q[3];
rz(-0.3347291) q[3];
sx q[3];
rz(1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(-2.4105371) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(-2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(-2.5336174) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(0.2342934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.308134) q[0];
sx q[0];
rz(-1.3123543) q[0];
sx q[0];
rz(1.0779557) q[0];
rz(-pi) q[1];
rz(-0.72006165) q[2];
sx q[2];
rz(-2.10663) q[2];
sx q[2];
rz(2.5502887) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.329511) q[1];
sx q[1];
rz(-0.78204621) q[1];
sx q[1];
rz(-1.8409607) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33741823) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(-2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12604788) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.4253433) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.5995246) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(1.19338) q[0];
rz(-1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017100447) q[0];
sx q[0];
rz(-1.955535) q[0];
sx q[0];
rz(0.7229294) q[0];
rz(-1.1759042) q[2];
sx q[2];
rz(-1.8443622) q[2];
sx q[2];
rz(3.0964031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8972842) q[1];
sx q[1];
rz(-0.329204) q[1];
sx q[1];
rz(-2.5876178) q[1];
x q[2];
rz(1.0778905) q[3];
sx q[3];
rz(-2.7422046) q[3];
sx q[3];
rz(-2.3139017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(0.28820583) q[2];
rz(-0.47973412) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(-1.5079927) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(2.7669725) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3044395) q[0];
sx q[0];
rz(-1.8189948) q[0];
sx q[0];
rz(0.34557839) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.183379) q[2];
sx q[2];
rz(-1.3365067) q[2];
sx q[2];
rz(2.2189552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76162321) q[1];
sx q[1];
rz(-0.39561158) q[1];
sx q[1];
rz(0.53669866) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6463514) q[3];
sx q[3];
rz(-1.75845) q[3];
sx q[3];
rz(-0.56896602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23641071) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(-2.6399844) q[2];
rz(-1.2891399) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(0.44617173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63507737) q[0];
sx q[0];
rz(-1.7000533) q[0];
sx q[0];
rz(0.62361367) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(1.7239465) q[2];
sx q[2];
rz(-1.3071123) q[2];
sx q[2];
rz(1.2109962) q[2];
rz(-0.049384762) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
