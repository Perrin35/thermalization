OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(0.88357893) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53579522) q[0];
sx q[0];
rz(-1.3942766) q[0];
sx q[0];
rz(1.4730886) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3281524) q[2];
sx q[2];
rz(-0.80291623) q[2];
sx q[2];
rz(0.83628718) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27904305) q[1];
sx q[1];
rz(-1.057813) q[1];
sx q[1];
rz(-1.1898477) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3084859) q[3];
sx q[3];
rz(-0.17671083) q[3];
sx q[3];
rz(1.7929329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(0.26120734) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(2.3852824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(0.38683495) q[0];
rz(-0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5997255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889266) q[0];
sx q[0];
rz(-0.98836556) q[0];
sx q[0];
rz(2.0684588) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16427152) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(2.3300366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6501573) q[1];
sx q[1];
rz(-3.1045034) q[1];
sx q[1];
rz(-2.8703719) q[1];
rz(-pi) q[2];
rz(-0.76962556) q[3];
sx q[3];
rz(-0.63109055) q[3];
sx q[3];
rz(-1.0069932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.60275045) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(-1.4146457) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(-1.458228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8787815) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(0.61022726) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(-2.1496225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45527601) q[0];
sx q[0];
rz(-2.1222097) q[0];
sx q[0];
rz(-2.9247012) q[0];
rz(-pi) q[1];
x q[1];
rz(2.473258) q[2];
sx q[2];
rz(-1.4112345) q[2];
sx q[2];
rz(-0.72232407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85769535) q[1];
sx q[1];
rz(-1.1313492) q[1];
sx q[1];
rz(1.6698014) q[1];
rz(0.45735995) q[3];
sx q[3];
rz(-1.7914346) q[3];
sx q[3];
rz(-1.8822576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(2.8733011) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(0.18850732) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(0.69008094) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(-0.69782034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989379) q[0];
sx q[0];
rz(-1.5473167) q[0];
sx q[0];
rz(-1.4809181) q[0];
x q[1];
rz(-2.5562416) q[2];
sx q[2];
rz(-0.62830892) q[2];
sx q[2];
rz(1.6274239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78742541) q[1];
sx q[1];
rz(-2.0560871) q[1];
sx q[1];
rz(-0.20809681) q[1];
rz(-0.97045578) q[3];
sx q[3];
rz(-0.49649039) q[3];
sx q[3];
rz(2.0127726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(1.6323803) q[2];
rz(2.737282) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(2.1160545) q[0];
rz(-2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(-0.62932032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6228018) q[0];
sx q[0];
rz(-2.6233221) q[0];
sx q[0];
rz(-2.5049514) q[0];
x q[1];
rz(1.8848558) q[2];
sx q[2];
rz(-0.70054189) q[2];
sx q[2];
rz(-2.8741921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1197966) q[1];
sx q[1];
rz(-2.6152059) q[1];
sx q[1];
rz(-2.9119133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26341565) q[3];
sx q[3];
rz(-2.1556971) q[3];
sx q[3];
rz(1.5358621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(-1.3458378) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(-0.19601823) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(2.956399) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(-1.3669744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0104116) q[0];
sx q[0];
rz(-0.62958065) q[0];
sx q[0];
rz(-2.7601943) q[0];
rz(2.30079) q[2];
sx q[2];
rz(-1.3422988) q[2];
sx q[2];
rz(-0.92323869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5806611) q[1];
sx q[1];
rz(-1.7226189) q[1];
sx q[1];
rz(0.092756943) q[1];
rz(-pi) q[2];
x q[2];
rz(2.572445) q[3];
sx q[3];
rz(-1.2043118) q[3];
sx q[3];
rz(-1.6397427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(2.2582167) q[2];
rz(1.7287438) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(-0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(3.1037519) q[1];
sx q[1];
rz(-0.81532878) q[1];
sx q[1];
rz(-1.3758804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517075) q[0];
sx q[0];
rz(-0.87806784) q[0];
sx q[0];
rz(-2.7566107) q[0];
rz(-pi) q[1];
rz(0.059869754) q[2];
sx q[2];
rz(-2.0436358) q[2];
sx q[2];
rz(-1.5945895) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.059767698) q[1];
sx q[1];
rz(-1.4653112) q[1];
sx q[1];
rz(0.11591537) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2405346) q[3];
sx q[3];
rz(-1.6263279) q[3];
sx q[3];
rz(-0.33952573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(-2.4105371) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(0.73356432) q[0];
rz(-0.60797524) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83345862) q[0];
sx q[0];
rz(-1.3123543) q[0];
sx q[0];
rz(1.0779557) q[0];
x q[1];
rz(-2.2393164) q[2];
sx q[2];
rz(-0.96792816) q[2];
sx q[2];
rz(0.55842802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1842321) q[1];
sx q[1];
rz(-0.82416526) q[1];
sx q[1];
rz(2.8824473) q[1];
rz(-pi) q[2];
rz(-2.8041744) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(-0.19677256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(-1.7162494) q[2];
rz(1.6783293) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(1.19338) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(0.9448005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1244922) q[0];
sx q[0];
rz(-1.955535) q[0];
sx q[0];
rz(-0.7229294) q[0];
rz(-pi) q[1];
rz(-1.9656885) q[2];
sx q[2];
rz(-1.8443622) q[2];
sx q[2];
rz(0.045189518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2855125) q[1];
sx q[1];
rz(-1.3998919) q[1];
sx q[1];
rz(0.28275615) q[1];
rz(-1.9267843) q[3];
sx q[3];
rz(-1.3857406) q[3];
sx q[3];
rz(-1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21645674) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(0.47973412) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.5079927) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(2.2286041) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3044395) q[0];
sx q[0];
rz(-1.8189948) q[0];
sx q[0];
rz(-0.34557839) q[0];
x q[1];
rz(0.28384039) q[2];
sx q[2];
rz(-2.1643057) q[2];
sx q[2];
rz(2.3317091) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8303463) q[1];
sx q[1];
rz(-1.7691358) q[1];
sx q[1];
rz(-2.7970008) q[1];
rz(-pi) q[2];
rz(1.6463514) q[3];
sx q[3];
rz(-1.3831426) q[3];
sx q[3];
rz(0.56896602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(-0.44617173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-1.1322017) q[1];
sx q[1];
rz(-2.3846346) q[1];
sx q[1];
rz(0.089288575) q[1];
rz(-1.7239465) q[2];
sx q[2];
rz(-1.8344804) q[2];
sx q[2];
rz(-1.9305965) q[2];
rz(0.049384762) q[3];
sx q[3];
rz(-1.0247083) q[3];
sx q[3];
rz(0.68030737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];