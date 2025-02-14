OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(-1.4077185) q[0];
sx q[0];
rz(-0.044535927) q[0];
rz(2.0308004) q[1];
sx q[1];
rz(-1.9171311) q[1];
sx q[1];
rz(0.78723025) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48101857) q[0];
sx q[0];
rz(-0.74587599) q[0];
sx q[0];
rz(-2.3700299) q[0];
rz(1.8960612) q[2];
sx q[2];
rz(-1.533154) q[2];
sx q[2];
rz(-1.7104488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2582142) q[1];
sx q[1];
rz(-1.5014663) q[1];
sx q[1];
rz(-0.99154559) q[1];
rz(-pi) q[2];
rz(2.1544632) q[3];
sx q[3];
rz(-2.7257724) q[3];
sx q[3];
rz(0.73279954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74701509) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(-0.75418312) q[2];
rz(-1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(-2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873782) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(-1.8074328) q[0];
rz(1.2913903) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(-2.2659567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9459045) q[0];
sx q[0];
rz(-2.2155846) q[0];
sx q[0];
rz(-0.88627215) q[0];
rz(-pi) q[1];
x q[1];
rz(0.062280999) q[2];
sx q[2];
rz(-1.0750293) q[2];
sx q[2];
rz(-0.77046662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7499979) q[1];
sx q[1];
rz(-1.1340967) q[1];
sx q[1];
rz(-2.8207247) q[1];
rz(1.7784714) q[3];
sx q[3];
rz(-1.8559578) q[3];
sx q[3];
rz(-0.36237291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4448173) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(-2.2689421) q[2];
rz(-0.78898346) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(2.9130329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2149684) q[0];
sx q[0];
rz(-1.3525532) q[0];
sx q[0];
rz(0.0080000814) q[0];
rz(-2.3232715) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(-1.2368894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43046666) q[0];
sx q[0];
rz(-0.9014117) q[0];
sx q[0];
rz(2.7276464) q[0];
x q[1];
rz(-2.2358782) q[2];
sx q[2];
rz(-2.6092898) q[2];
sx q[2];
rz(-1.5168845) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9489386) q[1];
sx q[1];
rz(-2.1100419) q[1];
sx q[1];
rz(1.9303028) q[1];
rz(-2.6872356) q[3];
sx q[3];
rz(-2.8792692) q[3];
sx q[3];
rz(1.232805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41039738) q[2];
sx q[2];
rz(-0.13430139) q[2];
sx q[2];
rz(-0.94089874) q[2];
rz(2.0375552) q[3];
sx q[3];
rz(-1.3528115) q[3];
sx q[3];
rz(-1.997939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10821548) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-2.9040842) q[0];
rz(-0.6595276) q[1];
sx q[1];
rz(-2.567465) q[1];
sx q[1];
rz(-3.1245756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68786808) q[0];
sx q[0];
rz(-2.9659418) q[0];
sx q[0];
rz(-1.2916628) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41827664) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(-1.0470225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3238195) q[1];
sx q[1];
rz(-0.57460472) q[1];
sx q[1];
rz(2.9486604) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7155101) q[3];
sx q[3];
rz(-2.1984716) q[3];
sx q[3];
rz(-1.1174517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90618769) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(0.91442937) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-0.80949628) q[3];
sx q[3];
rz(-0.58524281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.888716) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(-0.1121029) q[0];
rz(2.0762699) q[1];
sx q[1];
rz(-1.0708829) q[1];
sx q[1];
rz(-1.0692474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1493268) q[0];
sx q[0];
rz(-0.42851028) q[0];
sx q[0];
rz(-2.4403768) q[0];
rz(-pi) q[1];
rz(-0.5324131) q[2];
sx q[2];
rz(-1.1523917) q[2];
sx q[2];
rz(0.46506986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9578609) q[1];
sx q[1];
rz(-1.2535447) q[1];
sx q[1];
rz(0.31362335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97158708) q[3];
sx q[3];
rz(-1.6808482) q[3];
sx q[3];
rz(-2.9436265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7142882) q[2];
sx q[2];
rz(-2.4236743) q[2];
sx q[2];
rz(2.60738) q[2];
rz(0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(-1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7233906) q[0];
sx q[0];
rz(-1.284282) q[0];
sx q[0];
rz(0.89163017) q[0];
rz(1.3806237) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(-0.92323971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6398481) q[0];
sx q[0];
rz(-0.11127936) q[0];
sx q[0];
rz(0.48868816) q[0];
x q[1];
rz(2.7545287) q[2];
sx q[2];
rz(-1.825807) q[2];
sx q[2];
rz(-2.2789291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2378943) q[1];
sx q[1];
rz(-0.35939068) q[1];
sx q[1];
rz(0.87025799) q[1];
rz(-pi) q[2];
rz(-0.51839101) q[3];
sx q[3];
rz(-0.71037358) q[3];
sx q[3];
rz(-0.47391665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0455857) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(-1.879479) q[2];
rz(-1.1612085) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68814174) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(1.5334817) q[0];
rz(0.43117943) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(1.7327488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7745902) q[0];
sx q[0];
rz(-0.95159114) q[0];
sx q[0];
rz(1.9253233) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3546014) q[2];
sx q[2];
rz(-0.39297418) q[2];
sx q[2];
rz(2.4525688) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81804336) q[1];
sx q[1];
rz(-0.20092873) q[1];
sx q[1];
rz(-0.6769606) q[1];
rz(-pi) q[2];
rz(1.3471782) q[3];
sx q[3];
rz(-2.6287492) q[3];
sx q[3];
rz(-1.4033069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93671736) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(-1.3196866) q[2];
rz(-3.1070869) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(-1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40783229) q[0];
sx q[0];
rz(-2.5875081) q[0];
sx q[0];
rz(-1.2169417) q[0];
rz(-0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(1.9409174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9618098) q[0];
sx q[0];
rz(-2.2968074) q[0];
sx q[0];
rz(2.2711193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1633792) q[2];
sx q[2];
rz(-2.1317185) q[2];
sx q[2];
rz(1.4108059) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4958805) q[1];
sx q[1];
rz(-0.79949841) q[1];
sx q[1];
rz(1.6281566) q[1];
rz(-pi) q[2];
rz(2.7755599) q[3];
sx q[3];
rz(-1.1732374) q[3];
sx q[3];
rz(-1.1432858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.919148) q[2];
sx q[2];
rz(-1.6067952) q[2];
sx q[2];
rz(-2.0929125) q[2];
rz(2.9564296) q[3];
sx q[3];
rz(-1.1566297) q[3];
sx q[3];
rz(-0.10704253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.01934) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(-2.3308603) q[0];
rz(2.4108389) q[1];
sx q[1];
rz(-1.0565051) q[1];
sx q[1];
rz(0.96053851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.273868) q[0];
sx q[0];
rz(-1.5717713) q[0];
sx q[0];
rz(0.0033897059) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4308295) q[2];
sx q[2];
rz(-1.1053876) q[2];
sx q[2];
rz(-0.39164603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5170578) q[1];
sx q[1];
rz(-1.0427999) q[1];
sx q[1];
rz(2.8370884) q[1];
rz(-pi) q[2];
rz(-1.0455564) q[3];
sx q[3];
rz(-2.6572795) q[3];
sx q[3];
rz(-2.7698295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1824128) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(1.5394999) q[2];
rz(-2.0942073) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(-1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-0.86532101) q[0];
sx q[0];
rz(-0.58746946) q[0];
rz(-0.36382183) q[1];
sx q[1];
rz(-1.9963341) q[1];
sx q[1];
rz(-2.2344373) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.912187) q[0];
sx q[0];
rz(-2.7813781) q[0];
sx q[0];
rz(0.2304669) q[0];
rz(-0.67423363) q[2];
sx q[2];
rz(-0.57949726) q[2];
sx q[2];
rz(0.63636875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3876411) q[1];
sx q[1];
rz(-0.50261897) q[1];
sx q[1];
rz(0.44746621) q[1];
rz(-pi) q[2];
rz(-1.3499522) q[3];
sx q[3];
rz(-2.7691602) q[3];
sx q[3];
rz(2.7999634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5017447) q[2];
sx q[2];
rz(-3.0526243) q[2];
sx q[2];
rz(-2.7196344) q[2];
rz(-3.1320599) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137153) q[0];
sx q[0];
rz(-1.2953225) q[0];
sx q[0];
rz(-3.0184826) q[0];
rz(-0.70115024) q[1];
sx q[1];
rz(-1.2260561) q[1];
sx q[1];
rz(-1.4969926) q[1];
rz(2.7098165) q[2];
sx q[2];
rz(-0.74413055) q[2];
sx q[2];
rz(-1.6091138) q[2];
rz(-1.9811859) q[3];
sx q[3];
rz(-1.8262564) q[3];
sx q[3];
rz(-1.5953596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
