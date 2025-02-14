OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8799514) q[0];
sx q[0];
rz(-2.9961442) q[0];
sx q[0];
rz(2.8737336) q[0];
rz(-1.5675867) q[1];
sx q[1];
rz(-0.1685473) q[1];
sx q[1];
rz(0.57810098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.472653) q[0];
sx q[0];
rz(-2.2340206) q[0];
sx q[0];
rz(-0.35314631) q[0];
x q[1];
rz(0.16139754) q[2];
sx q[2];
rz(-0.70490743) q[2];
sx q[2];
rz(0.43625956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4207124) q[1];
sx q[1];
rz(-2.0053177) q[1];
sx q[1];
rz(-2.0220533) q[1];
rz(-pi) q[2];
rz(-1.7655444) q[3];
sx q[3];
rz(-0.71310242) q[3];
sx q[3];
rz(2.3793067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9428923) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(1.6655507) q[2];
rz(2.5623411) q[3];
sx q[3];
rz(-1.9871291) q[3];
sx q[3];
rz(0.66592413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2838374) q[0];
sx q[0];
rz(-2.0960161) q[0];
sx q[0];
rz(-0.055140821) q[0];
rz(-2.227123) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(-1.2158016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037579868) q[0];
sx q[0];
rz(-1.6336339) q[0];
sx q[0];
rz(1.5893905) q[0];
rz(-pi) q[1];
rz(2.397695) q[2];
sx q[2];
rz(-2.4716931) q[2];
sx q[2];
rz(-0.11775859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3926331) q[1];
sx q[1];
rz(-1.5523504) q[1];
sx q[1];
rz(-0.39876826) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8626446) q[3];
sx q[3];
rz(-0.46542811) q[3];
sx q[3];
rz(-0.60681776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9332283) q[2];
sx q[2];
rz(-0.51839447) q[2];
sx q[2];
rz(-0.62057692) q[2];
rz(-1.6879451) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(-2.0645963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8715912) q[0];
sx q[0];
rz(-0.33137614) q[0];
sx q[0];
rz(1.510386) q[0];
rz(-2.467678) q[1];
sx q[1];
rz(-1.2951415) q[1];
sx q[1];
rz(-1.6253701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5249859) q[0];
sx q[0];
rz(-1.5956889) q[0];
sx q[0];
rz(2.1014433) q[0];
x q[1];
rz(-1.0732365) q[2];
sx q[2];
rz(-2.4399827) q[2];
sx q[2];
rz(0.09532433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61805815) q[1];
sx q[1];
rz(-0.48413545) q[1];
sx q[1];
rz(0.50393288) q[1];
x q[2];
rz(0.55530352) q[3];
sx q[3];
rz(-1.318299) q[3];
sx q[3];
rz(0.50336526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7565833) q[2];
sx q[2];
rz(-1.5998806) q[2];
sx q[2];
rz(2.4760683) q[2];
rz(-2.3982128) q[3];
sx q[3];
rz(-1.1210818) q[3];
sx q[3];
rz(-2.9140748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7706364) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(1.3002522) q[0];
rz(-3.0216253) q[1];
sx q[1];
rz(-1.080546) q[1];
sx q[1];
rz(2.0948476) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9366953) q[0];
sx q[0];
rz(-0.012210695) q[0];
sx q[0];
rz(1.3365251) q[0];
x q[1];
rz(-0.8530944) q[2];
sx q[2];
rz(-1.8522369) q[2];
sx q[2];
rz(-2.6996343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6347329) q[1];
sx q[1];
rz(-1.6258937) q[1];
sx q[1];
rz(-1.0841838) q[1];
x q[2];
rz(0.98174121) q[3];
sx q[3];
rz(-1.1906173) q[3];
sx q[3];
rz(1.5693992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3469424) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(2.3700355) q[2];
rz(1.0509342) q[3];
sx q[3];
rz(-1.0335048) q[3];
sx q[3];
rz(-1.1933901) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13450384) q[0];
sx q[0];
rz(-1.6935885) q[0];
sx q[0];
rz(0.80528468) q[0];
rz(-1.6979506) q[1];
sx q[1];
rz(-0.81595683) q[1];
sx q[1];
rz(-3.1051292) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0092468) q[0];
sx q[0];
rz(-1.5224592) q[0];
sx q[0];
rz(-2.1383907) q[0];
rz(-3.0731673) q[2];
sx q[2];
rz(-2.1379316) q[2];
sx q[2];
rz(2.4829645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3111687) q[1];
sx q[1];
rz(-1.0941545) q[1];
sx q[1];
rz(-3.0817506) q[1];
rz(-2.8920435) q[3];
sx q[3];
rz(-0.59246906) q[3];
sx q[3];
rz(-2.1268001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61520758) q[2];
sx q[2];
rz(-2.0168763) q[2];
sx q[2];
rz(-0.57662326) q[2];
rz(1.5455101) q[3];
sx q[3];
rz(-2.2871064) q[3];
sx q[3];
rz(0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1103519) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(0.59979576) q[0];
rz(0.80232969) q[1];
sx q[1];
rz(-1.0917412) q[1];
sx q[1];
rz(-0.64782992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9947056) q[0];
sx q[0];
rz(-1.1297261) q[0];
sx q[0];
rz(-0.85710454) q[0];
rz(-pi) q[1];
rz(3.134583) q[2];
sx q[2];
rz(-1.3607549) q[2];
sx q[2];
rz(1.5633068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9279865) q[1];
sx q[1];
rz(-1.2756375) q[1];
sx q[1];
rz(0.61156433) q[1];
rz(0.69979005) q[3];
sx q[3];
rz(-2.095053) q[3];
sx q[3];
rz(-1.4211224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.99737793) q[2];
sx q[2];
rz(-0.37313676) q[2];
sx q[2];
rz(2.0443661) q[2];
rz(2.1413474) q[3];
sx q[3];
rz(-2.2179243) q[3];
sx q[3];
rz(1.3057115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0344987) q[0];
sx q[0];
rz(-1.1544363) q[0];
sx q[0];
rz(-2.9423998) q[0];
rz(-1.2983324) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(2.4995506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605153) q[0];
sx q[0];
rz(-2.2937728) q[0];
sx q[0];
rz(-2.8882746) q[0];
rz(-pi) q[1];
rz(-1.0429616) q[2];
sx q[2];
rz(-1.9342074) q[2];
sx q[2];
rz(0.052415457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7016546) q[1];
sx q[1];
rz(-1.3716193) q[1];
sx q[1];
rz(-1.4889731) q[1];
rz(2.9391791) q[3];
sx q[3];
rz(-2.0127986) q[3];
sx q[3];
rz(1.3998264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22884998) q[2];
sx q[2];
rz(-1.0656837) q[2];
sx q[2];
rz(-1.9672811) q[2];
rz(2.2187388) q[3];
sx q[3];
rz(-0.43761161) q[3];
sx q[3];
rz(0.16930425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787503) q[0];
sx q[0];
rz(-1.6599382) q[0];
sx q[0];
rz(-2.1424275) q[0];
rz(0.83813465) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(-2.1160486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2472256) q[0];
sx q[0];
rz(-2.0400088) q[0];
sx q[0];
rz(2.2870025) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3633358) q[2];
sx q[2];
rz(-2.7822251) q[2];
sx q[2];
rz(-3.0148413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22502314) q[1];
sx q[1];
rz(-0.65316641) q[1];
sx q[1];
rz(2.2520425) q[1];
rz(-pi) q[2];
rz(2.4225967) q[3];
sx q[3];
rz(-1.1482802) q[3];
sx q[3];
rz(-2.7718294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(-0.21719246) q[2];
rz(-1.1413261) q[3];
sx q[3];
rz(-2.7666028) q[3];
sx q[3];
rz(0.91514897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7104257) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(0.024209484) q[0];
rz(-1.0912033) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(-2.3275163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59210789) q[0];
sx q[0];
rz(-2.4257437) q[0];
sx q[0];
rz(0.13087337) q[0];
rz(2.5957554) q[2];
sx q[2];
rz(-1.0190735) q[2];
sx q[2];
rz(1.9474701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3755889) q[1];
sx q[1];
rz(-0.75889041) q[1];
sx q[1];
rz(1.5726456) q[1];
rz(-pi) q[2];
rz(1.9283251) q[3];
sx q[3];
rz(-0.24328624) q[3];
sx q[3];
rz(-1.4668087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39390627) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(-1.5300592) q[2];
rz(1.0511506) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(-0.10147258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53973389) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-0.47716004) q[0];
rz(-1.5513783) q[1];
sx q[1];
rz(-1.662622) q[1];
sx q[1];
rz(-1.8345376) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64335075) q[0];
sx q[0];
rz(-2.4242867) q[0];
sx q[0];
rz(-1.6397315) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3315694) q[2];
sx q[2];
rz(-0.56213435) q[2];
sx q[2];
rz(-0.13371828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.160678) q[1];
sx q[1];
rz(-1.0458993) q[1];
sx q[1];
rz(0.69634931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1496592) q[3];
sx q[3];
rz(-2.4483878) q[3];
sx q[3];
rz(-2.7424911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16605475) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(2.3892152) q[2];
rz(0.11263975) q[3];
sx q[3];
rz(-0.17387667) q[3];
sx q[3];
rz(1.9788474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024121506) q[0];
sx q[0];
rz(-1.385067) q[0];
sx q[0];
rz(-2.0140482) q[0];
rz(2.7680001) q[1];
sx q[1];
rz(-1.5126546) q[1];
sx q[1];
rz(1.0027813) q[1];
rz(-0.75171555) q[2];
sx q[2];
rz(-2.5645419) q[2];
sx q[2];
rz(0.6286055) q[2];
rz(-1.7711025) q[3];
sx q[3];
rz(-2.3080993) q[3];
sx q[3];
rz(-0.75626683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
