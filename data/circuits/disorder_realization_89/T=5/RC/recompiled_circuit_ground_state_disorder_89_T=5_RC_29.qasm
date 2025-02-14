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
rz(1.574006) q[1];
sx q[1];
rz(-2.9730453) q[1];
sx q[1];
rz(-0.57810098) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0165839) q[0];
sx q[0];
rz(-1.2947695) q[0];
sx q[0];
rz(2.2651423) q[0];
rz(-pi) q[1];
rz(0.69845818) q[2];
sx q[2];
rz(-1.6751117) q[2];
sx q[2];
rz(1.2579174) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0904304) q[1];
sx q[1];
rz(-1.9774861) q[1];
sx q[1];
rz(2.6654408) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9757479) q[3];
sx q[3];
rz(-2.2676617) q[3];
sx q[3];
rz(-2.6344185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9428923) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(-1.6655507) q[2];
rz(-2.5623411) q[3];
sx q[3];
rz(-1.1544635) q[3];
sx q[3];
rz(0.66592413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2838374) q[0];
sx q[0];
rz(-2.0960161) q[0];
sx q[0];
rz(-0.055140821) q[0];
rz(-0.91446963) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(-1.9257911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25032917) q[0];
sx q[0];
rz(-3.0760652) q[0];
sx q[0];
rz(0.28732462) q[0];
x q[1];
rz(-0.74389768) q[2];
sx q[2];
rz(-0.66989952) q[2];
sx q[2];
rz(-3.0238341) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2760073) q[1];
sx q[1];
rz(-0.39917183) q[1];
sx q[1];
rz(0.047476032) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4333822) q[3];
sx q[3];
rz(-2.0169037) q[3];
sx q[3];
rz(-0.29665369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2083644) q[2];
sx q[2];
rz(-0.51839447) q[2];
sx q[2];
rz(-2.5210157) q[2];
rz(1.4536475) q[3];
sx q[3];
rz(-1.0317289) q[3];
sx q[3];
rz(2.0645963) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8715912) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(-1.6312067) q[0];
rz(-0.67391467) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(-1.6253701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5249859) q[0];
sx q[0];
rz(-1.5956889) q[0];
sx q[0];
rz(-1.0401493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2095334) q[2];
sx q[2];
rz(-1.8839508) q[2];
sx q[2];
rz(1.0822288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61805815) q[1];
sx q[1];
rz(-2.6574572) q[1];
sx q[1];
rz(0.50393288) q[1];
rz(-0.45511873) q[3];
sx q[3];
rz(-0.60448217) q[3];
sx q[3];
rz(0.68463078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38500938) q[2];
sx q[2];
rz(-1.5998806) q[2];
sx q[2];
rz(0.6655244) q[2];
rz(-0.74337983) q[3];
sx q[3];
rz(-1.1210818) q[3];
sx q[3];
rz(2.9140748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7706364) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(1.8413405) q[0];
rz(3.0216253) q[1];
sx q[1];
rz(-2.0610466) q[1];
sx q[1];
rz(2.0948476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20489731) q[0];
sx q[0];
rz(-0.012210695) q[0];
sx q[0];
rz(-1.3365251) q[0];
x q[1];
rz(0.8530944) q[2];
sx q[2];
rz(-1.2893558) q[2];
sx q[2];
rz(0.44195833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0485301) q[1];
sx q[1];
rz(-1.0849864) q[1];
sx q[1];
rz(0.062315224) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94724058) q[3];
sx q[3];
rz(-0.6886607) q[3];
sx q[3];
rz(2.6332444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3469424) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(2.3700355) q[2];
rz(1.0509342) q[3];
sx q[3];
rz(-2.1080878) q[3];
sx q[3];
rz(1.1933901) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59234658) q[0];
sx q[0];
rz(-1.0039465) q[0];
sx q[0];
rz(-3.0842848) q[0];
rz(1.4638605) q[2];
sx q[2];
rz(-0.57079878) q[2];
sx q[2];
rz(-2.609848) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.441022) q[1];
sx q[1];
rz(-0.4800969) q[1];
sx q[1];
rz(-1.455485) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4060601) q[3];
sx q[3];
rz(-2.142557) q[3];
sx q[3];
rz(-0.71672196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61520758) q[2];
sx q[2];
rz(-1.1247164) q[2];
sx q[2];
rz(2.5649694) q[2];
rz(1.5960826) q[3];
sx q[3];
rz(-2.2871064) q[3];
sx q[3];
rz(-0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1103519) q[0];
sx q[0];
rz(-1.6540271) q[0];
sx q[0];
rz(-0.59979576) q[0];
rz(-0.80232969) q[1];
sx q[1];
rz(-2.0498514) q[1];
sx q[1];
rz(-0.64782992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9947056) q[0];
sx q[0];
rz(-2.0118666) q[0];
sx q[0];
rz(-0.85710454) q[0];
x q[1];
rz(1.3607499) q[2];
sx q[2];
rz(-1.5639407) q[2];
sx q[2];
rz(-3.1355646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55841724) q[1];
sx q[1];
rz(-2.1523471) q[1];
sx q[1];
rz(1.9263595) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4418026) q[3];
sx q[3];
rz(-1.0465396) q[3];
sx q[3];
rz(-1.4211224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1442147) q[2];
sx q[2];
rz(-2.7684559) q[2];
sx q[2];
rz(-2.0443661) q[2];
rz(2.1413474) q[3];
sx q[3];
rz(-0.92366832) q[3];
sx q[3];
rz(-1.3057115) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0344987) q[0];
sx q[0];
rz(-1.1544363) q[0];
sx q[0];
rz(2.9423998) q[0];
rz(1.2983324) q[1];
sx q[1];
rz(-2.7804764) q[1];
sx q[1];
rz(-0.64204204) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605153) q[0];
sx q[0];
rz(-2.2937728) q[0];
sx q[0];
rz(-0.25331803) q[0];
rz(2.2175353) q[2];
sx q[2];
rz(-2.5106608) q[2];
sx q[2];
rz(1.0756878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7016546) q[1];
sx q[1];
rz(-1.3716193) q[1];
sx q[1];
rz(-1.4889731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9391791) q[3];
sx q[3];
rz(-2.0127986) q[3];
sx q[3];
rz(1.3998264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9127427) q[2];
sx q[2];
rz(-1.0656837) q[2];
sx q[2];
rz(-1.1743116) q[2];
rz(2.2187388) q[3];
sx q[3];
rz(-0.43761161) q[3];
sx q[3];
rz(-2.9722884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6628424) q[0];
sx q[0];
rz(-1.4816544) q[0];
sx q[0];
rz(0.99916512) q[0];
rz(-2.303458) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(1.025544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89436707) q[0];
sx q[0];
rz(-1.1015838) q[0];
sx q[0];
rz(-2.2870025) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26141459) q[2];
sx q[2];
rz(-1.3213219) q[2];
sx q[2];
rz(0.69875137) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22502314) q[1];
sx q[1];
rz(-0.65316641) q[1];
sx q[1];
rz(0.8895501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5426264) q[3];
sx q[3];
rz(-2.3271797) q[3];
sx q[3];
rz(1.502445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.7842224) q[2];
sx q[2];
rz(0.21719246) q[2];
rz(-1.1413261) q[3];
sx q[3];
rz(-2.7666028) q[3];
sx q[3];
rz(0.91514897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4311669) q[0];
sx q[0];
rz(-1.8473666) q[0];
sx q[0];
rz(-0.024209484) q[0];
rz(1.0912033) q[1];
sx q[1];
rz(-2.0228491) q[1];
sx q[1];
rz(-2.3275163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494848) q[0];
sx q[0];
rz(-2.4257437) q[0];
sx q[0];
rz(-0.13087337) q[0];
rz(2.5957554) q[2];
sx q[2];
rz(-2.1225192) q[2];
sx q[2];
rz(1.1941225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1938656) q[1];
sx q[1];
rz(-1.5720688) q[1];
sx q[1];
rz(-0.81190677) q[1];
rz(-pi) q[2];
rz(-3.0549491) q[3];
sx q[3];
rz(-1.7984172) q[3];
sx q[3];
rz(1.3073453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7476864) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(-1.6115335) q[2];
rz(-2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(3.0401201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6018588) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-2.6644326) q[0];
rz(1.5513783) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(-1.8345376) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64335075) q[0];
sx q[0];
rz(-0.71730591) q[0];
sx q[0];
rz(1.5018612) q[0];
rz(-pi) q[1];
x q[1];
rz(2.731852) q[2];
sx q[2];
rz(-1.9671408) q[2];
sx q[2];
rz(0.71045638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8068731) q[1];
sx q[1];
rz(-2.1592437) q[1];
sx q[1];
rz(-0.92422331) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1496592) q[3];
sx q[3];
rz(-2.4483878) q[3];
sx q[3];
rz(-0.39910152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9755379) q[2];
sx q[2];
rz(-1.0155948) q[2];
sx q[2];
rz(2.3892152) q[2];
rz(3.0289529) q[3];
sx q[3];
rz(-2.967716) q[3];
sx q[3];
rz(1.9788474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024121506) q[0];
sx q[0];
rz(-1.385067) q[0];
sx q[0];
rz(-2.0140482) q[0];
rz(0.37359259) q[1];
sx q[1];
rz(-1.6289381) q[1];
sx q[1];
rz(-2.1388114) q[1];
rz(2.6977044) q[2];
sx q[2];
rz(-1.9525575) q[2];
sx q[2];
rz(-1.606751) q[2];
rz(0.21568372) q[3];
sx q[3];
rz(-2.3825323) q[3];
sx q[3];
rz(2.0920405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
