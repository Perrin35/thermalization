OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2616413) q[0];
sx q[0];
rz(-0.14544848) q[0];
sx q[0];
rz(0.26785904) q[0];
rz(-1.5675867) q[1];
sx q[1];
rz(-0.1685473) q[1];
sx q[1];
rz(0.57810098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0165839) q[0];
sx q[0];
rz(-1.8468231) q[0];
sx q[0];
rz(0.87645032) q[0];
rz(-pi) q[1];
rz(-0.16139754) q[2];
sx q[2];
rz(-2.4366852) q[2];
sx q[2];
rz(0.43625956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0511623) q[1];
sx q[1];
rz(-1.9774861) q[1];
sx q[1];
rz(-0.47615188) q[1];
x q[2];
rz(1.3760482) q[3];
sx q[3];
rz(-0.71310242) q[3];
sx q[3];
rz(-0.76228599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1987004) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(1.4760419) q[2];
rz(2.5623411) q[3];
sx q[3];
rz(-1.1544635) q[3];
sx q[3];
rz(-0.66592413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85775527) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(-3.0864518) q[0];
rz(-2.227123) q[1];
sx q[1];
rz(-1.4417459) q[1];
sx q[1];
rz(1.2158016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5343842) q[0];
sx q[0];
rz(-1.5522389) q[0];
sx q[0];
rz(0.062848363) q[0];
rz(-pi) q[1];
rz(-2.0631172) q[2];
sx q[2];
rz(-2.0452867) q[2];
sx q[2];
rz(-2.1585787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3926331) q[1];
sx q[1];
rz(-1.5892423) q[1];
sx q[1];
rz(0.39876826) q[1];
rz(-1.7082105) q[3];
sx q[3];
rz(-1.124689) q[3];
sx q[3];
rz(-2.844939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2083644) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(-0.62057692) q[2];
rz(-1.4536475) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(-1.0769963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8715912) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(-1.6312067) q[0];
rz(0.67391467) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(1.6253701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103863) q[0];
sx q[0];
rz(-1.0403311) q[0];
sx q[0];
rz(3.112733) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0683561) q[2];
sx q[2];
rz(-0.70161) q[2];
sx q[2];
rz(0.09532433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5235345) q[1];
sx q[1];
rz(-2.6574572) q[1];
sx q[1];
rz(0.50393288) q[1];
rz(-0.55530352) q[3];
sx q[3];
rz(-1.318299) q[3];
sx q[3];
rz(2.6382274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38500938) q[2];
sx q[2];
rz(-1.5417121) q[2];
sx q[2];
rz(-2.4760683) q[2];
rz(2.3982128) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(-2.9140748) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7706364) q[0];
sx q[0];
rz(-2.0891068) q[0];
sx q[0];
rz(-1.8413405) q[0];
rz(3.0216253) q[1];
sx q[1];
rz(-2.0610466) q[1];
sx q[1];
rz(2.0948476) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7024073) q[0];
sx q[0];
rz(-1.5826735) q[0];
sx q[0];
rz(-0.0028346527) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.98498) q[2];
sx q[2];
rz(-0.76167548) q[2];
sx q[2];
rz(1.7047395) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1018104) q[1];
sx q[1];
rz(-0.48947316) q[1];
sx q[1];
rz(1.4533978) q[1];
rz(-pi) q[2];
rz(-2.1598514) q[3];
sx q[3];
rz(-1.1906173) q[3];
sx q[3];
rz(1.5693992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3469424) q[2];
sx q[2];
rz(-1.0618989) q[2];
sx q[2];
rz(2.3700355) q[2];
rz(1.0509342) q[3];
sx q[3];
rz(-1.0335048) q[3];
sx q[3];
rz(1.9482025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0070888) q[0];
sx q[0];
rz(-1.4480042) q[0];
sx q[0];
rz(-2.336308) q[0];
rz(-1.6979506) q[1];
sx q[1];
rz(-0.81595683) q[1];
sx q[1];
rz(0.03646341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48591025) q[0];
sx q[0];
rz(-2.5721689) q[0];
sx q[0];
rz(-1.481056) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6777322) q[2];
sx q[2];
rz(-2.5707939) q[2];
sx q[2];
rz(2.609848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9094463) q[1];
sx q[1];
rz(-1.6239617) q[1];
sx q[1];
rz(1.0934238) q[1];
rz(-pi) q[2];
rz(1.7355326) q[3];
sx q[3];
rz(-2.142557) q[3];
sx q[3];
rz(-2.4248707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1103519) q[0];
sx q[0];
rz(-1.6540271) q[0];
sx q[0];
rz(2.5417969) q[0];
rz(2.339263) q[1];
sx q[1];
rz(-1.0917412) q[1];
sx q[1];
rz(-2.4937627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14688705) q[0];
sx q[0];
rz(-1.1297261) q[0];
sx q[0];
rz(2.2844881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3607499) q[2];
sx q[2];
rz(-1.5639407) q[2];
sx q[2];
rz(0.0060280212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5831754) q[1];
sx q[1];
rz(-0.98924556) q[1];
sx q[1];
rz(1.9263595) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69979005) q[3];
sx q[3];
rz(-2.095053) q[3];
sx q[3];
rz(1.7204703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1442147) q[2];
sx q[2];
rz(-2.7684559) q[2];
sx q[2];
rz(1.0972265) q[2];
rz(-1.0002452) q[3];
sx q[3];
rz(-2.2179243) q[3];
sx q[3];
rz(1.3057115) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(-0.19919285) q[0];
rz(1.8432603) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(-0.64204204) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10810971) q[0];
sx q[0];
rz(-2.3831522) q[0];
sx q[0];
rz(1.8475425) q[0];
rz(-pi) q[1];
rz(0.92405739) q[2];
sx q[2];
rz(-0.63093189) q[2];
sx q[2];
rz(-2.0659049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7016546) q[1];
sx q[1];
rz(-1.3716193) q[1];
sx q[1];
rz(1.6526196) q[1];
rz(-pi) q[2];
rz(1.1207668) q[3];
sx q[3];
rz(-1.3880669) q[3];
sx q[3];
rz(-0.25853005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9127427) q[2];
sx q[2];
rz(-2.0759089) q[2];
sx q[2];
rz(-1.1743116) q[2];
rz(2.2187388) q[3];
sx q[3];
rz(-2.703981) q[3];
sx q[3];
rz(2.9722884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787503) q[0];
sx q[0];
rz(-1.6599382) q[0];
sx q[0];
rz(-2.1424275) q[0];
rz(-2.303458) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(-2.1160486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8401153) q[0];
sx q[0];
rz(-0.94520926) q[0];
sx q[0];
rz(0.59176318) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8286669) q[2];
sx q[2];
rz(-1.317655) q[2];
sx q[2];
rz(-0.93800256) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9165695) q[1];
sx q[1];
rz(-0.65316641) q[1];
sx q[1];
rz(0.8895501) q[1];
rz(-pi) q[2];
rz(2.4225967) q[3];
sx q[3];
rz(-1.9933125) q[3];
sx q[3];
rz(2.7718294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(2.9244002) q[2];
rz(-2.0002666) q[3];
sx q[3];
rz(-2.7666028) q[3];
sx q[3];
rz(2.2264437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4311669) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-0.024209484) q[0];
rz(-1.0912033) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(-2.3275163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0639125) q[0];
sx q[0];
rz(-1.4850495) q[0];
sx q[0];
rz(0.7115988) q[0];
x q[1];
rz(-2.2714839) q[2];
sx q[2];
rz(-0.75557709) q[2];
sx q[2];
rz(0.33529624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7660038) q[1];
sx q[1];
rz(-0.75889041) q[1];
sx q[1];
rz(1.5726456) q[1];
rz(-pi) q[2];
rz(-0.086643593) q[3];
sx q[3];
rz(-1.3431755) q[3];
sx q[3];
rz(-1.8342474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7476864) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(1.5300592) q[2];
rz(-2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(3.0401201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6018588) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-0.47716004) q[0];
rz(1.5513783) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(1.307055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64335075) q[0];
sx q[0];
rz(-2.4242867) q[0];
sx q[0];
rz(1.5018612) q[0];
rz(0.81002323) q[2];
sx q[2];
rz(-2.5794583) q[2];
sx q[2];
rz(0.13371828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3347196) q[1];
sx q[1];
rz(-0.98234896) q[1];
sx q[1];
rz(-2.2173693) q[1];
rz(-pi) q[2];
rz(2.1784276) q[3];
sx q[3];
rz(-1.9279216) q[3];
sx q[3];
rz(2.4357093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16605475) q[2];
sx q[2];
rz(-1.0155948) q[2];
sx q[2];
rz(0.75237742) q[2];
rz(-3.0289529) q[3];
sx q[3];
rz(-2.967716) q[3];
sx q[3];
rz(1.1627452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(0.37359259) q[1];
sx q[1];
rz(-1.6289381) q[1];
sx q[1];
rz(-2.1388114) q[1];
rz(0.75171555) q[2];
sx q[2];
rz(-0.57705078) q[2];
sx q[2];
rz(-2.5129872) q[2];
rz(1.7711025) q[3];
sx q[3];
rz(-0.8334934) q[3];
sx q[3];
rz(2.3853258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
