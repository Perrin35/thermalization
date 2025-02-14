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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(1.8408096) q[1];
sx q[1];
rz(9.169133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545171) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(-2.4357027) q[0];
x q[1];
rz(-2.7529703) q[2];
sx q[2];
rz(-1.6291233) q[2];
sx q[2];
rz(-1.3921392) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5140443) q[1];
sx q[1];
rz(-1.5020292) q[1];
sx q[1];
rz(1.7538944) q[1];
rz(-pi) q[2];
rz(-2.643804) q[3];
sx q[3];
rz(-1.3827795) q[3];
sx q[3];
rz(1.9888442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(2.0261436) q[2];
rz(-3.1274146) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(-2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071335) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(1.8943262) q[0];
rz(-2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(2.3242548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8431664) q[0];
sx q[0];
rz(-2.0930556) q[0];
sx q[0];
rz(-2.9199615) q[0];
rz(-pi) q[1];
rz(2.8645682) q[2];
sx q[2];
rz(-1.1454586) q[2];
sx q[2];
rz(-2.9556731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0603588) q[1];
sx q[1];
rz(-2.4801755) q[1];
sx q[1];
rz(-2.8702186) q[1];
rz(0.67561356) q[3];
sx q[3];
rz(-2.2573376) q[3];
sx q[3];
rz(0.10304777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(-2.9212941) q[2];
rz(1.1822654) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050506266) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(2.128111) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0419755) q[0];
sx q[0];
rz(-1.8238153) q[0];
sx q[0];
rz(-1.5861041) q[0];
x q[1];
rz(-3.0092952) q[2];
sx q[2];
rz(-2.0084642) q[2];
sx q[2];
rz(-0.79298151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1057844) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(-0.99143274) q[1];
rz(2.2974422) q[3];
sx q[3];
rz(-1.4733847) q[3];
sx q[3];
rz(-2.0908434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1608405) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(-2.8694425) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(-3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18144064) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(-2.5312359) q[0];
rz(0.70862526) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(-1.5708539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7371668) q[0];
sx q[0];
rz(-0.20658399) q[0];
sx q[0];
rz(-2.4086359) q[0];
x q[1];
rz(-1.6897292) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(2.1708084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4492693) q[1];
sx q[1];
rz(-0.27579112) q[1];
sx q[1];
rz(0.17755601) q[1];
x q[2];
rz(1.1669772) q[3];
sx q[3];
rz(-0.81203264) q[3];
sx q[3];
rz(-0.62893516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(0.55541682) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(-0.24636191) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(-1.7074283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96077641) q[0];
sx q[0];
rz(-2.3276276) q[0];
sx q[0];
rz(1.5636982) q[0];
rz(-pi) q[1];
rz(1.3952012) q[2];
sx q[2];
rz(-1.1313952) q[2];
sx q[2];
rz(1.1281769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38404462) q[1];
sx q[1];
rz(-1.169012) q[1];
sx q[1];
rz(-1.5753788) q[1];
rz(-1.7982539) q[3];
sx q[3];
rz(-1.47336) q[3];
sx q[3];
rz(0.9909329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(2.0595713) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(1.2580416) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(0.24359447) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(-2.2056244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0662066) q[0];
sx q[0];
rz(-0.52156007) q[0];
sx q[0];
rz(-0.88514502) q[0];
rz(-pi) q[1];
rz(-0.82359969) q[2];
sx q[2];
rz(-1.8317501) q[2];
sx q[2];
rz(-1.5023155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0985581) q[1];
sx q[1];
rz(-1.9717311) q[1];
sx q[1];
rz(2.3844196) q[1];
x q[2];
rz(1.5395079) q[3];
sx q[3];
rz(-0.60006053) q[3];
sx q[3];
rz(-2.4158583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4794856) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(-2.7916419) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(-0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(-2.8217577) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(-2.005827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58590305) q[0];
sx q[0];
rz(-1.5169889) q[0];
sx q[0];
rz(2.4967628) q[0];
rz(-2.2280424) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(-0.28829703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37138501) q[1];
sx q[1];
rz(-1.4911629) q[1];
sx q[1];
rz(-0.97038986) q[1];
rz(-pi) q[2];
rz(1.4581751) q[3];
sx q[3];
rz(-2.3584692) q[3];
sx q[3];
rz(-0.52971958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(2.9988334) q[2];
rz(1.3879294) q[3];
sx q[3];
rz(-1.1889435) q[3];
sx q[3];
rz(-0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.9614354) q[0];
sx q[0];
rz(-3.1249983) q[0];
sx q[0];
rz(0.55602443) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(1.3831327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616104) q[0];
sx q[0];
rz(-1.5357657) q[0];
sx q[0];
rz(0.73765124) q[0];
rz(-1.6298619) q[2];
sx q[2];
rz(-1.1129654) q[2];
sx q[2];
rz(0.77564592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(2.8252557) q[1];
rz(-1.9541627) q[3];
sx q[3];
rz(-1.995242) q[3];
sx q[3];
rz(-2.2147873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(1.9005091) q[2];
rz(-0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0144219) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(-2.993809) q[0];
rz(-0.5087018) q[1];
sx q[1];
rz(-1.3721507) q[1];
sx q[1];
rz(0.37857372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258647) q[0];
sx q[0];
rz(-1.3368538) q[0];
sx q[0];
rz(-2.7153003) q[0];
x q[1];
rz(0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(2.1155807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55146688) q[1];
sx q[1];
rz(-2.6762137) q[1];
sx q[1];
rz(-0.82222934) q[1];
rz(3.0112565) q[3];
sx q[3];
rz(-2.3763083) q[3];
sx q[3];
rz(-1.9627067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(1.4281979) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(-2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7487504) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(3.0774935) q[0];
rz(-2.6129258) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(-2.166523) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1695568) q[0];
sx q[0];
rz(-2.0590879) q[0];
sx q[0];
rz(-1.2704865) q[0];
x q[1];
rz(-0.56985241) q[2];
sx q[2];
rz(-2.0790711) q[2];
sx q[2];
rz(-3.0076671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84555039) q[1];
sx q[1];
rz(-0.49440372) q[1];
sx q[1];
rz(0.47459666) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75210877) q[3];
sx q[3];
rz(-2.8127713) q[3];
sx q[3];
rz(-3.0260835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14572445) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(-2.8912344) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470227) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(2.6429214) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-2.278419) q[2];
sx q[2];
rz(-1.0959846) q[2];
sx q[2];
rz(-1.6352996) q[2];
rz(0.12228431) q[3];
sx q[3];
rz(-1.3663843) q[3];
sx q[3];
rz(-2.1257675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
