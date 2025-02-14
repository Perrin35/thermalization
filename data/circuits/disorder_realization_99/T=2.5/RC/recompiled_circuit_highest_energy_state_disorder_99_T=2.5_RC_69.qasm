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
rz(1.9361629) q[0];
sx q[0];
rz(3.087145) q[0];
sx q[0];
rz(8.5589391) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6893495) q[0];
sx q[0];
rz(-1.6177243) q[0];
sx q[0];
rz(-1.4062792) q[0];
rz(-pi) q[1];
rz(0.073926386) q[2];
sx q[2];
rz(-2.0652899) q[2];
sx q[2];
rz(-1.6520776) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4062351) q[1];
sx q[1];
rz(-2.3182097) q[1];
sx q[1];
rz(-0.75852345) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8179161) q[3];
sx q[3];
rz(-1.2403245) q[3];
sx q[3];
rz(0.068745384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3129348) q[2];
sx q[2];
rz(-2.2609495) q[2];
sx q[2];
rz(-0.9642967) q[2];
rz(3.0621081) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(-0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013414772) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(-1.6012023) q[0];
rz(-0.84114289) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(-0.85743633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6591561) q[0];
sx q[0];
rz(-1.1322339) q[0];
sx q[0];
rz(-0.0060122251) q[0];
rz(2.6688451) q[2];
sx q[2];
rz(-2.1523211) q[2];
sx q[2];
rz(0.66752226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5799478) q[1];
sx q[1];
rz(-2.0873293) q[1];
sx q[1];
rz(-1.0622611) q[1];
rz(-2.8624865) q[3];
sx q[3];
rz(-1.6630904) q[3];
sx q[3];
rz(-1.5444322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88605827) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(1.9319755) q[2];
rz(-1.7602734) q[3];
sx q[3];
rz(-1.2608903) q[3];
sx q[3];
rz(2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7968314) q[0];
sx q[0];
rz(-1.826257) q[0];
sx q[0];
rz(-1.8094081) q[0];
rz(1.4765129) q[1];
sx q[1];
rz(-1.3308728) q[1];
sx q[1];
rz(2.5217893) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22794321) q[0];
sx q[0];
rz(-1.648128) q[0];
sx q[0];
rz(1.5864775) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8263426) q[2];
sx q[2];
rz(-1.4076621) q[2];
sx q[2];
rz(0.28870764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42592749) q[1];
sx q[1];
rz(-1.8879379) q[1];
sx q[1];
rz(1.6460544) q[1];
x q[2];
rz(2.8768646) q[3];
sx q[3];
rz(-0.94784289) q[3];
sx q[3];
rz(0.61443751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0009813112) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(-0.78835431) q[2];
rz(1.2761448) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(0.86281323) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901504) q[0];
sx q[0];
rz(-1.7553512) q[0];
sx q[0];
rz(1.066712) q[0];
rz(1.7881296) q[1];
sx q[1];
rz(-2.2746494) q[1];
sx q[1];
rz(1.1563168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24329127) q[0];
sx q[0];
rz(-1.4292846) q[0];
sx q[0];
rz(-3.1226052) q[0];
x q[1];
rz(2.1422191) q[2];
sx q[2];
rz(-0.048437645) q[2];
sx q[2];
rz(-2.9349309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20322415) q[1];
sx q[1];
rz(-2.6611106) q[1];
sx q[1];
rz(-0.90684569) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66439028) q[3];
sx q[3];
rz(-1.5186465) q[3];
sx q[3];
rz(0.96168226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(3.0827674) q[2];
rz(2.5022653) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-2.2842469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0462129) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(-1.5785716) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(2.6522327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45324126) q[0];
sx q[0];
rz(-2.4972389) q[0];
sx q[0];
rz(2.8768566) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0496558) q[2];
sx q[2];
rz(-2.3917824) q[2];
sx q[2];
rz(-2.9887226) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.64397821) q[1];
sx q[1];
rz(-2.3757739) q[1];
sx q[1];
rz(-0.52180565) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030767083) q[3];
sx q[3];
rz(-1.5813785) q[3];
sx q[3];
rz(2.2969674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8420777) q[2];
sx q[2];
rz(-2.044951) q[2];
sx q[2];
rz(3.1411324) q[2];
rz(0.67356235) q[3];
sx q[3];
rz(-0.898415) q[3];
sx q[3];
rz(0.7091929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29086581) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(-1.3036183) q[0];
rz(-0.29779008) q[1];
sx q[1];
rz(-2.2460263) q[1];
sx q[1];
rz(-0.29388014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760982) q[0];
sx q[0];
rz(-2.2950603) q[0];
sx q[0];
rz(0.80418555) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0877327) q[2];
sx q[2];
rz(-2.392025) q[2];
sx q[2];
rz(2.7550239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9781295) q[1];
sx q[1];
rz(-2.4678964) q[1];
sx q[1];
rz(0.81836318) q[1];
rz(-2.8260293) q[3];
sx q[3];
rz(-1.7949008) q[3];
sx q[3];
rz(-1.699228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6728354) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(0.22077665) q[2];
rz(-1.2729493) q[3];
sx q[3];
rz(-1.7468529) q[3];
sx q[3];
rz(-1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64977589) q[0];
sx q[0];
rz(-0.62115541) q[0];
sx q[0];
rz(-2.5671) q[0];
rz(2.5993787) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(2.7508459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339392) q[0];
sx q[0];
rz(-1.4263194) q[0];
sx q[0];
rz(2.5424315) q[0];
rz(2.1152045) q[2];
sx q[2];
rz(-1.5261391) q[2];
sx q[2];
rz(0.80151973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.407651) q[1];
sx q[1];
rz(-1.9875257) q[1];
sx q[1];
rz(1.5768361) q[1];
rz(-0.096116738) q[3];
sx q[3];
rz(-0.36727723) q[3];
sx q[3];
rz(2.8678081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6812402) q[2];
sx q[2];
rz(-1.9443024) q[2];
sx q[2];
rz(2.528842) q[2];
rz(-0.98226205) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14923444) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(0.055334844) q[0];
rz(1.5845567) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(-2.7016644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377559) q[0];
sx q[0];
rz(-1.3246264) q[0];
sx q[0];
rz(-0.47292821) q[0];
rz(-pi) q[1];
rz(1.3662759) q[2];
sx q[2];
rz(-2.0530982) q[2];
sx q[2];
rz(-0.67687809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5422053) q[1];
sx q[1];
rz(-2.4894307) q[1];
sx q[1];
rz(2.60531) q[1];
x q[2];
rz(-1.1433268) q[3];
sx q[3];
rz(-2.5317268) q[3];
sx q[3];
rz(-2.0855869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7679193) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(0.31201735) q[2];
rz(2.2219374) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298518) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(-3.1134636) q[0];
rz(1.4460538) q[1];
sx q[1];
rz(-1.9606934) q[1];
sx q[1];
rz(-2.6820954) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7336373) q[0];
sx q[0];
rz(-1.6064294) q[0];
sx q[0];
rz(0.028546988) q[0];
x q[1];
rz(0.33513432) q[2];
sx q[2];
rz(-2.2907933) q[2];
sx q[2];
rz(-2.013226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0952009) q[1];
sx q[1];
rz(-1.715379) q[1];
sx q[1];
rz(0.88674366) q[1];
x q[2];
rz(1.2351843) q[3];
sx q[3];
rz(-2.7392355) q[3];
sx q[3];
rz(-1.622792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0405937) q[2];
sx q[2];
rz(-1.3446292) q[2];
sx q[2];
rz(-2.8864268) q[2];
rz(-2.1549639) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(-1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.7607018) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(-1.3421407) q[0];
rz(-3.1396719) q[1];
sx q[1];
rz(-1.0425967) q[1];
sx q[1];
rz(-2.0916746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494019) q[0];
sx q[0];
rz(-1.4225385) q[0];
sx q[0];
rz(-1.3561983) q[0];
rz(-pi) q[1];
rz(1.3948314) q[2];
sx q[2];
rz(-1.5719885) q[2];
sx q[2];
rz(-2.3228485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.695076) q[1];
sx q[1];
rz(-0.43276603) q[1];
sx q[1];
rz(-2.5997735) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9258419) q[3];
sx q[3];
rz(-0.85236824) q[3];
sx q[3];
rz(1.6529447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2125825) q[2];
sx q[2];
rz(-1.9276103) q[2];
sx q[2];
rz(0.49599656) q[2];
rz(-1.6811194) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(-1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37954189) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(2.5024391) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(1.8201309) q[2];
sx q[2];
rz(-1.1087266) q[2];
sx q[2];
rz(-2.5880819) q[2];
rz(0.88981723) q[3];
sx q[3];
rz(-0.87796904) q[3];
sx q[3];
rz(1.8129391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
