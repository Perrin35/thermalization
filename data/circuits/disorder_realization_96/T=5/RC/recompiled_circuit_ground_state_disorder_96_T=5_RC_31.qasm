OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(2.252993) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(-1.3808274) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.439405) q[0];
sx q[0];
rz(-2.8719465) q[0];
sx q[0];
rz(-0.7650956) q[0];
x q[1];
rz(-0.79973296) q[2];
sx q[2];
rz(-1.9603443) q[2];
sx q[2];
rz(-0.79546292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0128768) q[1];
sx q[1];
rz(-2.8059462) q[1];
sx q[1];
rz(0.0291834) q[1];
x q[2];
rz(1.2065776) q[3];
sx q[3];
rz(-1.3602339) q[3];
sx q[3];
rz(2.5488146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-2.3925609) q[2];
rz(-2.9876515) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(-3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.446949) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(1.9248167) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(0.42713508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0674853) q[0];
sx q[0];
rz(-0.68746131) q[0];
sx q[0];
rz(0.56337728) q[0];
x q[1];
rz(2.2160276) q[2];
sx q[2];
rz(-1.7665909) q[2];
sx q[2];
rz(-2.2044971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6500351) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(-0.14658714) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.322836) q[3];
sx q[3];
rz(-2.9313847) q[3];
sx q[3];
rz(1.601364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0624258) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.3986577) q[2];
rz(2.9178197) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021521213) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(3.0960826) q[0];
rz(1.9728164) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(-0.57560903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110014) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(-2.5024947) q[0];
rz(-pi) q[1];
rz(-1.5928629) q[2];
sx q[2];
rz(-0.66078545) q[2];
sx q[2];
rz(1.2118076) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3176885) q[1];
sx q[1];
rz(-1.0270734) q[1];
sx q[1];
rz(-0.25441092) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1166385) q[3];
sx q[3];
rz(-2.2234699) q[3];
sx q[3];
rz(-3.0252473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(0.95727813) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(-1.8804469) q[0];
rz(-3.0592697) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(1.3538768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.533723) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(-1.813867) q[0];
x q[1];
rz(2.655833) q[2];
sx q[2];
rz(-1.4076715) q[2];
sx q[2];
rz(0.32147929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30198797) q[1];
sx q[1];
rz(-2.0652373) q[1];
sx q[1];
rz(1.5093263) q[1];
x q[2];
rz(2.1534377) q[3];
sx q[3];
rz(-0.76863063) q[3];
sx q[3];
rz(2.5386794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-0.91840863) q[2];
sx q[2];
rz(-1.8850231) q[2];
rz(-2.4300857) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-2.7984483) q[0];
rz(-0.061773069) q[1];
sx q[1];
rz(-2.1715178) q[1];
sx q[1];
rz(1.4168581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34300229) q[0];
sx q[0];
rz(-1.8483642) q[0];
sx q[0];
rz(2.9647602) q[0];
rz(2.4409962) q[2];
sx q[2];
rz(-1.3178409) q[2];
sx q[2];
rz(2.7299936) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8843536) q[1];
sx q[1];
rz(-1.7063008) q[1];
sx q[1];
rz(2.0550904) q[1];
x q[2];
rz(-2.9067578) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(-1.0339586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(2.055638) q[2];
rz(-2.2222399) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(-2.5277188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(2.3486775) q[0];
rz(0.74053699) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-2.3103796) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342469) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(-2.9177279) q[0];
rz(-pi) q[1];
rz(-1.5969876) q[2];
sx q[2];
rz(-0.27715836) q[2];
sx q[2];
rz(1.3326534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38922849) q[1];
sx q[1];
rz(-0.78654754) q[1];
sx q[1];
rz(2.0081372) q[1];
rz(2.0415067) q[3];
sx q[3];
rz(-1.5892972) q[3];
sx q[3];
rz(-2.5262566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.071216019) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(2.0224723) q[2];
rz(-1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(-1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374461) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(-1.5420089) q[0];
rz(-1.0150602) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-2.9687845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6369404) q[0];
sx q[0];
rz(-1.0811792) q[0];
sx q[0];
rz(-2.2235653) q[0];
rz(-pi) q[1];
rz(1.9394933) q[2];
sx q[2];
rz(-2.2208344) q[2];
sx q[2];
rz(-1.7391313) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0513251) q[1];
sx q[1];
rz(-0.50053144) q[1];
sx q[1];
rz(2.1085018) q[1];
x q[2];
rz(0.5676078) q[3];
sx q[3];
rz(-2.2877573) q[3];
sx q[3];
rz(-2.6102118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.479594) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.9366555) q[3];
sx q[3];
rz(-1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(0.015722474) q[0];
rz(1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(1.1994919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2000113) q[0];
sx q[0];
rz(-1.1304378) q[0];
sx q[0];
rz(2.8811127) q[0];
rz(0.052414465) q[2];
sx q[2];
rz(-2.0611066) q[2];
sx q[2];
rz(-0.744396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3856628) q[1];
sx q[1];
rz(-2.0634868) q[1];
sx q[1];
rz(-1.0316959) q[1];
x q[2];
rz(-1.1114798) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1477995) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(0.6238474) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(1.0998211) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(-2.4483335) q[0];
rz(-0.26501003) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.6185435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023848195) q[0];
sx q[0];
rz(-1.5818412) q[0];
sx q[0];
rz(0.77356972) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9829282) q[2];
sx q[2];
rz(-1.0318021) q[2];
sx q[2];
rz(1.3608152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0654782) q[1];
sx q[1];
rz(-1.5468742) q[1];
sx q[1];
rz(-0.33965276) q[1];
rz(0.98966815) q[3];
sx q[3];
rz(-1.6299986) q[3];
sx q[3];
rz(-0.35602415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8388464) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(-1.4576853) q[2];
rz(-0.95528209) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569698) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(-0.48267522) q[0];
rz(-2.2118498) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(-2.7453056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.020655) q[0];
sx q[0];
rz(-1.2115941) q[0];
sx q[0];
rz(0.44393702) q[0];
rz(-2.9626289) q[2];
sx q[2];
rz(-1.4216627) q[2];
sx q[2];
rz(2.3863132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.86149) q[1];
sx q[1];
rz(-0.74833732) q[1];
sx q[1];
rz(-1.8381717) q[1];
rz(-pi) q[2];
x q[2];
rz(2.73783) q[3];
sx q[3];
rz(-1.9039394) q[3];
sx q[3];
rz(-3.1333095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3489909) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(1.6646632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(1.1031716) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(-0.36956638) q[1];
sx q[1];
rz(-0.8820487) q[1];
sx q[1];
rz(-0.65912156) q[1];
rz(-1.5255899) q[2];
sx q[2];
rz(-0.87799413) q[2];
sx q[2];
rz(-2.9064657) q[2];
rz(2.7575708) q[3];
sx q[3];
rz(-2.1441318) q[3];
sx q[3];
rz(-1.287788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
