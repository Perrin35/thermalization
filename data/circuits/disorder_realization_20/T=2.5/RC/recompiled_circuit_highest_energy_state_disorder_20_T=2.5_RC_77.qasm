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
rz(-0.72467726) q[0];
sx q[0];
rz(4.4530498) q[0];
sx q[0];
rz(11.761576) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(6.4810452) q[1];
sx q[1];
rz(5.0623464) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22053424) q[0];
sx q[0];
rz(-3.0166114) q[0];
sx q[0];
rz(0.28470953) q[0];
rz(-pi) q[1];
rz(2.8684565) q[2];
sx q[2];
rz(-0.47800666) q[2];
sx q[2];
rz(-0.82551685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3343494) q[1];
sx q[1];
rz(-0.74485272) q[1];
sx q[1];
rz(-3.0762385) q[1];
x q[2];
rz(-1.4335459) q[3];
sx q[3];
rz(-0.94232925) q[3];
sx q[3];
rz(-0.31627218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7120984) q[2];
sx q[2];
rz(-2.6338989) q[2];
sx q[2];
rz(1.3362159) q[2];
rz(1.7441689) q[3];
sx q[3];
rz(-1.6349399) q[3];
sx q[3];
rz(1.5857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812688) q[0];
sx q[0];
rz(-1.5070494) q[0];
sx q[0];
rz(-0.67833483) q[0];
rz(2.5256269) q[1];
sx q[1];
rz(-1.0266285) q[1];
sx q[1];
rz(-0.14332992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5472488) q[0];
sx q[0];
rz(-1.4137725) q[0];
sx q[0];
rz(-1.2444522) q[0];
rz(-pi) q[1];
rz(0.8524266) q[2];
sx q[2];
rz(-1.4839236) q[2];
sx q[2];
rz(-1.014682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.424984) q[1];
sx q[1];
rz(-2.0758481) q[1];
sx q[1];
rz(-1.6249955) q[1];
rz(-pi) q[2];
rz(-1.2202699) q[3];
sx q[3];
rz(-1.1683147) q[3];
sx q[3];
rz(0.19798812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6501288) q[2];
sx q[2];
rz(-2.3626707) q[2];
sx q[2];
rz(1.2196563) q[2];
rz(2.8236112) q[3];
sx q[3];
rz(-2.3173083) q[3];
sx q[3];
rz(0.73941755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3984482) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(2.5019124) q[0];
rz(-1.3321053) q[1];
sx q[1];
rz(-0.94146856) q[1];
sx q[1];
rz(2.926362) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71986854) q[0];
sx q[0];
rz(-2.6009237) q[0];
sx q[0];
rz(0.42406492) q[0];
x q[1];
rz(1.5528312) q[2];
sx q[2];
rz(-2.8572272) q[2];
sx q[2];
rz(0.21321061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58956026) q[1];
sx q[1];
rz(-1.5812687) q[1];
sx q[1];
rz(-1.1754491) q[1];
x q[2];
rz(-0.12037511) q[3];
sx q[3];
rz(-1.1167428) q[3];
sx q[3];
rz(-1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3053863) q[2];
sx q[2];
rz(-0.47577327) q[2];
sx q[2];
rz(-2.9212908) q[2];
rz(0.78271714) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(-0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7102605) q[0];
sx q[0];
rz(-0.98589698) q[0];
sx q[0];
rz(-1.8248935) q[0];
rz(-2.6915164) q[1];
sx q[1];
rz(-1.7066259) q[1];
sx q[1];
rz(-3.0302474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8985626) q[0];
sx q[0];
rz(-0.70625171) q[0];
sx q[0];
rz(-0.39470048) q[0];
x q[1];
rz(-1.8931023) q[2];
sx q[2];
rz(-1.4864902) q[2];
sx q[2];
rz(1.5706667) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.131241) q[1];
sx q[1];
rz(-1.2919754) q[1];
sx q[1];
rz(-2.1479021) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5841949) q[3];
sx q[3];
rz(-1.9950997) q[3];
sx q[3];
rz(-2.4941924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4075809) q[2];
sx q[2];
rz(-1.13052) q[2];
sx q[2];
rz(-3.0111175) q[2];
rz(2.981251) q[3];
sx q[3];
rz(-2.4800143) q[3];
sx q[3];
rz(-2.9714238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570356) q[0];
sx q[0];
rz(-2.9354876) q[0];
sx q[0];
rz(3.0233132) q[0];
rz(-2.2453399) q[1];
sx q[1];
rz(-1.6629442) q[1];
sx q[1];
rz(0.55312696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.129707) q[0];
sx q[0];
rz(-1.0965523) q[0];
sx q[0];
rz(1.9408731) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4055355) q[2];
sx q[2];
rz(-1.2345353) q[2];
sx q[2];
rz(-2.2081809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0748079) q[1];
sx q[1];
rz(-0.60943595) q[1];
sx q[1];
rz(1.8685196) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0703662) q[3];
sx q[3];
rz(-2.4132846) q[3];
sx q[3];
rz(-2.0508931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3502189) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(-2.1283894) q[2];
rz(-0.31747097) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1182227) q[0];
sx q[0];
rz(-2.2582) q[0];
sx q[0];
rz(0.50514847) q[0];
rz(-0.39816868) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(-1.8123951) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4666207) q[0];
sx q[0];
rz(-1.734991) q[0];
sx q[0];
rz(1.0999099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6077725) q[2];
sx q[2];
rz(-1.2580192) q[2];
sx q[2];
rz(0.61543265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8230629) q[1];
sx q[1];
rz(-1.6803871) q[1];
sx q[1];
rz(2.1353882) q[1];
rz(2.8839467) q[3];
sx q[3];
rz(-1.8339673) q[3];
sx q[3];
rz(-1.3130696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.087611) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(0.83016738) q[2];
rz(-1.6603445) q[3];
sx q[3];
rz(-2.0989428) q[3];
sx q[3];
rz(-1.164485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2723349) q[0];
sx q[0];
rz(-0.91599688) q[0];
sx q[0];
rz(1.2087615) q[0];
rz(2.8973978) q[1];
sx q[1];
rz(-1.1714275) q[1];
sx q[1];
rz(-0.29921439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41847412) q[0];
sx q[0];
rz(-1.5970073) q[0];
sx q[0];
rz(-1.3546726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1019548) q[2];
sx q[2];
rz(-2.1130145) q[2];
sx q[2];
rz(-1.7249267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4296885) q[1];
sx q[1];
rz(-1.619408) q[1];
sx q[1];
rz(0.91829388) q[1];
x q[2];
rz(-2.4387909) q[3];
sx q[3];
rz(-1.2497823) q[3];
sx q[3];
rz(1.7197844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19305912) q[2];
sx q[2];
rz(-1.4024573) q[2];
sx q[2];
rz(-2.194727) q[2];
rz(1.7638505) q[3];
sx q[3];
rz(-2.6008714) q[3];
sx q[3];
rz(-0.64259678) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4544025) q[0];
sx q[0];
rz(-0.43455046) q[0];
sx q[0];
rz(0.53892556) q[0];
rz(-2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(2.7573746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095188524) q[0];
sx q[0];
rz(-0.47391072) q[0];
sx q[0];
rz(-2.7907811) q[0];
rz(-pi) q[1];
rz(2.8649533) q[2];
sx q[2];
rz(-2.5211272) q[2];
sx q[2];
rz(2.7850604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3620262) q[1];
sx q[1];
rz(-3.0770731) q[1];
sx q[1];
rz(1.9883335) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0168276) q[3];
sx q[3];
rz(-1.5938543) q[3];
sx q[3];
rz(-0.35820828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82181278) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(-2.5452851) q[2];
rz(2.1042306) q[3];
sx q[3];
rz(-2.3682902) q[3];
sx q[3];
rz(1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7647112) q[0];
sx q[0];
rz(-2.3887971) q[0];
sx q[0];
rz(-2.9994614) q[0];
rz(1.1171974) q[1];
sx q[1];
rz(-1.486472) q[1];
sx q[1];
rz(1.2275009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1466897) q[0];
sx q[0];
rz(-2.0462543) q[0];
sx q[0];
rz(-1.2593996) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1255163) q[2];
sx q[2];
rz(-1.4332891) q[2];
sx q[2];
rz(-1.2924271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0490659) q[1];
sx q[1];
rz(-2.5124308) q[1];
sx q[1];
rz(1.4607304) q[1];
x q[2];
rz(2.3120425) q[3];
sx q[3];
rz(-1.2365474) q[3];
sx q[3];
rz(0.087141872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.050798945) q[2];
sx q[2];
rz(-0.24238452) q[2];
sx q[2];
rz(2.3432815) q[2];
rz(1.5618886) q[3];
sx q[3];
rz(-1.7589933) q[3];
sx q[3];
rz(-0.17670512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.056219) q[0];
sx q[0];
rz(-2.4985785) q[0];
sx q[0];
rz(0.0052848919) q[0];
rz(3.058744) q[1];
sx q[1];
rz(-2.0382035) q[1];
sx q[1];
rz(-2.8740035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2076715) q[0];
sx q[0];
rz(-0.90822847) q[0];
sx q[0];
rz(-1.1712058) q[0];
rz(1.1830953) q[2];
sx q[2];
rz(-1.8418717) q[2];
sx q[2];
rz(2.3623206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5983683) q[1];
sx q[1];
rz(-2.7562047) q[1];
sx q[1];
rz(2.1699127) q[1];
rz(-pi) q[2];
rz(2.0844719) q[3];
sx q[3];
rz(-1.1506211) q[3];
sx q[3];
rz(0.021180245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.573632) q[2];
sx q[2];
rz(-0.48144123) q[2];
sx q[2];
rz(1.2130223) q[2];
rz(-1.2447478) q[3];
sx q[3];
rz(-0.48143482) q[3];
sx q[3];
rz(1.9511694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5382814) q[0];
sx q[0];
rz(-0.9700226) q[0];
sx q[0];
rz(0.097401311) q[0];
rz(0.010460336) q[1];
sx q[1];
rz(-2.2757826) q[1];
sx q[1];
rz(-0.59715685) q[1];
rz(2.6059628) q[2];
sx q[2];
rz(-1.6471072) q[2];
sx q[2];
rz(-2.9181388) q[2];
rz(0.98164557) q[3];
sx q[3];
rz(-1.7170148) q[3];
sx q[3];
rz(-0.96819504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
