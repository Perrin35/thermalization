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
rz(-0.25804582) q[0];
rz(-0.99335042) q[1];
sx q[1];
rz(-1.8408096) q[1];
sx q[1];
rz(2.8859477) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50798049) q[0];
sx q[0];
rz(-0.88417142) q[0];
sx q[0];
rz(-0.79721862) q[0];
rz(-pi) q[1];
rz(-1.5077816) q[2];
sx q[2];
rz(-1.1828701) q[2];
sx q[2];
rz(0.20252075) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5140443) q[1];
sx q[1];
rz(-1.5020292) q[1];
sx q[1];
rz(-1.7538944) q[1];
rz(1.7840476) q[3];
sx q[3];
rz(-1.0825601) q[3];
sx q[3];
rz(-0.31682107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6306182) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(2.0261436) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(2.7052963) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2071335) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(-2.3242548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674448) q[0];
sx q[0];
rz(-0.56330452) q[0];
sx q[0];
rz(1.9356217) q[0];
rz(-pi) q[1];
rz(1.1306612) q[2];
sx q[2];
rz(-1.8225742) q[2];
sx q[2];
rz(1.8735069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8838447) q[1];
sx q[1];
rz(-0.93755975) q[1];
sx q[1];
rz(-1.7764938) q[1];
rz(-pi) q[2];
rz(2.3807008) q[3];
sx q[3];
rz(-2.0756654) q[3];
sx q[3];
rz(0.99772108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77357972) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(-2.9212941) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(2.7298722) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(1.0134816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47501144) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(0.25304742) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2960245) q[2];
sx q[2];
rz(-0.45598511) q[2];
sx q[2];
rz(-0.48874654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0358082) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(-2.1501599) q[1];
rz(-pi) q[2];
rz(-3.0115836) q[3];
sx q[3];
rz(-2.293236) q[3];
sx q[3];
rz(-0.60628451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98075214) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(-2.8694425) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(-0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(-2.5312359) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(1.5707387) q[1];
rz(-pi) q[2];
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
rz(-1.8557252) q[2];
sx q[2];
rz(0.97078427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8766986) q[1];
sx q[1];
rz(-1.8421409) q[1];
sx q[1];
rz(1.5208526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39289423) q[3];
sx q[3];
rz(-0.84019444) q[3];
sx q[3];
rz(-0.073176633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2107971) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(0.24636191) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(-1.7074283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1808162) q[0];
sx q[0];
rz(-0.81396507) q[0];
sx q[0];
rz(-1.5636982) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7463914) q[2];
sx q[2];
rz(-1.1313952) q[2];
sx q[2];
rz(-2.0134157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38404462) q[1];
sx q[1];
rz(-1.169012) q[1];
sx q[1];
rz(1.5753788) q[1];
rz(-pi) q[2];
rz(-1.7982539) q[3];
sx q[3];
rz(-1.47336) q[3];
sx q[3];
rz(0.9909329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65115923) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(-1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(-1.2580416) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27570462) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(0.24359447) q[0];
rz(-0.98681915) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(2.2056244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87847951) q[0];
sx q[0];
rz(-1.2498444) q[0];
sx q[0];
rz(1.1522989) q[0];
x q[1];
rz(-2.317993) q[2];
sx q[2];
rz(-1.3098426) q[2];
sx q[2];
rz(-1.5023155) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0985581) q[1];
sx q[1];
rz(-1.1698616) q[1];
sx q[1];
rz(-2.3844196) q[1];
rz(3.1201911) q[3];
sx q[3];
rz(-2.170522) q[3];
sx q[3];
rz(0.76364005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.662107) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(2.7916419) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.9316659) q[3];
sx q[3];
rz(0.85132712) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418942) q[0];
sx q[0];
rz(-0.63755578) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(-2.8217577) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(1.1357657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58590305) q[0];
sx q[0];
rz(-1.6246038) q[0];
sx q[0];
rz(-2.4967628) q[0];
rz(-pi) q[1];
rz(0.93923969) q[2];
sx q[2];
rz(-1.9983091) q[2];
sx q[2];
rz(-1.7981426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7702076) q[1];
sx q[1];
rz(-1.6504297) q[1];
sx q[1];
rz(-0.97038986) q[1];
rz(-pi) q[2];
rz(0.11140996) q[3];
sx q[3];
rz(-2.3476331) q[3];
sx q[3];
rz(-0.68796989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(-2.9988334) q[2];
rz(1.3879294) q[3];
sx q[3];
rz(-1.1889435) q[3];
sx q[3];
rz(2.4587542) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1801572) q[0];
sx q[0];
rz(-3.1249983) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(1.75846) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55229462) q[0];
sx q[0];
rz(-2.4032666) q[0];
sx q[0];
rz(3.0895322) q[0];
rz(-3.0223614) q[2];
sx q[2];
rz(-0.46135739) q[2];
sx q[2];
rz(-2.2329494) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8361159) q[1];
sx q[1];
rz(-2.6295942) q[1];
sx q[1];
rz(-2.192537) q[1];
x q[2];
rz(-1.1874299) q[3];
sx q[3];
rz(-1.995242) q[3];
sx q[3];
rz(2.2147873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(-1.2410835) q[2];
rz(-3.0789913) q[3];
sx q[3];
rz(-1.7187985) q[3];
sx q[3];
rz(2.3144498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(-0.14778368) q[0];
rz(-0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157279) q[0];
sx q[0];
rz(-1.3368538) q[0];
sx q[0];
rz(0.42629231) q[0];
rz(0.36208533) q[2];
sx q[2];
rz(-2.1323983) q[2];
sx q[2];
rz(-0.74383273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7853773) q[1];
sx q[1];
rz(-1.9058203) q[1];
sx q[1];
rz(0.32932333) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13033615) q[3];
sx q[3];
rz(-2.3763083) q[3];
sx q[3];
rz(1.9627067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.2790722) q[2];
rz(1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(3.0774935) q[0];
rz(0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(-2.166523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720358) q[0];
sx q[0];
rz(-2.0590879) q[0];
sx q[0];
rz(-1.2704865) q[0];
rz(-pi) q[1];
rz(0.80143546) q[2];
sx q[2];
rz(-2.3972627) q[2];
sx q[2];
rz(-2.3542885) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31723695) q[1];
sx q[1];
rz(-2.0065161) q[1];
sx q[1];
rz(-1.3292666) q[1];
rz(-2.3894839) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(-0.11550918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9958682) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(-0.25035826) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(-2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.09457) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-0.90169883) q[2];
sx q[2];
rz(-2.31291) q[2];
sx q[2];
rz(-0.55533021) q[2];
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
