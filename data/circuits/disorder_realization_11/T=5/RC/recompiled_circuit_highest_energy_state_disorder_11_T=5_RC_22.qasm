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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(-2.8859477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545171) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(0.70588995) q[0];
rz(-pi) q[1];
rz(-2.9886888) q[2];
sx q[2];
rz(-0.39275482) q[2];
sx q[2];
rz(-0.037234779) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5140443) q[1];
sx q[1];
rz(-1.6395634) q[1];
sx q[1];
rz(-1.3876983) q[1];
rz(-pi) q[2];
x q[2];
rz(2.643804) q[3];
sx q[3];
rz(-1.3827795) q[3];
sx q[3];
rz(1.1527485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6306182) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(-1.1154491) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(-2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9344591) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(-1.2472664) q[0];
rz(-2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(2.3242548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674448) q[0];
sx q[0];
rz(-0.56330452) q[0];
sx q[0];
rz(1.9356217) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1139718) q[2];
sx q[2];
rz(-2.6386542) q[2];
sx q[2];
rz(-2.3523112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25774792) q[1];
sx q[1];
rz(-0.93755975) q[1];
sx q[1];
rz(1.3650989) q[1];
rz(0.76089184) q[3];
sx q[3];
rz(-1.0659273) q[3];
sx q[3];
rz(0.99772108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3680129) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(-0.22029857) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050506266) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(-1.1385981) q[0];
rz(-2.7298722) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(-1.0134816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0419755) q[0];
sx q[0];
rz(-1.8238153) q[0];
sx q[0];
rz(1.5861041) q[0];
rz(3.0092952) q[2];
sx q[2];
rz(-2.0084642) q[2];
sx q[2];
rz(-2.3486111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0358082) q[1];
sx q[1];
rz(-2.2757451) q[1];
sx q[1];
rz(-2.1501599) q[1];
rz(-pi) q[2];
rz(-0.8441505) q[3];
sx q[3];
rz(-1.4733847) q[3];
sx q[3];
rz(-2.0908434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.960152) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(0.61035672) q[0];
rz(0.70862526) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(1.5707387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529396) q[0];
sx q[0];
rz(-1.4331237) q[0];
sx q[0];
rz(0.15451365) q[0];
rz(-pi) q[1];
rz(1.4518634) q[2];
sx q[2];
rz(-1.8557252) q[2];
sx q[2];
rz(0.97078427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.264894) q[1];
sx q[1];
rz(-1.8421409) q[1];
sx q[1];
rz(1.62074) q[1];
rz(-pi) q[2];
rz(0.39289423) q[3];
sx q[3];
rz(-2.3013982) q[3];
sx q[3];
rz(0.073176633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2107971) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50437462) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(0.24636191) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(-1.7074283) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.526699) q[0];
sx q[0];
rz(-1.5759567) q[0];
sx q[0];
rz(2.3847488) q[0];
x q[1];
rz(-0.44539659) q[2];
sx q[2];
rz(-1.7295618) q[2];
sx q[2];
rz(2.7743055) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38404462) q[1];
sx q[1];
rz(-1.9725807) q[1];
sx q[1];
rz(1.5662138) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1617817) q[3];
sx q[3];
rz(-0.24711619) q[3];
sx q[3];
rz(-0.97774001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.3110524) q[2];
sx q[2];
rz(-1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.865888) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-2.8979982) q[0];
rz(-0.98681915) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(-0.93596828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075386062) q[0];
sx q[0];
rz(-0.52156007) q[0];
sx q[0];
rz(-2.2564476) q[0];
rz(-pi) q[1];
rz(-0.34910874) q[2];
sx q[2];
rz(-0.85452467) q[2];
sx q[2];
rz(2.8385065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0430345) q[1];
sx q[1];
rz(-1.1698616) q[1];
sx q[1];
rz(-2.3844196) q[1];
x q[2];
rz(1.6020847) q[3];
sx q[3];
rz(-0.60006053) q[3];
sx q[3];
rz(-0.72573435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(-2.7916419) q[2];
rz(2.5838666) q[3];
sx q[3];
rz(-1.9316659) q[3];
sx q[3];
rz(-0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418942) q[0];
sx q[0];
rz(-0.63755578) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(2.8217577) q[1];
sx q[1];
rz(-1.6022857) q[1];
sx q[1];
rz(1.1357657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1162735) q[0];
sx q[0];
rz(-0.92705446) q[0];
sx q[0];
rz(-1.5035065) q[0];
rz(-2.6276845) q[2];
sx q[2];
rz(-2.137988) q[2];
sx q[2];
rz(2.6197768) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3150683) q[1];
sx q[1];
rz(-0.60501912) q[1];
sx q[1];
rz(-1.7111163) q[1];
rz(-pi) q[2];
rz(-2.3507422) q[3];
sx q[3];
rz(-1.6501657) q[3];
sx q[3];
rz(0.96109238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4096421) q[2];
sx q[2];
rz(-2.4739517) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(0.55602443) q[0];
rz(0.22932886) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(-1.75846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616104) q[0];
sx q[0];
rz(-1.605827) q[0];
sx q[0];
rz(2.4039414) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5117308) q[2];
sx q[2];
rz(-1.1129654) q[2];
sx q[2];
rz(-0.77564592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(-0.31633693) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9541627) q[3];
sx q[3];
rz(-1.1463506) q[3];
sx q[3];
rz(2.2147873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(-1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.7187985) q[3];
sx q[3];
rz(-0.82714287) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(0.14778368) q[0];
rz(2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145549) q[0];
sx q[0];
rz(-2.6588024) q[0];
sx q[0];
rz(-2.618769) q[0];
rz(-2.0835593) q[2];
sx q[2];
rz(-2.4840925) q[2];
sx q[2];
rz(3.0160144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3264824) q[1];
sx q[1];
rz(-1.8811783) q[1];
sx q[1];
rz(1.2182477) q[1];
rz(-pi) q[2];
rz(1.4465973) q[3];
sx q[3];
rz(-2.327965) q[3];
sx q[3];
rz(-0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96616894) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(-1.8625205) q[2];
rz(-1.4281979) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3928423) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(-3.0774935) q[0];
rz(-0.52866689) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(0.97506964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25697432) q[0];
sx q[0];
rz(-1.3064837) q[0];
sx q[0];
rz(-0.5075016) q[0];
rz(-pi) q[1];
rz(2.1553173) q[2];
sx q[2];
rz(-2.0615675) q[2];
sx q[2];
rz(-1.7391537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9916315) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(-0.44709713) q[1];
rz(-pi) q[2];
rz(1.7998135) q[3];
sx q[3];
rz(-1.8088565) q[3];
sx q[3];
rz(2.2464667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14572445) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(0.60689849) q[2];
rz(2.8912344) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(0.50601602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470227) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(-0.49867123) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-2.2398938) q[2];
sx q[2];
rz(-0.82868262) q[2];
sx q[2];
rz(2.5862624) q[2];
rz(-1.7767033) q[3];
sx q[3];
rz(-1.4510703) q[3];
sx q[3];
rz(2.5616796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
