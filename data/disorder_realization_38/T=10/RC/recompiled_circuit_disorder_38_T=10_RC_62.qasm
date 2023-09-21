OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.881641) q[0];
sx q[0];
rz(-1.5437484) q[0];
sx q[0];
rz(-1.3511488) q[0];
rz(-pi) q[1];
rz(1.0371738) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(-1.2984315) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7632335) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(-0.50904973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.689765) q[3];
sx q[3];
rz(-0.61421466) q[3];
sx q[3];
rz(3.0005232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(-0.63981167) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(0.38133347) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249811) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(0.38831098) q[0];
x q[1];
rz(0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-2.6736205) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7035547) q[1];
sx q[1];
rz(-0.8982656) q[1];
sx q[1];
rz(2.9707675) q[1];
rz(-1.8566425) q[3];
sx q[3];
rz(-1.8010406) q[3];
sx q[3];
rz(0.86499351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(2.846068) q[0];
rz(0.23513901) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.495372) q[0];
sx q[0];
rz(-1.6618068) q[0];
sx q[0];
rz(-1.3489086) q[0];
rz(0.27772851) q[2];
sx q[2];
rz(-2.5979497) q[2];
sx q[2];
rz(3.1090528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5902139) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(-2.6022634) q[1];
rz(-pi) q[2];
rz(-1.3219222) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(-0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(0.27329683) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.8434974) q[0];
rz(-pi) q[1];
rz(0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(0.7960745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30444333) q[1];
sx q[1];
rz(-1.7177561) q[1];
sx q[1];
rz(1.3346346) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23394211) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(-3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894831) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(2.2228918) q[0];
rz(0.012533112) q[2];
sx q[2];
rz(-0.64385157) q[2];
sx q[2];
rz(-0.24521337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(2.4815464) q[1];
rz(-1.9636642) q[3];
sx q[3];
rz(-1.6097027) q[3];
sx q[3];
rz(1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50514454) q[0];
sx q[0];
rz(-0.2821869) q[0];
sx q[0];
rz(1.5178773) q[0];
rz(-pi) q[1];
rz(2.7369667) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(-0.65587703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2020859) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(0.034394666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2313314) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(-1.4753301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-0.0035704426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033747) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(0.72871448) q[0];
rz(-pi) q[1];
rz(1.1793169) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(-1.6187402) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4184985) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(2.198092) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3725029) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(-1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(-0.40813804) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5355276) q[0];
sx q[0];
rz(-1.3074271) q[0];
sx q[0];
rz(-2.8019964) q[0];
rz(2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(-1.0516143) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39290909) q[1];
sx q[1];
rz(-2.6137685) q[1];
sx q[1];
rz(-1.4455568) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7292761) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(0.52694595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066756847) q[0];
sx q[0];
rz(-1.2196676) q[0];
sx q[0];
rz(-1.8021276) q[0];
rz(-0.91295816) q[2];
sx q[2];
rz(-1.3557938) q[2];
sx q[2];
rz(-1.022033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.770551) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(0.71055926) q[1];
rz(-2.5832289) q[3];
sx q[3];
rz(-0.51279587) q[3];
sx q[3];
rz(-2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(0.3607761) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4496778) q[0];
sx q[0];
rz(-1.5458108) q[0];
sx q[0];
rz(-1.4033532) q[0];
x q[1];
rz(-1.785727) q[2];
sx q[2];
rz(-1.4533918) q[2];
sx q[2];
rz(1.4431151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68041486) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(-1.7635036) q[1];
x q[2];
rz(0.71513128) q[3];
sx q[3];
rz(-0.22278856) q[3];
sx q[3];
rz(0.95105329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5671134) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-2.0104682) q[2];
sx q[2];
rz(-2.9154073) q[2];
sx q[2];
rz(0.67868457) q[2];
rz(0.73137024) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
