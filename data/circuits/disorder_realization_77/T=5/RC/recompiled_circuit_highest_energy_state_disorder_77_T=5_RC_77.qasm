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
rz(2.0692985) q[0];
sx q[0];
rz(-2.9579853) q[0];
sx q[0];
rz(-2.973383) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(1.1163196) q[1];
sx q[1];
rz(7.4012227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7752286) q[0];
sx q[0];
rz(-0.34465955) q[0];
sx q[0];
rz(2.8135236) q[0];
rz(-pi) q[1];
rz(1.9498212) q[2];
sx q[2];
rz(-2.2281335) q[2];
sx q[2];
rz(2.71703) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.019844481) q[1];
sx q[1];
rz(-2.3308874) q[1];
sx q[1];
rz(1.5365063) q[1];
x q[2];
rz(1.8421451) q[3];
sx q[3];
rz(-1.4381471) q[3];
sx q[3];
rz(2.1675183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85718021) q[2];
sx q[2];
rz(-0.99049157) q[2];
sx q[2];
rz(-1.8641137) q[2];
rz(0.60753456) q[3];
sx q[3];
rz(-1.7888864) q[3];
sx q[3];
rz(-1.7551306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4472537) q[0];
sx q[0];
rz(-0.67830938) q[0];
sx q[0];
rz(0.69764486) q[0];
rz(2.8088226) q[1];
sx q[1];
rz(-2.0336626) q[1];
sx q[1];
rz(0.3322126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8150114) q[0];
sx q[0];
rz(-1.7012811) q[0];
sx q[0];
rz(-1.6925808) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23770419) q[2];
sx q[2];
rz(-2.1344886) q[2];
sx q[2];
rz(1.7089882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.028113508) q[1];
sx q[1];
rz(-1.3225127) q[1];
sx q[1];
rz(-2.0015654) q[1];
rz(-1.148613) q[3];
sx q[3];
rz(-2.2530972) q[3];
sx q[3];
rz(-2.9158031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.701) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(2.8458703) q[2];
rz(-0.010585636) q[3];
sx q[3];
rz(-0.9897832) q[3];
sx q[3];
rz(-2.4081965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6538786) q[0];
sx q[0];
rz(-2.6831477) q[0];
sx q[0];
rz(2.606875) q[0];
rz(2.0405105) q[1];
sx q[1];
rz(-1.705876) q[1];
sx q[1];
rz(-2.2284177) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951626) q[0];
sx q[0];
rz(-2.5341153) q[0];
sx q[0];
rz(0.98275252) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1712672) q[2];
sx q[2];
rz(-2.5472884) q[2];
sx q[2];
rz(-0.28501302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5194381) q[1];
sx q[1];
rz(-1.039912) q[1];
sx q[1];
rz(0.77150788) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2997268) q[3];
sx q[3];
rz(-1.6252691) q[3];
sx q[3];
rz(-1.5845875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59213006) q[2];
sx q[2];
rz(-0.89207804) q[2];
sx q[2];
rz(2.7590052) q[2];
rz(-1.654918) q[3];
sx q[3];
rz(-0.29522172) q[3];
sx q[3];
rz(0.91992831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.53052467) q[0];
sx q[0];
rz(-1.4570197) q[0];
sx q[0];
rz(-1.0462421) q[0];
rz(2.5776082) q[1];
sx q[1];
rz(-2.1395186) q[1];
sx q[1];
rz(-1.8246338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3377258) q[0];
sx q[0];
rz(-0.90529862) q[0];
sx q[0];
rz(-2.4212877) q[0];
x q[1];
rz(-0.91805075) q[2];
sx q[2];
rz(-1.3155283) q[2];
sx q[2];
rz(1.4087311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7822588) q[1];
sx q[1];
rz(-0.31681199) q[1];
sx q[1];
rz(-3.003939) q[1];
x q[2];
rz(2.7183771) q[3];
sx q[3];
rz(-0.82941705) q[3];
sx q[3];
rz(-2.2256782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4468533) q[2];
sx q[2];
rz(-1.4444192) q[2];
sx q[2];
rz(-0.99855885) q[2];
rz(2.5616732) q[3];
sx q[3];
rz(-1.6603575) q[3];
sx q[3];
rz(-1.8370321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21438433) q[0];
sx q[0];
rz(-1.2283607) q[0];
sx q[0];
rz(2.4217915) q[0];
rz(0.7431227) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(-0.057083759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.538495) q[0];
sx q[0];
rz(-0.39347005) q[0];
sx q[0];
rz(-0.35719494) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.916831) q[2];
sx q[2];
rz(-2.2979696) q[2];
sx q[2];
rz(0.42842254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47810208) q[1];
sx q[1];
rz(-1.2896207) q[1];
sx q[1];
rz(-1.4040787) q[1];
rz(-pi) q[2];
rz(-2.1674183) q[3];
sx q[3];
rz(-1.2032185) q[3];
sx q[3];
rz(-1.1413464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4382978) q[2];
sx q[2];
rz(-1.9652099) q[2];
sx q[2];
rz(-1.2698184) q[2];
rz(-2.5979009) q[3];
sx q[3];
rz(-0.68259493) q[3];
sx q[3];
rz(1.773275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.066034) q[0];
sx q[0];
rz(-2.7664001) q[0];
sx q[0];
rz(0.22848465) q[0];
rz(-0.98555073) q[1];
sx q[1];
rz(-1.5551714) q[1];
sx q[1];
rz(-1.9353345) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0096480308) q[0];
sx q[0];
rz(-1.4033269) q[0];
sx q[0];
rz(0.33958774) q[0];
rz(-2.9183675) q[2];
sx q[2];
rz(-1.8656261) q[2];
sx q[2];
rz(-2.1708174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6724193) q[1];
sx q[1];
rz(-2.1757728) q[1];
sx q[1];
rz(0.14766001) q[1];
rz(-pi) q[2];
rz(-2.966622) q[3];
sx q[3];
rz(-0.55040437) q[3];
sx q[3];
rz(-3.1249078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1365405) q[2];
sx q[2];
rz(-1.8661934) q[2];
sx q[2];
rz(1.8955815) q[2];
rz(0.15240845) q[3];
sx q[3];
rz(-2.992575) q[3];
sx q[3];
rz(-0.77308956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3054672) q[0];
sx q[0];
rz(-1.9274599) q[0];
sx q[0];
rz(0.25492302) q[0];
rz(2.1615324) q[1];
sx q[1];
rz(-1.7995116) q[1];
sx q[1];
rz(0.23475383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8844719) q[0];
sx q[0];
rz(-2.5354374) q[0];
sx q[0];
rz(-1.1324833) q[0];
x q[1];
rz(1.8053986) q[2];
sx q[2];
rz(-1.4035181) q[2];
sx q[2];
rz(-1.7931615) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3666881) q[1];
sx q[1];
rz(-1.8340381) q[1];
sx q[1];
rz(0.4270435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8942306) q[3];
sx q[3];
rz(-1.4787592) q[3];
sx q[3];
rz(2.2313554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0945101) q[2];
sx q[2];
rz(-0.42028704) q[2];
sx q[2];
rz(-1.9455998) q[2];
rz(0.99810537) q[3];
sx q[3];
rz(-1.2752897) q[3];
sx q[3];
rz(-0.78817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.3274479) q[0];
sx q[0];
rz(-0.87004167) q[0];
sx q[0];
rz(-2.5092464) q[0];
rz(-1.5898204) q[1];
sx q[1];
rz(-1.7951782) q[1];
sx q[1];
rz(-1.8428892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017865861) q[0];
sx q[0];
rz(-1.3204823) q[0];
sx q[0];
rz(1.1160442) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66565158) q[2];
sx q[2];
rz(-1.5857664) q[2];
sx q[2];
rz(2.7240208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3159193) q[1];
sx q[1];
rz(-0.50858595) q[1];
sx q[1];
rz(-3.0561165) q[1];
rz(-0.18142314) q[3];
sx q[3];
rz(-2.9876935) q[3];
sx q[3];
rz(-1.5937663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.27114) q[2];
sx q[2];
rz(-2.2282034) q[2];
sx q[2];
rz(-2.3823628) q[2];
rz(0.40902725) q[3];
sx q[3];
rz(-1.7274011) q[3];
sx q[3];
rz(-0.11212382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8522393) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(2.7759743) q[0];
rz(-0.002937142) q[1];
sx q[1];
rz(-1.5298839) q[1];
sx q[1];
rz(-2.073435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.070052) q[0];
sx q[0];
rz(-2.5134633) q[0];
sx q[0];
rz(-0.48372901) q[0];
rz(0.67068213) q[2];
sx q[2];
rz(-2.1895231) q[2];
sx q[2];
rz(0.77341341) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0823487) q[1];
sx q[1];
rz(-1.6089393) q[1];
sx q[1];
rz(1.1831985) q[1];
x q[2];
rz(1.5434331) q[3];
sx q[3];
rz(-1.7723736) q[3];
sx q[3];
rz(1.8044173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3785765) q[2];
sx q[2];
rz(-0.91332674) q[2];
sx q[2];
rz(2.1545289) q[2];
rz(2.3297564) q[3];
sx q[3];
rz(-0.87849179) q[3];
sx q[3];
rz(-1.2442376) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064747485) q[0];
sx q[0];
rz(-0.61926121) q[0];
sx q[0];
rz(1.6774119) q[0];
rz(-2.1773416) q[1];
sx q[1];
rz(-1.3136656) q[1];
sx q[1];
rz(-1.7589232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9486757) q[0];
sx q[0];
rz(-2.3423704) q[0];
sx q[0];
rz(-2.1069631) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8455454) q[2];
sx q[2];
rz(-1.2128069) q[2];
sx q[2];
rz(-1.6894947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5556293) q[1];
sx q[1];
rz(-2.4307287) q[1];
sx q[1];
rz(2.2522624) q[1];
rz(-pi) q[2];
rz(-0.85282214) q[3];
sx q[3];
rz(-0.75197619) q[3];
sx q[3];
rz(-0.2096006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1796403) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(-1.4947653) q[2];
rz(1.1163813) q[3];
sx q[3];
rz(-0.79018441) q[3];
sx q[3];
rz(-2.5377929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8981001) q[0];
sx q[0];
rz(-1.6742764) q[0];
sx q[0];
rz(1.8754638) q[0];
rz(2.948214) q[1];
sx q[1];
rz(-1.0693751) q[1];
sx q[1];
rz(1.3042915) q[1];
rz(-1.8776352) q[2];
sx q[2];
rz(-0.86021341) q[2];
sx q[2];
rz(-2.025069) q[2];
rz(-2.4520649) q[3];
sx q[3];
rz(-1.22682) q[3];
sx q[3];
rz(0.20709909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
