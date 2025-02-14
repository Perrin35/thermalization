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
rz(2.430001) q[0];
sx q[0];
rz(-1.1049668) q[0];
sx q[0];
rz(2.0018863) q[0];
rz(-0.79144129) q[1];
sx q[1];
rz(-1.0410407) q[1];
sx q[1];
rz(2.0120373) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.043046) q[0];
sx q[0];
rz(-2.401863) q[0];
sx q[0];
rz(-1.84124) q[0];
x q[1];
rz(3.1055345) q[2];
sx q[2];
rz(-2.4781422) q[2];
sx q[2];
rz(-0.61095881) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1892231) q[1];
sx q[1];
rz(-0.86878759) q[1];
sx q[1];
rz(2.7527806) q[1];
rz(3.0409052) q[3];
sx q[3];
rz(-0.33439454) q[3];
sx q[3];
rz(-1.3833801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1019885) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(0.25092009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(3.0372341) q[0];
rz(1.9508427) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(-2.6844535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16743827) q[0];
sx q[0];
rz(-1.3773241) q[0];
sx q[0];
rz(-3.1102212) q[0];
rz(-pi) q[1];
rz(-2.1188583) q[2];
sx q[2];
rz(-0.75621683) q[2];
sx q[2];
rz(-0.05847419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7105661) q[1];
sx q[1];
rz(-1.2557286) q[1];
sx q[1];
rz(-3.0418212) q[1];
rz(-0.13808226) q[3];
sx q[3];
rz(-2.2347817) q[3];
sx q[3];
rz(0.21747227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55598688) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(0.15698329) q[2];
rz(0.5213151) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(-3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3954725) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-2.795862) q[0];
rz(-1.455447) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.5706496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5356876) q[0];
sx q[0];
rz(-0.55258026) q[0];
sx q[0];
rz(-1.7406169) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5657896) q[2];
sx q[2];
rz(-2.5189159) q[2];
sx q[2];
rz(-1.1345991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6863043) q[1];
sx q[1];
rz(-1.6197816) q[1];
sx q[1];
rz(2.9279885) q[1];
rz(2.219958) q[3];
sx q[3];
rz(-1.7329669) q[3];
sx q[3];
rz(2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7741144) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(2.1841614) q[2];
rz(-2.0188324) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(-1.3721589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51266176) q[0];
sx q[0];
rz(-1.7191732) q[0];
sx q[0];
rz(1.5976394) q[0];
rz(1.3900025) q[1];
sx q[1];
rz(-1.8252239) q[1];
sx q[1];
rz(-1.66473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4721253) q[0];
sx q[0];
rz(-0.60807121) q[0];
sx q[0];
rz(-2.6345992) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0576287) q[2];
sx q[2];
rz(-1.3332841) q[2];
sx q[2];
rz(0.31667865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37472607) q[1];
sx q[1];
rz(-2.0927168) q[1];
sx q[1];
rz(2.7609227) q[1];
rz(-2.8119726) q[3];
sx q[3];
rz(-1.8330036) q[3];
sx q[3];
rz(0.58673687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59820286) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(0.99679917) q[2];
rz(1.8761084) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0604414) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(-1.0720217) q[0];
rz(2.187166) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(-1.4929474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38781751) q[0];
sx q[0];
rz(-0.79083453) q[0];
sx q[0];
rz(-2.8641939) q[0];
rz(-pi) q[1];
rz(2.4193108) q[2];
sx q[2];
rz(-1.9123291) q[2];
sx q[2];
rz(2.3010437) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6750329) q[1];
sx q[1];
rz(-1.7132812) q[1];
sx q[1];
rz(-2.811048) q[1];
x q[2];
rz(0.32060726) q[3];
sx q[3];
rz(-1.6658452) q[3];
sx q[3];
rz(-2.7966201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2299049) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(-1.9913199) q[2];
rz(3.1005499) q[3];
sx q[3];
rz(-1.5464455) q[3];
sx q[3];
rz(-0.40157792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3732442) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.1596229) q[0];
rz(1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.9083091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4077743) q[0];
sx q[0];
rz(-2.3253157) q[0];
sx q[0];
rz(0.4305779) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6242113) q[2];
sx q[2];
rz(-1.6279977) q[2];
sx q[2];
rz(-0.085669667) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6813876) q[1];
sx q[1];
rz(-2.3355995) q[1];
sx q[1];
rz(0.56350033) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43253492) q[3];
sx q[3];
rz(-2.5721241) q[3];
sx q[3];
rz(-1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50501862) q[2];
sx q[2];
rz(-2.7754112) q[2];
sx q[2];
rz(-1.150307) q[2];
rz(-2.4023174) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4742541) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(2.6928103) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-2.0346784) q[1];
sx q[1];
rz(-1.1610228) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5509508) q[0];
sx q[0];
rz(-1.30473) q[0];
sx q[0];
rz(-2.7177285) q[0];
rz(-pi) q[1];
rz(2.3660701) q[2];
sx q[2];
rz(-2.4897414) q[2];
sx q[2];
rz(-1.8076123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22299448) q[1];
sx q[1];
rz(-1.9004993) q[1];
sx q[1];
rz(-1.2597643) q[1];
rz(-pi) q[2];
rz(1.0169898) q[3];
sx q[3];
rz(-1.7968868) q[3];
sx q[3];
rz(0.46976058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9728969) q[2];
sx q[2];
rz(-1.8099512) q[2];
sx q[2];
rz(-1.1807582) q[2];
rz(2.0598038) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(2.7150174) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54843724) q[0];
sx q[0];
rz(-1.3189545) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(-1.1697191) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(0.39264548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8317345) q[0];
sx q[0];
rz(-1.4295973) q[0];
sx q[0];
rz(-1.0910499) q[0];
x q[1];
rz(-2.8296109) q[2];
sx q[2];
rz(-1.4380961) q[2];
sx q[2];
rz(-0.080094425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.709014) q[1];
sx q[1];
rz(-0.48332941) q[1];
sx q[1];
rz(2.9523115) q[1];
rz(1.6527376) q[3];
sx q[3];
rz(-1.7590766) q[3];
sx q[3];
rz(0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9795867) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-0.77097547) q[3];
sx q[3];
rz(0.60105598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9252121) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(1.3339169) q[0];
rz(1.1434309) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(-0.7935895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5549907) q[0];
sx q[0];
rz(-1.7287917) q[0];
sx q[0];
rz(1.642307) q[0];
rz(-pi) q[1];
rz(-0.066448224) q[2];
sx q[2];
rz(-2.1796103) q[2];
sx q[2];
rz(-2.7087536) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55793437) q[1];
sx q[1];
rz(-1.0685295) q[1];
sx q[1];
rz(1.7947547) q[1];
rz(-1.3511485) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.8203514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44635832) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(1.0515155) q[2];
rz(0.71600437) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-2.9606384) q[0];
rz(1.1894233) q[1];
sx q[1];
rz(-1.1782497) q[1];
sx q[1];
rz(-2.1612371) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086209379) q[0];
sx q[0];
rz(-2.6329941) q[0];
sx q[0];
rz(2.9706756) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3848214) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(2.7293918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7804523) q[1];
sx q[1];
rz(-1.3899511) q[1];
sx q[1];
rz(2.3922582) q[1];
rz(1.7813563) q[3];
sx q[3];
rz(-1.9580152) q[3];
sx q[3];
rz(0.59084597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.02562) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.8170961) q[2];
rz(-2.4888511) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(0.038221922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41435913) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(-0.7242135) q[1];
sx q[1];
rz(-2.0661294) q[1];
sx q[1];
rz(-2.9758458) q[1];
rz(1.3158348) q[2];
sx q[2];
rz(-2.795965) q[2];
sx q[2];
rz(-0.74898455) q[2];
rz(1.6809088) q[3];
sx q[3];
rz(-0.80169818) q[3];
sx q[3];
rz(2.0924951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
