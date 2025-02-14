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
rz(-0.71159166) q[0];
sx q[0];
rz(-2.0366259) q[0];
sx q[0];
rz(-2.0018863) q[0];
rz(-0.79144129) q[1];
sx q[1];
rz(-1.0410407) q[1];
sx q[1];
rz(2.0120373) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.043046) q[0];
sx q[0];
rz(-0.73972964) q[0];
sx q[0];
rz(1.3003527) q[0];
rz(-pi) q[1];
rz(-0.66313498) q[2];
sx q[2];
rz(-1.5929993) q[2];
sx q[2];
rz(0.93142366) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.754555) q[1];
sx q[1];
rz(-2.3554152) q[1];
sx q[1];
rz(1.992186) q[1];
rz(-0.10068746) q[3];
sx q[3];
rz(-0.33439454) q[3];
sx q[3];
rz(1.7582126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0396042) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(-2.2339036) q[2];
rz(-3.0229819) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-0.66221607) q[0];
sx q[0];
rz(0.10435852) q[0];
rz(1.9508427) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(0.45713919) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9741544) q[0];
sx q[0];
rz(-1.7642685) q[0];
sx q[0];
rz(-3.1102212) q[0];
rz(-pi) q[1];
rz(1.0227343) q[2];
sx q[2];
rz(-2.3853758) q[2];
sx q[2];
rz(0.05847419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7434843) q[1];
sx q[1];
rz(-2.8116075) q[1];
sx q[1];
rz(-1.8673926) q[1];
rz(0.13808226) q[3];
sx q[3];
rz(-0.90681091) q[3];
sx q[3];
rz(0.21747227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(-0.15698329) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(0.014001525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(0.34573063) q[0];
rz(1.455447) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(-1.5706496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5356876) q[0];
sx q[0];
rz(-0.55258026) q[0];
sx q[0];
rz(-1.4009757) q[0];
rz(-pi) q[1];
rz(-1.5657896) q[2];
sx q[2];
rz(-2.5189159) q[2];
sx q[2];
rz(-1.1345991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6863043) q[1];
sx q[1];
rz(-1.6197816) q[1];
sx q[1];
rz(-0.21360417) q[1];
x q[2];
rz(-0.9216347) q[3];
sx q[3];
rz(-1.4086257) q[3];
sx q[3];
rz(-2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7741144) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(-0.95743123) q[2];
rz(-1.1227603) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(-1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51266176) q[0];
sx q[0];
rz(-1.7191732) q[0];
sx q[0];
rz(1.5439532) q[0];
rz(1.7515901) q[1];
sx q[1];
rz(-1.8252239) q[1];
sx q[1];
rz(-1.4768627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8771956) q[0];
sx q[0];
rz(-2.0937284) q[0];
sx q[0];
rz(-1.2448798) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0576287) q[2];
sx q[2];
rz(-1.8083085) q[2];
sx q[2];
rz(-0.31667865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99914306) q[1];
sx q[1];
rz(-1.2428741) q[1];
sx q[1];
rz(1.0161922) q[1];
rz(-0.32962004) q[3];
sx q[3];
rz(-1.8330036) q[3];
sx q[3];
rz(2.5548558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5433898) q[2];
sx q[2];
rz(-1.2849839) q[2];
sx q[2];
rz(-0.99679917) q[2];
rz(1.2654842) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0811512) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(-2.0695709) q[0];
rz(2.187166) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(-1.4929474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98534855) q[0];
sx q[0];
rz(-1.3748504) q[0];
sx q[0];
rz(-0.77134706) q[0];
rz(-pi) q[1];
rz(-0.72228186) q[2];
sx q[2];
rz(-1.2292635) q[2];
sx q[2];
rz(0.84054899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.055549057) q[1];
sx q[1];
rz(-1.2437268) q[1];
sx q[1];
rz(-1.7213166) q[1];
x q[2];
rz(0.32060726) q[3];
sx q[3];
rz(-1.6658452) q[3];
sx q[3];
rz(0.34497258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91168779) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(-1.9913199) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(-1.9819697) q[0];
rz(1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.9083091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.999812) q[0];
sx q[0];
rz(-1.8798057) q[0];
sx q[0];
rz(-2.3731493) q[0];
x q[1];
rz(-0.057282863) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.48207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.94168866) q[1];
sx q[1];
rz(-2.2268128) q[1];
sx q[1];
rz(1.0628878) q[1];
rz(-0.52652518) q[3];
sx q[3];
rz(-1.3428146) q[3];
sx q[3];
rz(0.15074061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50501862) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(-1.150307) q[2];
rz(0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(-0.44481835) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66733852) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(-2.6928103) q[0];
rz(1.601864) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(1.1610228) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338116) q[0];
sx q[0];
rz(-2.6454661) q[0];
sx q[0];
rz(-0.58519848) q[0];
rz(-pi) q[1];
x q[1];
rz(2.061474) q[2];
sx q[2];
rz(-1.1227692) q[2];
sx q[2];
rz(-0.44448822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6900857) q[1];
sx q[1];
rz(-1.2770318) q[1];
sx q[1];
rz(2.7965332) q[1];
rz(1.0169898) q[3];
sx q[3];
rz(-1.3447059) q[3];
sx q[3];
rz(-0.46976058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(1.1807582) q[2];
rz(-2.0598038) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(-0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5931554) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(-0.68761188) q[0];
rz(-1.9718735) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(0.39264548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3340281) q[0];
sx q[0];
rz(-1.0962209) q[0];
sx q[0];
rz(0.15888283) q[0];
x q[1];
rz(-2.8296109) q[2];
sx q[2];
rz(-1.7034966) q[2];
sx q[2];
rz(-3.0614982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029812977) q[1];
sx q[1];
rz(-1.4832442) q[1];
sx q[1];
rz(0.47595166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4057377) q[3];
sx q[3];
rz(-0.20514454) q[3];
sx q[3];
rz(2.1577378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16200599) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(-1.6188072) q[3];
sx q[3];
rz(-0.77097547) q[3];
sx q[3];
rz(2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9252121) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(1.3339169) q[0];
rz(1.1434309) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(2.3480031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1370572) q[0];
sx q[0];
rz(-1.5001778) q[0];
sx q[0];
rz(-2.9831992) q[0];
rz(-0.066448224) q[2];
sx q[2];
rz(-2.1796103) q[2];
sx q[2];
rz(0.43283909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55793437) q[1];
sx q[1];
rz(-2.0730632) q[1];
sx q[1];
rz(1.3468379) q[1];
rz(-pi) q[2];
rz(-0.10754866) q[3];
sx q[3];
rz(-2.0189773) q[3];
sx q[3];
rz(-1.5761689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6952343) q[2];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-0.85549131) q[0];
sx q[0];
rz(0.18095428) q[0];
rz(-1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-2.1612371) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086209379) q[0];
sx q[0];
rz(-2.6329941) q[0];
sx q[0];
rz(-0.17091708) q[0];
rz(2.0043401) q[2];
sx q[2];
rz(-1.1522276) q[2];
sx q[2];
rz(0.43805447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1229062) q[1];
sx q[1];
rz(-0.76670206) q[1];
sx q[1];
rz(0.26224978) q[1];
rz(2.6679831) q[3];
sx q[3];
rz(-0.43821143) q[3];
sx q[3];
rz(2.0357063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.02562) q[2];
sx q[2];
rz(-2.9774234) q[2];
sx q[2];
rz(-1.3244965) q[2];
rz(-0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41435913) q[0];
sx q[0];
rz(-1.4215195) q[0];
sx q[0];
rz(-0.97378578) q[0];
rz(2.4173792) q[1];
sx q[1];
rz(-2.0661294) q[1];
sx q[1];
rz(-2.9758458) q[1];
rz(-3.0510256) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(0.11304819) q[3];
sx q[3];
rz(-2.3662576) q[3];
sx q[3];
rz(1.934847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
