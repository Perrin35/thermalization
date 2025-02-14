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
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(-2.0120373) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72973824) q[0];
sx q[0];
rz(-1.751873) q[0];
sx q[0];
rz(0.84946655) q[0];
x q[1];
rz(3.1055345) q[2];
sx q[2];
rz(-2.4781422) q[2];
sx q[2];
rz(-0.61095881) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3870377) q[1];
sx q[1];
rz(-2.3554152) q[1];
sx q[1];
rz(-1.1494067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10068746) q[3];
sx q[3];
rz(-0.33439454) q[3];
sx q[3];
rz(-1.3833801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0396042) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(-3.0229819) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(0.25092009) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94659656) q[0];
sx q[0];
rz(-0.66221607) q[0];
sx q[0];
rz(-3.0372341) q[0];
rz(1.9508427) q[1];
sx q[1];
rz(-1.933681) q[1];
sx q[1];
rz(2.6844535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3973244) q[0];
sx q[0];
rz(-1.6015823) q[0];
sx q[0];
rz(1.7643614) q[0];
rz(-2.6847849) q[2];
sx q[2];
rz(-0.94508445) q[2];
sx q[2];
rz(-2.5018951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7105661) q[1];
sx q[1];
rz(-1.2557286) q[1];
sx q[1];
rz(3.0418212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7449154) q[3];
sx q[3];
rz(-2.4655364) q[3];
sx q[3];
rz(3.1372748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55598688) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-2.9846094) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3954725) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(2.795862) q[0];
rz(1.455447) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(-1.5706496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80469823) q[0];
sx q[0];
rz(-2.11453) q[0];
sx q[0];
rz(3.0377484) q[0];
x q[1];
rz(-1.5758031) q[2];
sx q[2];
rz(-2.5189159) q[2];
sx q[2];
rz(1.1345991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8040672) q[1];
sx q[1];
rz(-0.21906549) q[1];
sx q[1];
rz(0.22727025) q[1];
x q[2];
rz(1.3064874) q[3];
sx q[3];
rz(-0.66625957) q[3];
sx q[3];
rz(-2.8961757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3674783) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(0.95743123) q[2];
rz(1.1227603) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(1.3721589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51266176) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5976394) q[0];
rz(-1.7515901) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(1.66473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2643971) q[0];
sx q[0];
rz(-1.0478643) q[0];
sx q[0];
rz(-1.8967129) q[0];
rz(-pi) q[1];
rz(1.3324758) q[2];
sx q[2];
rz(-1.4891948) q[2];
sx q[2];
rz(1.2343182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0904729) q[1];
sx q[1];
rz(-0.63544151) q[1];
sx q[1];
rz(-0.99721554) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2943997) q[3];
sx q[3];
rz(-1.8887465) q[3];
sx q[3];
rz(-0.89561392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59820286) q[2];
sx q[2];
rz(-1.2849839) q[2];
sx q[2];
rz(-0.99679917) q[2];
rz(-1.2654842) q[3];
sx q[3];
rz(-2.4235453) q[3];
sx q[3];
rz(0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.0811512) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(-2.0695709) q[0];
rz(-2.187166) q[1];
sx q[1];
rz(-1.9642893) q[1];
sx q[1];
rz(-1.4929474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98534855) q[0];
sx q[0];
rz(-1.3748504) q[0];
sx q[0];
rz(2.3702456) q[0];
rz(0.72228186) q[2];
sx q[2];
rz(-1.2292635) q[2];
sx q[2];
rz(2.3010437) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6449738) q[1];
sx q[1];
rz(-2.7826834) q[1];
sx q[1];
rz(0.41618698) q[1];
x q[2];
rz(-2.8478283) q[3];
sx q[3];
rz(-0.33393327) q[3];
sx q[3];
rz(1.5042083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2299049) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(0.041042717) q[3];
sx q[3];
rz(-1.5464455) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.9819697) q[0];
rz(1.6805964) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.2332835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4077743) q[0];
sx q[0];
rz(-2.3253157) q[0];
sx q[0];
rz(-2.7110148) q[0];
rz(-0.057282863) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.48207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4602051) q[1];
sx q[1];
rz(-2.3355995) q[1];
sx q[1];
rz(0.56350033) q[1];
x q[2];
rz(-2.7090577) q[3];
sx q[3];
rz(-0.56946856) q[3];
sx q[3];
rz(-1.3506952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50501862) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(1.9912857) q[2];
rz(2.4023174) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(-2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4742541) q[0];
sx q[0];
rz(-0.69304729) q[0];
sx q[0];
rz(-0.44878238) q[0];
rz(1.5397286) q[1];
sx q[1];
rz(-2.0346784) q[1];
sx q[1];
rz(1.1610228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13793531) q[0];
sx q[0];
rz(-1.9788392) q[0];
sx q[0];
rz(-1.8613226) q[0];
rz(-0.77552253) q[2];
sx q[2];
rz(-0.65185129) q[2];
sx q[2];
rz(1.8076123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55884493) q[1];
sx q[1];
rz(-0.44932355) q[1];
sx q[1];
rz(0.72968633) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26412873) q[3];
sx q[3];
rz(-2.1089411) q[3];
sx q[3];
rz(1.2387741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9728969) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(1.9608344) q[2];
rz(-2.0598038) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(-2.7150174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54843724) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(1.1697191) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(-2.7489472) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8317345) q[0];
sx q[0];
rz(-1.7119954) q[0];
sx q[0];
rz(-2.0505428) q[0];
x q[1];
rz(-0.41021014) q[2];
sx q[2];
rz(-0.33818076) q[2];
sx q[2];
rz(-1.8800125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4959405) q[1];
sx q[1];
rz(-2.0447746) q[1];
sx q[1];
rz(-1.4723634) q[1];
rz(-1.6527376) q[3];
sx q[3];
rz(-1.382516) q[3];
sx q[3];
rz(0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9795867) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(-1.5227854) q[3];
sx q[3];
rz(-0.77097547) q[3];
sx q[3];
rz(-2.5405367) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2163806) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(-1.8076757) q[0];
rz(-1.1434309) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(-2.3480031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.586602) q[0];
sx q[0];
rz(-1.7287917) q[0];
sx q[0];
rz(1.4992856) q[0];
x q[1];
rz(-0.96094544) q[2];
sx q[2];
rz(-1.5163002) q[2];
sx q[2];
rz(1.1759963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55793437) q[1];
sx q[1];
rz(-1.0685295) q[1];
sx q[1];
rz(-1.7947547) q[1];
x q[2];
rz(-0.10754866) q[3];
sx q[3];
rz(-2.0189773) q[3];
sx q[3];
rz(1.5654237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6952343) q[2];
sx q[2];
rz(-1.0506722) q[2];
sx q[2];
rz(2.0900772) q[2];
rz(0.71600437) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-0.85549131) q[0];
sx q[0];
rz(-2.9606384) q[0];
rz(-1.9521693) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(2.1612371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0553833) q[0];
sx q[0];
rz(-0.50859857) q[0];
sx q[0];
rz(2.9706756) q[0];
rz(-2.3848214) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(0.41220081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3754543) q[1];
sx q[1];
rz(-0.836529) q[1];
sx q[1];
rz(-1.8155273) q[1];
rz(-2.7465024) q[3];
sx q[3];
rz(-1.3760341) q[3];
sx q[3];
rz(2.2421746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1159726) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.3244965) q[2];
rz(-2.4888511) q[3];
sx q[3];
rz(-0.96941152) q[3];
sx q[3];
rz(3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7272335) q[0];
sx q[0];
rz(-1.4215195) q[0];
sx q[0];
rz(-0.97378578) q[0];
rz(-0.7242135) q[1];
sx q[1];
rz(-2.0661294) q[1];
sx q[1];
rz(-2.9758458) q[1];
rz(0.090567055) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(0.77213415) q[3];
sx q[3];
rz(-1.6498389) q[3];
sx q[3];
rz(-2.6966358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
