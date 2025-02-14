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
rz(1.1049668) q[0];
sx q[0];
rz(8.2850716) q[0];
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(1.1295553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4576042) q[0];
sx q[0];
rz(-2.2778371) q[0];
sx q[0];
rz(2.9024505) q[0];
rz(-pi) q[1];
rz(-1.5426251) q[2];
sx q[2];
rz(-2.2337388) q[2];
sx q[2];
rz(2.4848795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9523695) q[1];
sx q[1];
rz(-0.86878759) q[1];
sx q[1];
rz(2.7527806) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33282354) q[3];
sx q[3];
rz(-1.6037919) q[3];
sx q[3];
rz(-2.8590315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0396042) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-2.2339036) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(0.25092009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-1.933681) q[1];
sx q[1];
rz(2.6844535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32923082) q[0];
sx q[0];
rz(-0.19596772) q[0];
sx q[0];
rz(-1.4120483) q[0];
rz(0.45680775) q[2];
sx q[2];
rz(-0.94508445) q[2];
sx q[2];
rz(0.63969757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7105661) q[1];
sx q[1];
rz(-1.2557286) q[1];
sx q[1];
rz(-0.099771413) q[1];
rz(-2.2394286) q[3];
sx q[3];
rz(-1.6794101) q[3];
sx q[3];
rz(1.4387552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55598688) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-2.9846094) q[2];
rz(-0.5213151) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(-0.014001525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3954725) q[0];
sx q[0];
rz(-0.59500256) q[0];
sx q[0];
rz(2.795862) q[0];
rz(-1.6861457) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.570943) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60590505) q[0];
sx q[0];
rz(-2.5890124) q[0];
sx q[0];
rz(1.7406169) q[0];
rz(-pi) q[1];
rz(3.137998) q[2];
sx q[2];
rz(-2.1934641) q[2];
sx q[2];
rz(-1.1284356) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3375255) q[1];
sx q[1];
rz(-0.21906549) q[1];
sx q[1];
rz(-2.9143224) q[1];
rz(-pi) q[2];
rz(1.3064874) q[3];
sx q[3];
rz(-2.4753331) q[3];
sx q[3];
rz(-0.245417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7741144) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(-0.95743123) q[2];
rz(-1.1227603) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51266176) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(1.3900025) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(-1.4768627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47361237) q[0];
sx q[0];
rz(-1.851871) q[0];
sx q[0];
rz(-2.5949508) q[0];
x q[1];
rz(-1.8091168) q[2];
sx q[2];
rz(-1.4891948) q[2];
sx q[2];
rz(-1.9072744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0511198) q[1];
sx q[1];
rz(-0.63544151) q[1];
sx q[1];
rz(0.99721554) q[1];
x q[2];
rz(-1.847193) q[3];
sx q[3];
rz(-1.8887465) q[3];
sx q[3];
rz(-0.89561392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0604414) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(-2.0695709) q[0];
rz(-2.187166) q[1];
sx q[1];
rz(-1.9642893) q[1];
sx q[1];
rz(-1.4929474) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1562441) q[0];
sx q[0];
rz(-1.3748504) q[0];
sx q[0];
rz(2.3702456) q[0];
rz(1.1283595) q[2];
sx q[2];
rz(-0.89833288) q[2];
sx q[2];
rz(-2.6983124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4665597) q[1];
sx q[1];
rz(-1.7132812) q[1];
sx q[1];
rz(0.33054466) q[1];
x q[2];
rz(0.32060726) q[3];
sx q[3];
rz(-1.4757475) q[3];
sx q[3];
rz(2.7966201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2299049) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(-1.9819697) q[0];
rz(1.6805964) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.2332835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430968) q[0];
sx q[0];
rz(-0.84718207) q[0];
sx q[0];
rz(1.1529403) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5173814) q[2];
sx q[2];
rz(-1.513595) q[2];
sx q[2];
rz(3.055923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94168866) q[1];
sx q[1];
rz(-2.2268128) q[1];
sx q[1];
rz(-1.0628878) q[1];
rz(0.43253492) q[3];
sx q[3];
rz(-0.56946856) q[3];
sx q[3];
rz(1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.636574) q[2];
sx q[2];
rz(-2.7754112) q[2];
sx q[2];
rz(1.9912857) q[2];
rz(-0.73927528) q[3];
sx q[3];
rz(-1.247568) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4742541) q[0];
sx q[0];
rz(-0.69304729) q[0];
sx q[0];
rz(2.6928103) q[0];
rz(1.601864) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(1.1610228) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036573) q[0];
sx q[0];
rz(-1.9788392) q[0];
sx q[0];
rz(-1.8613226) q[0];
rz(-pi) q[1];
rz(2.3660701) q[2];
sx q[2];
rz(-2.4897414) q[2];
sx q[2];
rz(1.3339804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22299448) q[1];
sx q[1];
rz(-1.9004993) q[1];
sx q[1];
rz(-1.8818284) q[1];
rz(0.26412873) q[3];
sx q[3];
rz(-1.0326516) q[3];
sx q[3];
rz(1.2387741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9728969) q[2];
sx q[2];
rz(-1.8099512) q[2];
sx q[2];
rz(1.1807582) q[2];
rz(-2.0598038) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54843724) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(-1.9718735) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(-0.39264548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30985818) q[0];
sx q[0];
rz(-1.7119954) q[0];
sx q[0];
rz(2.0505428) q[0];
rz(-pi) q[1];
rz(-1.7101426) q[2];
sx q[2];
rz(-1.2616488) q[2];
sx q[2];
rz(1.4480556) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4959405) q[1];
sx q[1];
rz(-1.0968181) q[1];
sx q[1];
rz(1.6692293) q[1];
rz(2.9526934) q[3];
sx q[3];
rz(-1.6512863) q[3];
sx q[3];
rz(-0.98505011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9795867) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(1.6458192) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2163806) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(1.3339169) q[0];
rz(-1.1434309) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(2.3480031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5549907) q[0];
sx q[0];
rz(-1.412801) q[0];
sx q[0];
rz(-1.642307) q[0];
x q[1];
rz(1.6657532) q[2];
sx q[2];
rz(-2.5296202) q[2];
sx q[2];
rz(-2.8245935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90364218) q[1];
sx q[1];
rz(-1.3748843) q[1];
sx q[1];
rz(0.51301051) q[1];
rz(-3.034044) q[3];
sx q[3];
rz(-1.1226153) q[3];
sx q[3];
rz(-1.5761689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44635832) q[2];
sx q[2];
rz(-1.0506722) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(-2.4255883) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(-0.21014617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(2.9606384) q[0];
rz(1.9521693) q[1];
sx q[1];
rz(-1.1782497) q[1];
sx q[1];
rz(2.1612371) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10889036) q[0];
sx q[0];
rz(-2.0712896) q[0];
sx q[0];
rz(1.6653401) q[0];
rz(-0.75677121) q[2];
sx q[2];
rz(-2.5483661) q[2];
sx q[2];
rz(0.41220081) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3754543) q[1];
sx q[1];
rz(-2.3050637) q[1];
sx q[1];
rz(1.8155273) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7813563) q[3];
sx q[3];
rz(-1.9580152) q[3];
sx q[3];
rz(2.5507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.02562) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(-1.8170961) q[2];
rz(2.4888511) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.090567055) q[2];
sx q[2];
rz(-1.236785) q[2];
sx q[2];
rz(-1.0192237) q[2];
rz(-1.4606838) q[3];
sx q[3];
rz(-0.80169818) q[3];
sx q[3];
rz(2.0924951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
