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
rz(1.1295553) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0985467) q[0];
sx q[0];
rz(-0.73972964) q[0];
sx q[0];
rz(-1.3003527) q[0];
rz(-0.66313498) q[2];
sx q[2];
rz(-1.5485933) q[2];
sx q[2];
rz(-0.93142366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1892231) q[1];
sx q[1];
rz(-0.86878759) q[1];
sx q[1];
rz(-2.7527806) q[1];
x q[2];
rz(-2.8087691) q[3];
sx q[3];
rz(-1.5378008) q[3];
sx q[3];
rz(2.8590315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1019885) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(2.2339036) q[2];
rz(3.0229819) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-0.25092009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(-0.10435852) q[0];
rz(-1.1907499) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(-2.6844535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16743827) q[0];
sx q[0];
rz(-1.7642685) q[0];
sx q[0];
rz(-0.03137145) q[0];
rz(-pi) q[1];
rz(-2.6847849) q[2];
sx q[2];
rz(-2.1965082) q[2];
sx q[2];
rz(-0.63969757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7434843) q[1];
sx q[1];
rz(-0.32998514) q[1];
sx q[1];
rz(-1.2742001) q[1];
x q[2];
rz(3.0035104) q[3];
sx q[3];
rz(-2.2347817) q[3];
sx q[3];
rz(-2.9241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5856058) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-0.15698329) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(-3.1275911) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3954725) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-2.795862) q[0];
rz(1.455447) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.570943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3368944) q[0];
sx q[0];
rz(-1.0270627) q[0];
sx q[0];
rz(3.0377484) q[0];
rz(2.1934671) q[2];
sx q[2];
rz(-1.5678763) q[2];
sx q[2];
rz(-0.44026431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0367056) q[1];
sx q[1];
rz(-1.3574523) q[1];
sx q[1];
rz(-1.5206737) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9390252) q[3];
sx q[3];
rz(-0.93155471) q[3];
sx q[3];
rz(0.57697421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7741144) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(-0.95743123) q[2];
rz(-2.0188324) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(1.7694337) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51266176) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(1.7515901) q[1];
sx q[1];
rz(-1.8252239) q[1];
sx q[1];
rz(-1.4768627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6679803) q[0];
sx q[0];
rz(-1.2897217) q[0];
sx q[0];
rz(-0.54664183) q[0];
rz(1.8091168) q[2];
sx q[2];
rz(-1.6523978) q[2];
sx q[2];
rz(1.2343182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37472607) q[1];
sx q[1];
rz(-2.0927168) q[1];
sx q[1];
rz(0.38066997) q[1];
rz(-pi) q[2];
rz(1.2943997) q[3];
sx q[3];
rz(-1.2528462) q[3];
sx q[3];
rz(-2.2459787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5433898) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(0.99679917) q[2];
rz(1.2654842) q[3];
sx q[3];
rz(-2.4235453) q[3];
sx q[3];
rz(-0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0604414) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(1.0720217) q[0];
rz(-2.187166) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(1.4929474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77253262) q[0];
sx q[0];
rz(-0.81789069) q[0];
sx q[0];
rz(-1.8408804) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0132332) q[2];
sx q[2];
rz(-2.2432598) q[2];
sx q[2];
rz(-0.44328025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.055549057) q[1];
sx q[1];
rz(-1.2437268) q[1];
sx q[1];
rz(-1.7213166) q[1];
rz(0.29376438) q[3];
sx q[3];
rz(-2.8076594) q[3];
sx q[3];
rz(-1.5042083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91168779) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(-1.9913199) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(0.40157792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7683485) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.9819697) q[0];
rz(1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.9083091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984959) q[0];
sx q[0];
rz(-2.2944106) q[0];
sx q[0];
rz(-1.1529403) q[0];
x q[1];
rz(3.0843098) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.48207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4602051) q[1];
sx q[1];
rz(-2.3355995) q[1];
sx q[1];
rz(-2.5780923) q[1];
rz(-pi) q[2];
rz(-2.7090577) q[3];
sx q[3];
rz(-0.56946856) q[3];
sx q[3];
rz(1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.636574) q[2];
sx q[2];
rz(-2.7754112) q[2];
sx q[2];
rz(-1.150307) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.4897414) q[2];
sx q[2];
rz(1.3339804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4515069) q[1];
sx q[1];
rz(-1.2770318) q[1];
sx q[1];
rz(-0.34505941) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8774639) q[3];
sx q[3];
rz(-1.0326516) q[3];
sx q[3];
rz(1.9028185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(-1.1807582) q[2];
rz(2.0598038) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54843724) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(-2.4539808) q[0];
rz(1.9718735) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(-2.7489472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8075645) q[0];
sx q[0];
rz(-1.0962209) q[0];
sx q[0];
rz(-0.15888283) q[0];
rz(-0.41021014) q[2];
sx q[2];
rz(-2.8034119) q[2];
sx q[2];
rz(-1.2615801) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.709014) q[1];
sx q[1];
rz(-0.48332941) q[1];
sx q[1];
rz(0.18928115) q[1];
x q[2];
rz(-0.4057377) q[3];
sx q[3];
rz(-2.9364481) q[3];
sx q[3];
rz(-2.1577378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9795867) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(1.4957734) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(-0.60105598) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2163806) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(-1.8076757) q[0];
rz(1.9981617) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(-0.7935895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5549907) q[0];
sx q[0];
rz(-1.412801) q[0];
sx q[0];
rz(-1.642307) q[0];
rz(0.066448224) q[2];
sx q[2];
rz(-2.1796103) q[2];
sx q[2];
rz(2.7087536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55793437) q[1];
sx q[1];
rz(-1.0685295) q[1];
sx q[1];
rz(-1.7947547) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7904441) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.3212412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44635832) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(-0.71600437) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(-2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.42651549) q[0];
sx q[0];
rz(-0.85549131) q[0];
sx q[0];
rz(-2.9606384) q[0];
rz(-1.9521693) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-0.98035556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086209379) q[0];
sx q[0];
rz(-2.6329941) q[0];
sx q[0];
rz(0.17091708) q[0];
x q[1];
rz(0.75677121) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(0.41220081) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36114039) q[1];
sx q[1];
rz(-1.3899511) q[1];
sx q[1];
rz(-0.7493345) q[1];
rz(-pi) q[2];
rz(1.7813563) q[3];
sx q[3];
rz(-1.9580152) q[3];
sx q[3];
rz(-2.5507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.02562) q[2];
sx q[2];
rz(-2.9774234) q[2];
sx q[2];
rz(-1.3244965) q[2];
rz(0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(-3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7272335) q[0];
sx q[0];
rz(-1.4215195) q[0];
sx q[0];
rz(-0.97378578) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(-1.235511) q[2];
sx q[2];
rz(-1.4852471) q[2];
sx q[2];
rz(0.5813364) q[2];
rz(-0.11304819) q[3];
sx q[3];
rz(-0.77533508) q[3];
sx q[3];
rz(-1.2067457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
