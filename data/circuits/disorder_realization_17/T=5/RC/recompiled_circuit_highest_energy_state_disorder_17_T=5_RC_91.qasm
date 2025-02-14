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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(-2.4315779) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(1.8097872) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074866991) q[0];
sx q[0];
rz(-0.35860379) q[0];
sx q[0];
rz(-0.68645607) q[0];
x q[1];
rz(-1.3178145) q[2];
sx q[2];
rz(-1.2377146) q[2];
sx q[2];
rz(-2.2673504) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.097447473) q[1];
sx q[1];
rz(-0.65939476) q[1];
sx q[1];
rz(1.161518) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.136257) q[3];
sx q[3];
rz(-2.4255803) q[3];
sx q[3];
rz(0.93040066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72267246) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(-2.326272) q[2];
rz(2.3319862) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(2.1965006) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5792712) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54540578) q[0];
sx q[0];
rz(-0.89841398) q[0];
sx q[0];
rz(0.59242146) q[0];
x q[1];
rz(0.42065545) q[2];
sx q[2];
rz(-1.9845982) q[2];
sx q[2];
rz(0.54383531) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89797276) q[1];
sx q[1];
rz(-1.0884411) q[1];
sx q[1];
rz(1.6887168) q[1];
rz(-pi) q[2];
rz(0.48247997) q[3];
sx q[3];
rz(-2.2594995) q[3];
sx q[3];
rz(2.3107993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-2.6914524) q[2];
rz(-2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-2.1639737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728773) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(-1.3866792) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(0.6001572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068014) q[0];
sx q[0];
rz(-1.8610483) q[0];
sx q[0];
rz(-1.6192379) q[0];
x q[1];
rz(-1.6515031) q[2];
sx q[2];
rz(-2.5332402) q[2];
sx q[2];
rz(-2.0474912) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0827209) q[1];
sx q[1];
rz(-1.8959672) q[1];
sx q[1];
rz(0.60589686) q[1];
rz(-pi) q[2];
rz(-2.8801877) q[3];
sx q[3];
rz(-2.5525744) q[3];
sx q[3];
rz(-2.1847069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0967789) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(-2.0951994) q[2];
rz(-0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1716877) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(2.0890253) q[0];
rz(-1.9453847) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(-0.41950163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7561121) q[0];
sx q[0];
rz(-1.0051553) q[0];
sx q[0];
rz(-1.2570981) q[0];
x q[1];
rz(0.93650903) q[2];
sx q[2];
rz(-1.7542931) q[2];
sx q[2];
rz(1.2681792) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0191325) q[1];
sx q[1];
rz(-1.9472709) q[1];
sx q[1];
rz(-2.2458162) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6330248) q[3];
sx q[3];
rz(-1.8265822) q[3];
sx q[3];
rz(-3.0884107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(1.0901964) q[2];
rz(-2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5721385) q[0];
sx q[0];
rz(1.74362) q[0];
rz(-3.0086503) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(0.001210777) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0069790445) q[0];
sx q[0];
rz(-1.9513357) q[0];
sx q[0];
rz(-1.5534205) q[0];
rz(1.8563849) q[2];
sx q[2];
rz(-1.4175709) q[2];
sx q[2];
rz(2.0375348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9922793) q[1];
sx q[1];
rz(-2.1985801) q[1];
sx q[1];
rz(-1.1924442) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7453543) q[3];
sx q[3];
rz(-0.15197411) q[3];
sx q[3];
rz(-1.1308972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2289537) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(2.9310628) q[2];
rz(1.0362961) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.2150631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-2.2192661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2012537) q[0];
sx q[0];
rz(-2.8329599) q[0];
sx q[0];
rz(-1.2135394) q[0];
rz(-pi) q[1];
rz(0.90409235) q[2];
sx q[2];
rz(-1.4897926) q[2];
sx q[2];
rz(-1.9503649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74444425) q[1];
sx q[1];
rz(-1.6849298) q[1];
sx q[1];
rz(-0.95790095) q[1];
rz(-pi) q[2];
rz(-0.57862307) q[3];
sx q[3];
rz(-0.64217438) q[3];
sx q[3];
rz(0.35143055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-0.40194884) q[2];
rz(-1.2881783) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.2791546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61889082) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(3.0772305) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(1.3305957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42444077) q[0];
sx q[0];
rz(-1.0006051) q[0];
sx q[0];
rz(1.4217522) q[0];
rz(-0.66489545) q[2];
sx q[2];
rz(-2.4192296) q[2];
sx q[2];
rz(1.5255873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45695282) q[1];
sx q[1];
rz(-1.7906396) q[1];
sx q[1];
rz(2.6067642) q[1];
rz(-pi) q[2];
rz(-0.17223151) q[3];
sx q[3];
rz(-2.1759998) q[3];
sx q[3];
rz(-1.8021405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4055206) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-0.68354496) q[2];
rz(2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(-2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(1.1731359) q[0];
rz(0.37551156) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(-2.1048022) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3507337) q[0];
sx q[0];
rz(-3.1111801) q[0];
sx q[0];
rz(0.68002547) q[0];
rz(2.9762245) q[2];
sx q[2];
rz(-1.8435119) q[2];
sx q[2];
rz(2.337526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97809659) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(2.3771493) q[1];
rz(-pi) q[2];
rz(1.0408786) q[3];
sx q[3];
rz(-1.138759) q[3];
sx q[3];
rz(1.0505067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5658687) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(0.73053989) q[2];
rz(2.9324487) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(-1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65649477) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(3.094161) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(2.2304631) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1909433) q[0];
sx q[0];
rz(-1.5438269) q[0];
sx q[0];
rz(3.1080676) q[0];
rz(0.92911559) q[2];
sx q[2];
rz(-1.5589899) q[2];
sx q[2];
rz(-0.75071834) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92420134) q[1];
sx q[1];
rz(-1.5581521) q[1];
sx q[1];
rz(-0.36751698) q[1];
rz(-pi) q[2];
rz(-2.9373475) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(1.1413407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(0.16858777) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-2.2897661) q[3];
sx q[3];
rz(-1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.99437) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(2.6841573) q[0];
rz(1.9153197) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(1.7785243) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03662388) q[0];
sx q[0];
rz(-1.5104806) q[0];
sx q[0];
rz(-0.6078267) q[0];
x q[1];
rz(-0.3293475) q[2];
sx q[2];
rz(-0.65635704) q[2];
sx q[2];
rz(-2.658297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0911965) q[1];
sx q[1];
rz(-1.1938057) q[1];
sx q[1];
rz(-1.7524377) q[1];
rz(-pi) q[2];
rz(1.4665589) q[3];
sx q[3];
rz(-0.6851894) q[3];
sx q[3];
rz(-0.70885056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3739796) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(0.45200959) q[2];
rz(0.40397817) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6954738) q[0];
sx q[0];
rz(-1.6071381) q[0];
sx q[0];
rz(0.18679609) q[0];
rz(-0.52234621) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(0.28235565) q[2];
sx q[2];
rz(-1.3236125) q[2];
sx q[2];
rz(-0.57409928) q[2];
rz(-0.50857827) q[3];
sx q[3];
rz(-2.8420035) q[3];
sx q[3];
rz(-3.0723078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
