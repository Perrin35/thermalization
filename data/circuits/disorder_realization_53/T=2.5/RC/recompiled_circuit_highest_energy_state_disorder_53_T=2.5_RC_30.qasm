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
rz(-1.1797969) q[0];
sx q[0];
rz(-1.0581764) q[0];
sx q[0];
rz(-0.41184586) q[0];
rz(0.62077275) q[1];
sx q[1];
rz(-0.70561916) q[1];
sx q[1];
rz(-0.128428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0043853) q[0];
sx q[0];
rz(-2.5126728) q[0];
sx q[0];
rz(-0.0034751277) q[0];
x q[1];
rz(-1.0780222) q[2];
sx q[2];
rz(-1.2076305) q[2];
sx q[2];
rz(2.9271379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.4219749) q[1];
sx q[1];
rz(-2.1901844) q[1];
sx q[1];
rz(-1.8108111) q[1];
rz(-pi) q[2];
rz(3.0459636) q[3];
sx q[3];
rz(-1.9668035) q[3];
sx q[3];
rz(-1.2850266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8992369) q[2];
sx q[2];
rz(-0.54348102) q[2];
sx q[2];
rz(-0.71329722) q[2];
rz(1.9244309) q[3];
sx q[3];
rz(-1.4274012) q[3];
sx q[3];
rz(-3.0349019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941968) q[0];
sx q[0];
rz(-2.1155745) q[0];
sx q[0];
rz(0.062653616) q[0];
rz(-2.2368536) q[1];
sx q[1];
rz(-0.86974564) q[1];
sx q[1];
rz(-3.0135801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63120121) q[0];
sx q[0];
rz(-1.2423853) q[0];
sx q[0];
rz(-1.2174219) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4887183) q[2];
sx q[2];
rz(-1.0316398) q[2];
sx q[2];
rz(1.7918794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3941895) q[1];
sx q[1];
rz(-0.48626562) q[1];
sx q[1];
rz(-3.1267605) q[1];
rz(-pi) q[2];
rz(-0.73776695) q[3];
sx q[3];
rz(-1.4633661) q[3];
sx q[3];
rz(-0.28001912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5855828) q[2];
sx q[2];
rz(-1.9364919) q[2];
sx q[2];
rz(-1.0949868) q[2];
rz(-0.60947841) q[3];
sx q[3];
rz(-0.947781) q[3];
sx q[3];
rz(1.3013209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-0.13780093) q[0];
sx q[0];
rz(-2.4584558) q[0];
sx q[0];
rz(-1.1411427) q[0];
rz(-1.2992651) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(3.0019147) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3052854) q[0];
sx q[0];
rz(-1.7392312) q[0];
sx q[0];
rz(1.3091716) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9150117) q[2];
sx q[2];
rz(-1.7936631) q[2];
sx q[2];
rz(1.8755166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91979474) q[1];
sx q[1];
rz(-1.1573024) q[1];
sx q[1];
rz(-1.7899124) q[1];
rz(-2.1909149) q[3];
sx q[3];
rz(-1.0178627) q[3];
sx q[3];
rz(1.3554791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37207681) q[2];
sx q[2];
rz(-2.3774827) q[2];
sx q[2];
rz(-2.8238943) q[2];
rz(2.5959173) q[3];
sx q[3];
rz(-0.92394865) q[3];
sx q[3];
rz(-1.8359449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9724156) q[0];
sx q[0];
rz(-1.7298052) q[0];
sx q[0];
rz(-2.8929945) q[0];
rz(1.2480805) q[1];
sx q[1];
rz(-1.0634407) q[1];
sx q[1];
rz(0.6691106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4658858) q[0];
sx q[0];
rz(-2.0629333) q[0];
sx q[0];
rz(-1.6791315) q[0];
rz(-pi) q[1];
rz(-1.8774154) q[2];
sx q[2];
rz(-0.82984041) q[2];
sx q[2];
rz(-1.6774943) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65655639) q[1];
sx q[1];
rz(-1.4244231) q[1];
sx q[1];
rz(1.778641) q[1];
x q[2];
rz(0.44336278) q[3];
sx q[3];
rz(-2.0896455) q[3];
sx q[3];
rz(-0.31135294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.756838) q[2];
sx q[2];
rz(-2.8110795) q[2];
sx q[2];
rz(-3.030704) q[2];
rz(0.97892654) q[3];
sx q[3];
rz(-1.3803692) q[3];
sx q[3];
rz(-3.0389649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0454309) q[0];
sx q[0];
rz(-1.0935421) q[0];
sx q[0];
rz(0.1732711) q[0];
rz(0.551956) q[1];
sx q[1];
rz(-2.5872784) q[1];
sx q[1];
rz(-0.90131235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1364897) q[0];
sx q[0];
rz(-2.1629961) q[0];
sx q[0];
rz(-0.96969684) q[0];
rz(-pi) q[1];
rz(-2.1036673) q[2];
sx q[2];
rz(-2.0882389) q[2];
sx q[2];
rz(-2.8252476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3696987) q[1];
sx q[1];
rz(-1.7189184) q[1];
sx q[1];
rz(3.1225608) q[1];
x q[2];
rz(-0.41660033) q[3];
sx q[3];
rz(-1.9754656) q[3];
sx q[3];
rz(-2.4373151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1426455) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(2.9360845) q[2];
rz(2.5213126) q[3];
sx q[3];
rz(-0.55448237) q[3];
sx q[3];
rz(-1.2834572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.810629) q[0];
sx q[0];
rz(-2.0331419) q[0];
sx q[0];
rz(-0.95322815) q[0];
rz(0.78298059) q[1];
sx q[1];
rz(-0.51376659) q[1];
sx q[1];
rz(-0.89797529) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1785068) q[0];
sx q[0];
rz(-2.0971812) q[0];
sx q[0];
rz(0.25776074) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9421517) q[2];
sx q[2];
rz(-1.1716649) q[2];
sx q[2];
rz(-2.6836306) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0197438) q[1];
sx q[1];
rz(-1.6770096) q[1];
sx q[1];
rz(-1.386045) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98358752) q[3];
sx q[3];
rz(-2.7977571) q[3];
sx q[3];
rz(1.6621906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7440685) q[2];
sx q[2];
rz(-2.4351138) q[2];
sx q[2];
rz(1.2495329) q[2];
rz(-1.9254098) q[3];
sx q[3];
rz(-1.3666697) q[3];
sx q[3];
rz(1.7887615) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3406496) q[0];
sx q[0];
rz(-0.018434374) q[0];
sx q[0];
rz(1.8311485) q[0];
rz(-0.75366968) q[1];
sx q[1];
rz(-0.33652702) q[1];
sx q[1];
rz(2.9897142) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57275332) q[0];
sx q[0];
rz(-1.4985159) q[0];
sx q[0];
rz(-0.092503017) q[0];
rz(-2.8114258) q[2];
sx q[2];
rz(-1.6474198) q[2];
sx q[2];
rz(-0.1506131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2703524) q[1];
sx q[1];
rz(-2.0899822) q[1];
sx q[1];
rz(-3.0445208) q[1];
rz(-2.2521654) q[3];
sx q[3];
rz(-1.424289) q[3];
sx q[3];
rz(1.5582939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.001363) q[2];
sx q[2];
rz(-1.9915853) q[2];
sx q[2];
rz(-1.4353732) q[2];
rz(2.5367149) q[3];
sx q[3];
rz(-1.4698272) q[3];
sx q[3];
rz(0.83355561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75989938) q[0];
sx q[0];
rz(-0.23867358) q[0];
sx q[0];
rz(0.47644404) q[0];
rz(2.6051688) q[1];
sx q[1];
rz(-1.8073578) q[1];
sx q[1];
rz(0.36472067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3140253) q[0];
sx q[0];
rz(-0.3985346) q[0];
sx q[0];
rz(0.78273313) q[0];
rz(-pi) q[1];
rz(0.1171072) q[2];
sx q[2];
rz(-2.6016781) q[2];
sx q[2];
rz(1.2123666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9832416) q[1];
sx q[1];
rz(-2.4102306) q[1];
sx q[1];
rz(-1.8476827) q[1];
rz(-pi) q[2];
rz(-2.164806) q[3];
sx q[3];
rz(-1.6315809) q[3];
sx q[3];
rz(1.2244299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81858188) q[2];
sx q[2];
rz(-1.8135169) q[2];
sx q[2];
rz(-2.2200457) q[2];
rz(1.796683) q[3];
sx q[3];
rz(-1.4321046) q[3];
sx q[3];
rz(0.014605453) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0607249) q[0];
sx q[0];
rz(-2.9089071) q[0];
sx q[0];
rz(0.099207148) q[0];
rz(-1.5946782) q[1];
sx q[1];
rz(-0.64259905) q[1];
sx q[1];
rz(3.0329472) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54421595) q[0];
sx q[0];
rz(-1.6013494) q[0];
sx q[0];
rz(1.0511418) q[0];
x q[1];
rz(-1.0119087) q[2];
sx q[2];
rz(-2.4299653) q[2];
sx q[2];
rz(0.54838442) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80444944) q[1];
sx q[1];
rz(-1.2986309) q[1];
sx q[1];
rz(2.0963349) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9401764) q[3];
sx q[3];
rz(-1.7699852) q[3];
sx q[3];
rz(-1.3808113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8607311) q[2];
sx q[2];
rz(-2.9555369) q[2];
sx q[2];
rz(0.59530386) q[2];
rz(0.91296494) q[3];
sx q[3];
rz(-0.76415092) q[3];
sx q[3];
rz(1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9106134) q[0];
sx q[0];
rz(-1.9934886) q[0];
sx q[0];
rz(-0.30368152) q[0];
rz(-0.23823711) q[1];
sx q[1];
rz(-0.99681011) q[1];
sx q[1];
rz(2.4441267) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.781002) q[0];
sx q[0];
rz(-1.7502898) q[0];
sx q[0];
rz(0.54474564) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4749305) q[2];
sx q[2];
rz(-1.2156475) q[2];
sx q[2];
rz(-0.42527664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94983236) q[1];
sx q[1];
rz(-1.0737541) q[1];
sx q[1];
rz(-2.0326896) q[1];
rz(-1.0806709) q[3];
sx q[3];
rz(-2.625386) q[3];
sx q[3];
rz(1.4913477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4136228) q[2];
sx q[2];
rz(-2.0406849) q[2];
sx q[2];
rz(-3.0481763) q[2];
rz(1.5591239) q[3];
sx q[3];
rz(-3.0714572) q[3];
sx q[3];
rz(0.087433405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14417905) q[0];
sx q[0];
rz(-1.1638673) q[0];
sx q[0];
rz(-1.5735004) q[0];
rz(-2.6724124) q[1];
sx q[1];
rz(-2.4041685) q[1];
sx q[1];
rz(2.266177) q[1];
rz(-0.45892528) q[2];
sx q[2];
rz(-2.060072) q[2];
sx q[2];
rz(1.2960808) q[2];
rz(-0.31147891) q[3];
sx q[3];
rz(-2.3337915) q[3];
sx q[3];
rz(1.3561377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
