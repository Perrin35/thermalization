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
rz(1.9617957) q[0];
sx q[0];
rz(-2.0834162) q[0];
sx q[0];
rz(0.41184586) q[0];
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
rz(2.1415048) q[0];
sx q[0];
rz(-0.94188085) q[0];
sx q[0];
rz(-1.5682683) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0635705) q[2];
sx q[2];
rz(-1.9339621) q[2];
sx q[2];
rz(-2.9271379) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8516396) q[1];
sx q[1];
rz(-1.3760097) q[1];
sx q[1];
rz(-2.508393) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9684359) q[3];
sx q[3];
rz(-1.4825882) q[3];
sx q[3];
rz(-0.24878605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2423557) q[2];
sx q[2];
rz(-2.5981116) q[2];
sx q[2];
rz(-0.71329722) q[2];
rz(-1.9244309) q[3];
sx q[3];
rz(-1.7141914) q[3];
sx q[3];
rz(-3.0349019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2941968) q[0];
sx q[0];
rz(-2.1155745) q[0];
sx q[0];
rz(3.078939) q[0];
rz(2.2368536) q[1];
sx q[1];
rz(-0.86974564) q[1];
sx q[1];
rz(-0.12801257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5103914) q[0];
sx q[0];
rz(-1.8992074) q[0];
sx q[0];
rz(-1.9241707) q[0];
rz(-0.92526154) q[2];
sx q[2];
rz(-1.0224258) q[2];
sx q[2];
rz(-2.546371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.377413) q[1];
sx q[1];
rz(-1.0845889) q[1];
sx q[1];
rz(1.5629565) q[1];
x q[2];
rz(2.9826134) q[3];
sx q[3];
rz(-2.3975054) q[3];
sx q[3];
rz(-1.1733623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5855828) q[2];
sx q[2];
rz(-1.9364919) q[2];
sx q[2];
rz(2.0466059) q[2];
rz(-2.5321142) q[3];
sx q[3];
rz(-2.1938117) q[3];
sx q[3];
rz(1.3013209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.13780093) q[0];
sx q[0];
rz(-0.68313685) q[0];
sx q[0];
rz(1.1411427) q[0];
rz(-1.2992651) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(3.0019147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2938625) q[0];
sx q[0];
rz(-2.8314857) q[0];
sx q[0];
rz(-0.98921138) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3422826) q[2];
sx q[2];
rz(-1.3499134) q[2];
sx q[2];
rz(-0.2538089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5798333) q[1];
sx q[1];
rz(-1.7711825) q[1];
sx q[1];
rz(0.42247129) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4927718) q[3];
sx q[3];
rz(-1.0534957) q[3];
sx q[3];
rz(0.14347883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37207681) q[2];
sx q[2];
rz(-0.76411) q[2];
sx q[2];
rz(2.8238943) q[2];
rz(-0.54567537) q[3];
sx q[3];
rz(-2.217644) q[3];
sx q[3];
rz(1.8359449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9724156) q[0];
sx q[0];
rz(-1.4117874) q[0];
sx q[0];
rz(-0.2485982) q[0];
rz(1.8935122) q[1];
sx q[1];
rz(-2.0781519) q[1];
sx q[1];
rz(-2.4724821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6921223) q[0];
sx q[0];
rz(-0.50296268) q[0];
sx q[0];
rz(0.19900222) q[0];
rz(1.2641773) q[2];
sx q[2];
rz(-0.82984041) q[2];
sx q[2];
rz(-1.6774943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94498881) q[1];
sx q[1];
rz(-1.3652063) q[1];
sx q[1];
rz(2.9920472) q[1];
rz(-pi) q[2];
rz(-2.1345244) q[3];
sx q[3];
rz(-1.9525212) q[3];
sx q[3];
rz(1.4907211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.756838) q[2];
sx q[2];
rz(-2.8110795) q[2];
sx q[2];
rz(-0.11088863) q[2];
rz(-2.1626661) q[3];
sx q[3];
rz(-1.3803692) q[3];
sx q[3];
rz(-3.0389649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0454309) q[0];
sx q[0];
rz(-2.0480506) q[0];
sx q[0];
rz(0.1732711) q[0];
rz(0.551956) q[1];
sx q[1];
rz(-0.55431429) q[1];
sx q[1];
rz(-2.2402803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93126585) q[0];
sx q[0];
rz(-1.0823609) q[0];
sx q[0];
rz(2.4573321) q[0];
x q[1];
rz(-0.58392151) q[2];
sx q[2];
rz(-1.1135227) q[2];
sx q[2];
rz(1.6032796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3696987) q[1];
sx q[1];
rz(-1.4226743) q[1];
sx q[1];
rz(0.019031899) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1327733) q[3];
sx q[3];
rz(-1.1896648) q[3];
sx q[3];
rz(1.0390067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1426455) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(0.20550814) q[2];
rz(0.62028003) q[3];
sx q[3];
rz(-0.55448237) q[3];
sx q[3];
rz(-1.8581355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3309636) q[0];
sx q[0];
rz(-1.1084508) q[0];
sx q[0];
rz(-0.95322815) q[0];
rz(0.78298059) q[1];
sx q[1];
rz(-2.6278261) q[1];
sx q[1];
rz(-2.2436174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47602874) q[0];
sx q[0];
rz(-1.7930287) q[0];
sx q[0];
rz(2.1118947) q[0];
rz(2.7165604) q[2];
sx q[2];
rz(-1.2298744) q[2];
sx q[2];
rz(-1.2630315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6728357) q[1];
sx q[1];
rz(-1.3870981) q[1];
sx q[1];
rz(0.10803799) q[1];
rz(-pi) q[2];
rz(-0.98358752) q[3];
sx q[3];
rz(-2.7977571) q[3];
sx q[3];
rz(1.4794021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39752412) q[2];
sx q[2];
rz(-0.70647883) q[2];
sx q[2];
rz(1.8920598) q[2];
rz(-1.2161829) q[3];
sx q[3];
rz(-1.3666697) q[3];
sx q[3];
rz(1.3528311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3406496) q[0];
sx q[0];
rz(-3.1231583) q[0];
sx q[0];
rz(-1.8311485) q[0];
rz(-0.75366968) q[1];
sx q[1];
rz(-2.8050656) q[1];
sx q[1];
rz(-2.9897142) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57275332) q[0];
sx q[0];
rz(-1.4985159) q[0];
sx q[0];
rz(3.0490896) q[0];
x q[1];
rz(-1.6517761) q[2];
sx q[2];
rz(-1.2416349) q[2];
sx q[2];
rz(1.6951814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8712403) q[1];
sx q[1];
rz(-2.0899822) q[1];
sx q[1];
rz(3.0445208) q[1];
x q[2];
rz(0.88942727) q[3];
sx q[3];
rz(-1.7173036) q[3];
sx q[3];
rz(1.5832988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.001363) q[2];
sx q[2];
rz(-1.9915853) q[2];
sx q[2];
rz(1.4353732) q[2];
rz(-2.5367149) q[3];
sx q[3];
rz(-1.6717654) q[3];
sx q[3];
rz(-2.308037) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75989938) q[0];
sx q[0];
rz(-2.9029191) q[0];
sx q[0];
rz(-2.6651486) q[0];
rz(2.6051688) q[1];
sx q[1];
rz(-1.8073578) q[1];
sx q[1];
rz(-2.776872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49052377) q[0];
sx q[0];
rz(-1.8495275) q[0];
sx q[0];
rz(1.8594478) q[0];
rz(-1.6407058) q[2];
sx q[2];
rz(-1.0349816) q[2];
sx q[2];
rz(1.792921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1583511) q[1];
sx q[1];
rz(-2.4102306) q[1];
sx q[1];
rz(-1.29391) q[1];
x q[2];
rz(2.164806) q[3];
sx q[3];
rz(-1.6315809) q[3];
sx q[3];
rz(-1.2244299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81858188) q[2];
sx q[2];
rz(-1.8135169) q[2];
sx q[2];
rz(2.2200457) q[2];
rz(-1.796683) q[3];
sx q[3];
rz(-1.709488) q[3];
sx q[3];
rz(0.014605453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0808678) q[0];
sx q[0];
rz(-0.23268555) q[0];
sx q[0];
rz(-0.099207148) q[0];
rz(-1.5469145) q[1];
sx q[1];
rz(-0.64259905) q[1];
sx q[1];
rz(-3.0329472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0799262) q[0];
sx q[0];
rz(-2.6211229) q[0];
sx q[0];
rz(-1.6322648) q[0];
rz(-0.93946879) q[2];
sx q[2];
rz(-1.924404) q[2];
sx q[2];
rz(-0.58009321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3371432) q[1];
sx q[1];
rz(-1.2986309) q[1];
sx q[1];
rz(-1.0452578) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2014162) q[3];
sx q[3];
rz(-1.7699852) q[3];
sx q[3];
rz(-1.3808113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28086153) q[2];
sx q[2];
rz(-0.18605575) q[2];
sx q[2];
rz(2.5462888) q[2];
rz(-2.2286277) q[3];
sx q[3];
rz(-2.3774417) q[3];
sx q[3];
rz(-1.3179774) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2309793) q[0];
sx q[0];
rz(-1.148104) q[0];
sx q[0];
rz(-2.8379111) q[0];
rz(0.23823711) q[1];
sx q[1];
rz(-0.99681011) q[1];
sx q[1];
rz(0.69746596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36059067) q[0];
sx q[0];
rz(-1.7502898) q[0];
sx q[0];
rz(-2.596847) q[0];
x q[1];
rz(-0.2525783) q[2];
sx q[2];
rz(-0.36732964) q[2];
sx q[2];
rz(-2.9861116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7567544) q[1];
sx q[1];
rz(-2.4765444) q[1];
sx q[1];
rz(-0.68772067) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8805299) q[3];
sx q[3];
rz(-1.1202284) q[3];
sx q[3];
rz(1.0999668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72796983) q[2];
sx q[2];
rz(-2.0406849) q[2];
sx q[2];
rz(3.0481763) q[2];
rz(-1.5591239) q[3];
sx q[3];
rz(-0.070135442) q[3];
sx q[3];
rz(-3.0541592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9974136) q[0];
sx q[0];
rz(-1.1638673) q[0];
sx q[0];
rz(-1.5735004) q[0];
rz(0.46918029) q[1];
sx q[1];
rz(-2.4041685) q[1];
sx q[1];
rz(2.266177) q[1];
rz(-0.45892528) q[2];
sx q[2];
rz(-2.060072) q[2];
sx q[2];
rz(1.2960808) q[2];
rz(2.3584414) q[3];
sx q[3];
rz(-1.7941536) q[3];
sx q[3];
rz(-0.4335946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
