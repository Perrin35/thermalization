OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90836877) q[0];
sx q[0];
rz(4.6908661) q[0];
sx q[0];
rz(9.6022845) q[0];
rz(1.5597875) q[1];
sx q[1];
rz(-1.8145365) q[1];
sx q[1];
rz(-1.1787193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8713407) q[0];
sx q[0];
rz(-0.9862904) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(-pi) q[1];
rz(0.92843797) q[2];
sx q[2];
rz(-1.8132767) q[2];
sx q[2];
rz(1.2280304) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5753769) q[1];
sx q[1];
rz(-2.1388091) q[1];
sx q[1];
rz(0.52277188) q[1];
rz(1.7186382) q[3];
sx q[3];
rz(-2.4018722) q[3];
sx q[3];
rz(1.5802368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.061964758) q[2];
sx q[2];
rz(-0.92755932) q[2];
sx q[2];
rz(-0.24688841) q[2];
rz(0.43761349) q[3];
sx q[3];
rz(-2.103431) q[3];
sx q[3];
rz(2.8810697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0831809) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(-2.0368982) q[0];
rz(0.80766922) q[1];
sx q[1];
rz(-2.3876815) q[1];
sx q[1];
rz(0.63899904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09428962) q[0];
sx q[0];
rz(-1.2942137) q[0];
sx q[0];
rz(-1.7533789) q[0];
rz(0.17484395) q[2];
sx q[2];
rz(-0.27894601) q[2];
sx q[2];
rz(-2.1871559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5377755) q[1];
sx q[1];
rz(-0.81934419) q[1];
sx q[1];
rz(0.50756331) q[1];
x q[2];
rz(-1.9023537) q[3];
sx q[3];
rz(-1.8428118) q[3];
sx q[3];
rz(2.5355458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69940007) q[2];
sx q[2];
rz(-1.3905808) q[2];
sx q[2];
rz(-0.81713444) q[2];
rz(2.6160431) q[3];
sx q[3];
rz(-2.3216129) q[3];
sx q[3];
rz(1.9513963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007871) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(0.15980414) q[0];
rz(2.6218759) q[1];
sx q[1];
rz(-1.0451885) q[1];
sx q[1];
rz(2.2249178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5938497) q[0];
sx q[0];
rz(-1.5446988) q[0];
sx q[0];
rz(2.5711077) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26600142) q[2];
sx q[2];
rz(-0.75845462) q[2];
sx q[2];
rz(0.11125362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8868272) q[1];
sx q[1];
rz(-1.6187833) q[1];
sx q[1];
rz(1.3814203) q[1];
x q[2];
rz(2.5289747) q[3];
sx q[3];
rz(-1.6809941) q[3];
sx q[3];
rz(-1.2904905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89797574) q[2];
sx q[2];
rz(-0.50668442) q[2];
sx q[2];
rz(-1.766073) q[2];
rz(0.033179387) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(-0.83056617) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0631977) q[0];
sx q[0];
rz(-2.4802408) q[0];
sx q[0];
rz(0.98264328) q[0];
rz(0.64535087) q[1];
sx q[1];
rz(-1.9159562) q[1];
sx q[1];
rz(0.4172079) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010734) q[0];
sx q[0];
rz(-2.5357781) q[0];
sx q[0];
rz(0.54744253) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.190416) q[2];
sx q[2];
rz(-1.5065168) q[2];
sx q[2];
rz(2.1716139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55298333) q[1];
sx q[1];
rz(-0.89834301) q[1];
sx q[1];
rz(-1.6469127) q[1];
rz(2.6989815) q[3];
sx q[3];
rz(-2.5472982) q[3];
sx q[3];
rz(-1.3292918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.30503201) q[2];
sx q[2];
rz(-2.8053668) q[2];
sx q[2];
rz(-1.9545371) q[2];
rz(0.120397) q[3];
sx q[3];
rz(-1.7359366) q[3];
sx q[3];
rz(-2.9482237) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3099986) q[0];
sx q[0];
rz(-0.11045063) q[0];
sx q[0];
rz(-0.67685342) q[0];
rz(1.5325158) q[1];
sx q[1];
rz(-2.0307505) q[1];
sx q[1];
rz(-2.9422876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520357) q[0];
sx q[0];
rz(-1.3180719) q[0];
sx q[0];
rz(0.093555497) q[0];
rz(-pi) q[1];
rz(-0.3466269) q[2];
sx q[2];
rz(-1.4557654) q[2];
sx q[2];
rz(2.4644867) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3518954) q[1];
sx q[1];
rz(-1.9421862) q[1];
sx q[1];
rz(-1.8892193) q[1];
x q[2];
rz(1.8477047) q[3];
sx q[3];
rz(-1.5770638) q[3];
sx q[3];
rz(-2.4771877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6469595) q[2];
sx q[2];
rz(-2.5854752) q[2];
sx q[2];
rz(-0.49590084) q[2];
rz(0.9209218) q[3];
sx q[3];
rz(-1.0420739) q[3];
sx q[3];
rz(0.31057772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344236) q[0];
sx q[0];
rz(-1.2260219) q[0];
sx q[0];
rz(-2.7219211) q[0];
rz(0.93027973) q[1];
sx q[1];
rz(-1.0134965) q[1];
sx q[1];
rz(2.1297661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8735805) q[0];
sx q[0];
rz(-2.4873864) q[0];
sx q[0];
rz(-0.081120567) q[0];
rz(0.65134766) q[2];
sx q[2];
rz(-0.22964165) q[2];
sx q[2];
rz(2.1852213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7020445) q[1];
sx q[1];
rz(-1.3633894) q[1];
sx q[1];
rz(-2.1316567) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18475549) q[3];
sx q[3];
rz(-2.1596381) q[3];
sx q[3];
rz(2.7240804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4100627) q[2];
sx q[2];
rz(-0.56649929) q[2];
sx q[2];
rz(-1.1140964) q[2];
rz(1.3725494) q[3];
sx q[3];
rz(-2.1362344) q[3];
sx q[3];
rz(-2.4323288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1019185) q[0];
sx q[0];
rz(-1.2693951) q[0];
sx q[0];
rz(-1.6495548) q[0];
rz(2.6986625) q[1];
sx q[1];
rz(-1.2716753) q[1];
sx q[1];
rz(0.84239229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32552606) q[0];
sx q[0];
rz(-0.93579817) q[0];
sx q[0];
rz(0.86845894) q[0];
rz(-2.2350253) q[2];
sx q[2];
rz(-1.3593055) q[2];
sx q[2];
rz(1.4480526) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8460626) q[1];
sx q[1];
rz(-1.1605442) q[1];
sx q[1];
rz(2.2890291) q[1];
rz(1.674323) q[3];
sx q[3];
rz(-1.6134422) q[3];
sx q[3];
rz(1.9089684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3288021) q[2];
sx q[2];
rz(-2.2219358) q[2];
sx q[2];
rz(-0.021050464) q[2];
rz(0.99714315) q[3];
sx q[3];
rz(-2.5992584) q[3];
sx q[3];
rz(1.5890582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7954623) q[0];
sx q[0];
rz(-1.2393476) q[0];
sx q[0];
rz(1.1783266) q[0];
rz(-1.0109674) q[1];
sx q[1];
rz(-1.6594454) q[1];
sx q[1];
rz(2.7706026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3431992) q[0];
sx q[0];
rz(-1.4552552) q[0];
sx q[0];
rz(-0.2327524) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.095711555) q[2];
sx q[2];
rz(-1.7388027) q[2];
sx q[2];
rz(-1.8883349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64196649) q[1];
sx q[1];
rz(-1.2807944) q[1];
sx q[1];
rz(2.9846342) q[1];
x q[2];
rz(0.99222547) q[3];
sx q[3];
rz(-0.5360837) q[3];
sx q[3];
rz(1.2257731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79169881) q[2];
sx q[2];
rz(-2.0588304) q[2];
sx q[2];
rz(0.91378158) q[2];
rz(-0.53798109) q[3];
sx q[3];
rz(-0.83298433) q[3];
sx q[3];
rz(-2.7624281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051788483) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(2.9096933) q[0];
rz(-1.6357577) q[1];
sx q[1];
rz(-1.6927405) q[1];
sx q[1];
rz(0.64461446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1132677) q[0];
sx q[0];
rz(-1.9132943) q[0];
sx q[0];
rz(-1.2069846) q[0];
x q[1];
rz(0.03657954) q[2];
sx q[2];
rz(-2.5365005) q[2];
sx q[2];
rz(-0.99672752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7773103) q[1];
sx q[1];
rz(-0.16723086) q[1];
sx q[1];
rz(1.689453) q[1];
rz(-2.4394611) q[3];
sx q[3];
rz(-1.9117179) q[3];
sx q[3];
rz(-3.0263755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5171648) q[2];
sx q[2];
rz(-1.6021873) q[2];
sx q[2];
rz(-1.2702764) q[2];
rz(-2.7020448) q[3];
sx q[3];
rz(-2.5054273) q[3];
sx q[3];
rz(2.6308681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9257833) q[0];
sx q[0];
rz(-1.2405115) q[0];
sx q[0];
rz(-1.0260169) q[0];
rz(-2.6846474) q[1];
sx q[1];
rz(-0.89070717) q[1];
sx q[1];
rz(-0.93422008) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608604) q[0];
sx q[0];
rz(-0.44507682) q[0];
sx q[0];
rz(0.31240065) q[0];
x q[1];
rz(-1.9719741) q[2];
sx q[2];
rz(-0.27571644) q[2];
sx q[2];
rz(1.8091615) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6187541) q[1];
sx q[1];
rz(-1.4013023) q[1];
sx q[1];
rz(2.3970669) q[1];
rz(-pi) q[2];
rz(-2.281743) q[3];
sx q[3];
rz(-2.352664) q[3];
sx q[3];
rz(-3.0956059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43941867) q[2];
sx q[2];
rz(-0.69434035) q[2];
sx q[2];
rz(-1.0609421) q[2];
rz(-2.151978) q[3];
sx q[3];
rz(-1.3720392) q[3];
sx q[3];
rz(0.26452479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2117352) q[0];
sx q[0];
rz(-1.095357) q[0];
sx q[0];
rz(2.3768421) q[0];
rz(1.9542971) q[1];
sx q[1];
rz(-1.4359513) q[1];
sx q[1];
rz(-0.5618701) q[1];
rz(-1.622772) q[2];
sx q[2];
rz(-1.5511654) q[2];
sx q[2];
rz(1.8775107) q[2];
rz(0.76239794) q[3];
sx q[3];
rz(-0.60653209) q[3];
sx q[3];
rz(-1.7222986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
