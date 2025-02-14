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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(-2.5831232) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0315379) q[0];
sx q[0];
rz(-1.0194091) q[0];
sx q[0];
rz(-0.64161513) q[0];
rz(-pi) q[1];
rz(-2.9869674) q[2];
sx q[2];
rz(-2.4569656) q[2];
sx q[2];
rz(-1.3473617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3529571) q[1];
sx q[1];
rz(-0.26600263) q[1];
sx q[1];
rz(2.8255844) q[1];
x q[2];
rz(-0.73689245) q[3];
sx q[3];
rz(-2.6813563) q[3];
sx q[3];
rz(0.62084197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89995304) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-2.107035) q[2];
rz(3.0311846) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-0.30606562) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050874) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(-1.2677445) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.3892106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0076375) q[0];
sx q[0];
rz(-1.2651772) q[0];
sx q[0];
rz(-1.0108749) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6600962) q[2];
sx q[2];
rz(-1.6698368) q[2];
sx q[2];
rz(-1.5426829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25575614) q[1];
sx q[1];
rz(-0.44609082) q[1];
sx q[1];
rz(2.2913833) q[1];
rz(1.4382382) q[3];
sx q[3];
rz(-2.5131559) q[3];
sx q[3];
rz(-0.6836764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9691951) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(-1.4364852) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(3.1148541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7773975) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.2478849) q[1];
sx q[1];
rz(-1.4395813) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091451784) q[0];
sx q[0];
rz(-2.2268852) q[0];
sx q[0];
rz(1.4138282) q[0];
x q[1];
rz(1.4051938) q[2];
sx q[2];
rz(-0.51711997) q[2];
sx q[2];
rz(-1.0696326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90104687) q[1];
sx q[1];
rz(-0.26296294) q[1];
sx q[1];
rz(-2.0790624) q[1];
rz(-pi) q[2];
x q[2];
rz(2.923449) q[3];
sx q[3];
rz(-1.3036332) q[3];
sx q[3];
rz(1.4845136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8000468) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(2.3254507) q[2];
rz(-0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5661904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.4207882) q[0];
sx q[0];
rz(-1.9705462) q[0];
sx q[0];
rz(0.26914832) q[0];
rz(1.7587657) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-2.6699578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1203493) q[0];
sx q[0];
rz(-0.85130748) q[0];
sx q[0];
rz(2.9553652) q[0];
rz(-0.083375562) q[2];
sx q[2];
rz(-0.40503392) q[2];
sx q[2];
rz(-2.0138559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2216404) q[1];
sx q[1];
rz(-0.62949179) q[1];
sx q[1];
rz(-2.5821177) q[1];
rz(-0.90507953) q[3];
sx q[3];
rz(-1.2879914) q[3];
sx q[3];
rz(-2.1259411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4857594) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(1.887623) q[2];
rz(0.29821011) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(0.75697672) q[0];
rz(2.0825999) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(2.7991378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22638527) q[0];
sx q[0];
rz(-1.3801551) q[0];
sx q[0];
rz(-1.9739499) q[0];
x q[1];
rz(-1.9497237) q[2];
sx q[2];
rz(-0.73816002) q[2];
sx q[2];
rz(-1.9110695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1736502) q[1];
sx q[1];
rz(-1.0999803) q[1];
sx q[1];
rz(1.985926) q[1];
rz(1.2737688) q[3];
sx q[3];
rz(-1.9701463) q[3];
sx q[3];
rz(3.0210329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(2.5214419) q[2];
rz(0.54689637) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(-2.0641522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3243489) q[0];
sx q[0];
rz(-2.9381848) q[0];
sx q[0];
rz(-1.2212344) q[0];
rz(-2.0381894) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(2.3308636) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9098783) q[0];
sx q[0];
rz(-2.8389769) q[0];
sx q[0];
rz(2.581654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2037494) q[2];
sx q[2];
rz(-1.019291) q[2];
sx q[2];
rz(-2.1126975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.394388) q[1];
sx q[1];
rz(-1.516368) q[1];
sx q[1];
rz(1.7096976) q[1];
x q[2];
rz(2.8437988) q[3];
sx q[3];
rz(-1.9341365) q[3];
sx q[3];
rz(-2.9693574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3409884) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(-1.1654589) q[2];
rz(-3.0321339) q[3];
sx q[3];
rz(-2.1200924) q[3];
sx q[3];
rz(3.0808595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3133746) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(-0.24318801) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(-0.64116716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.015804) q[0];
sx q[0];
rz(-2.1282853) q[0];
sx q[0];
rz(-0.95161887) q[0];
rz(1.2095318) q[2];
sx q[2];
rz(-0.85339499) q[2];
sx q[2];
rz(-3.1255869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4879681) q[1];
sx q[1];
rz(-1.89978) q[1];
sx q[1];
rz(2.3423561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0585467) q[3];
sx q[3];
rz(-0.44971684) q[3];
sx q[3];
rz(-1.095677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13333653) q[2];
sx q[2];
rz(-1.1587605) q[2];
sx q[2];
rz(2.1978417) q[2];
rz(0.64822316) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7159395) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(0.42914036) q[0];
rz(-0.40402135) q[1];
sx q[1];
rz(-0.41530135) q[1];
sx q[1];
rz(-0.9187575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8480534) q[0];
sx q[0];
rz(-1.812879) q[0];
sx q[0];
rz(-2.5771067) q[0];
rz(-pi) q[1];
rz(1.8317675) q[2];
sx q[2];
rz(-1.6918285) q[2];
sx q[2];
rz(2.3132035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4705953) q[1];
sx q[1];
rz(-1.8806702) q[1];
sx q[1];
rz(-2.9122374) q[1];
rz(-pi) q[2];
rz(-1.2567344) q[3];
sx q[3];
rz(-1.5005642) q[3];
sx q[3];
rz(-2.8633871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1227485) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(2.8090779) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(1.4211753) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(0.054904003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47726705) q[0];
sx q[0];
rz(-0.83819637) q[0];
sx q[0];
rz(-0.038756702) q[0];
rz(2.730092) q[2];
sx q[2];
rz(-1.9691281) q[2];
sx q[2];
rz(2.8702132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0753638) q[1];
sx q[1];
rz(-1.7875331) q[1];
sx q[1];
rz(1.073887) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6087481) q[3];
sx q[3];
rz(-2.5379764) q[3];
sx q[3];
rz(-2.5127905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2533337) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(-2.9448275) q[2];
rz(-2.0678068) q[3];
sx q[3];
rz(-2.00627) q[3];
sx q[3];
rz(-2.4055068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9842904) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(-0.89944696) q[0];
rz(2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.4003632) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6104855) q[0];
sx q[0];
rz(-0.60043469) q[0];
sx q[0];
rz(1.9470293) q[0];
rz(-pi) q[1];
rz(2.0092404) q[2];
sx q[2];
rz(-2.6137527) q[2];
sx q[2];
rz(0.49293016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77620095) q[1];
sx q[1];
rz(-1.9579871) q[1];
sx q[1];
rz(2.4939818) q[1];
x q[2];
rz(-2.4254341) q[3];
sx q[3];
rz(-1.0993488) q[3];
sx q[3];
rz(-0.52385073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0805936) q[2];
sx q[2];
rz(-2.521535) q[2];
sx q[2];
rz(1.1617804) q[2];
rz(1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494736) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(-3.0081765) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-1.5897122) q[2];
sx q[2];
rz(-0.053413548) q[2];
sx q[2];
rz(1.5656768) q[2];
rz(2.6399163) q[3];
sx q[3];
rz(-1.6335084) q[3];
sx q[3];
rz(-0.15214534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
