OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1342993) q[0];
sx q[0];
rz(-1.1987885) q[0];
sx q[0];
rz(-0.78517908) q[0];
rz(2.5685318) q[1];
sx q[1];
rz(2.1905724) q[1];
sx q[1];
rz(11.934927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19772675) q[0];
sx q[0];
rz(-0.80538865) q[0];
sx q[0];
rz(0.95279537) q[0];
rz(-pi) q[1];
rz(2.3518121) q[2];
sx q[2];
rz(-2.1319816) q[2];
sx q[2];
rz(-1.1488069) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3062485) q[1];
sx q[1];
rz(-2.1901096) q[1];
sx q[1];
rz(1.033551) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20578014) q[3];
sx q[3];
rz(-0.45318402) q[3];
sx q[3];
rz(-0.024621016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.637392) q[2];
sx q[2];
rz(-1.4376419) q[2];
sx q[2];
rz(0.27153095) q[2];
rz(-3.0951989) q[3];
sx q[3];
rz(-1.7312739) q[3];
sx q[3];
rz(-0.36762777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72535998) q[0];
sx q[0];
rz(-0.027149057) q[0];
sx q[0];
rz(-0.44560462) q[0];
rz(0.36000571) q[1];
sx q[1];
rz(-1.5841192) q[1];
sx q[1];
rz(2.1645567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038526857) q[0];
sx q[0];
rz(-1.0630106) q[0];
sx q[0];
rz(0.610791) q[0];
rz(-pi) q[1];
rz(-2.2266892) q[2];
sx q[2];
rz(-0.96157461) q[2];
sx q[2];
rz(0.34944361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9865446) q[1];
sx q[1];
rz(-2.1594271) q[1];
sx q[1];
rz(0.18208336) q[1];
x q[2];
rz(0.70313022) q[3];
sx q[3];
rz(-1.0323914) q[3];
sx q[3];
rz(-2.3634274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6147324) q[2];
sx q[2];
rz(-1.7975668) q[2];
sx q[2];
rz(2.1993707) q[2];
rz(2.6451536) q[3];
sx q[3];
rz(-1.6783291) q[3];
sx q[3];
rz(-2.0739323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577069) q[0];
sx q[0];
rz(-1.3297465) q[0];
sx q[0];
rz(0.65621334) q[0];
rz(-2.6612142) q[1];
sx q[1];
rz(-0.96157938) q[1];
sx q[1];
rz(2.5286455) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543535) q[0];
sx q[0];
rz(-1.4738074) q[0];
sx q[0];
rz(1.9796275) q[0];
x q[1];
rz(-1.5781338) q[2];
sx q[2];
rz(-2.3068608) q[2];
sx q[2];
rz(-1.8328345) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9631098) q[1];
sx q[1];
rz(-1.6043538) q[1];
sx q[1];
rz(1.5419649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49895371) q[3];
sx q[3];
rz(-1.6315579) q[3];
sx q[3];
rz(-2.2419597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2059242) q[2];
sx q[2];
rz(-2.084145) q[2];
sx q[2];
rz(0.52440468) q[2];
rz(-2.2762401) q[3];
sx q[3];
rz(-2.083358) q[3];
sx q[3];
rz(-0.66044468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.653729) q[0];
sx q[0];
rz(-1.8346584) q[0];
sx q[0];
rz(-3.0453239) q[0];
rz(3.1277711) q[1];
sx q[1];
rz(-2.6410069) q[1];
sx q[1];
rz(1.0523419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88599624) q[0];
sx q[0];
rz(-1.4491206) q[0];
sx q[0];
rz(-0.047329655) q[0];
rz(2.0180447) q[2];
sx q[2];
rz(-2.5442225) q[2];
sx q[2];
rz(1.1505466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4197156) q[1];
sx q[1];
rz(-1.8977563) q[1];
sx q[1];
rz(1.378233) q[1];
x q[2];
rz(0.035798612) q[3];
sx q[3];
rz(-2.0960663) q[3];
sx q[3];
rz(0.73204277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72157613) q[2];
sx q[2];
rz(-2.6105011) q[2];
sx q[2];
rz(2.6225923) q[2];
rz(1.0985993) q[3];
sx q[3];
rz(-1.4562621) q[3];
sx q[3];
rz(2.4115653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62243432) q[0];
sx q[0];
rz(-2.9485478) q[0];
sx q[0];
rz(0.73424196) q[0];
rz(-1.2892067) q[1];
sx q[1];
rz(-1.0451008) q[1];
sx q[1];
rz(2.2330914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168981) q[0];
sx q[0];
rz(-2.1229345) q[0];
sx q[0];
rz(-0.79663251) q[0];
x q[1];
rz(0.20573767) q[2];
sx q[2];
rz(-2.5605953) q[2];
sx q[2];
rz(-0.28827039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3215103) q[1];
sx q[1];
rz(-1.501298) q[1];
sx q[1];
rz(0.837079) q[1];
rz(-pi) q[2];
rz(-1.8374422) q[3];
sx q[3];
rz(-1.5943267) q[3];
sx q[3];
rz(-1.2465828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.044540731) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(-0.68361863) q[2];
rz(-1.094007) q[3];
sx q[3];
rz(-1.8180327) q[3];
sx q[3];
rz(-1.5723642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053726824) q[0];
sx q[0];
rz(-1.6931067) q[0];
sx q[0];
rz(0.42561612) q[0];
rz(0.60676891) q[1];
sx q[1];
rz(-0.68777045) q[1];
sx q[1];
rz(1.2026131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2736061) q[0];
sx q[0];
rz(-2.3487954) q[0];
sx q[0];
rz(0.85203627) q[0];
x q[1];
rz(-1.556237) q[2];
sx q[2];
rz(-1.5049767) q[2];
sx q[2];
rz(-1.3902537) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7640671) q[1];
sx q[1];
rz(-2.0646618) q[1];
sx q[1];
rz(2.4000044) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8011064) q[3];
sx q[3];
rz(-1.7969683) q[3];
sx q[3];
rz(-0.61128547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6732424) q[2];
sx q[2];
rz(-2.612256) q[2];
sx q[2];
rz(-0.62874111) q[2];
rz(2.8041503) q[3];
sx q[3];
rz(-1.3597666) q[3];
sx q[3];
rz(0.76702816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.5479946) q[0];
sx q[0];
rz(-3.0715521) q[0];
sx q[0];
rz(-1.8016169) q[0];
rz(-3.0361259) q[1];
sx q[1];
rz(-0.98140812) q[1];
sx q[1];
rz(-2.9973105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893989) q[0];
sx q[0];
rz(-1.1849863) q[0];
sx q[0];
rz(-1.2560638) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3445129) q[2];
sx q[2];
rz(-0.72090518) q[2];
sx q[2];
rz(3.0932567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49715091) q[1];
sx q[1];
rz(-1.9195917) q[1];
sx q[1];
rz(-1.7921053) q[1];
x q[2];
rz(-1.027368) q[3];
sx q[3];
rz(-2.1680834) q[3];
sx q[3];
rz(1.3940108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1062801) q[2];
sx q[2];
rz(-1.9929746) q[2];
sx q[2];
rz(2.3036352) q[2];
rz(-0.49977866) q[3];
sx q[3];
rz(-0.79334799) q[3];
sx q[3];
rz(2.4404073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4226294) q[0];
sx q[0];
rz(-2.0327649) q[0];
sx q[0];
rz(-0.75123373) q[0];
rz(-1.2535837) q[1];
sx q[1];
rz(-1.3026078) q[1];
sx q[1];
rz(1.7734897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7984895) q[0];
sx q[0];
rz(-1.1310078) q[0];
sx q[0];
rz(-0.90476157) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31603916) q[2];
sx q[2];
rz(-1.2908986) q[2];
sx q[2];
rz(0.47376493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66361278) q[1];
sx q[1];
rz(-2.3638032) q[1];
sx q[1];
rz(1.6915093) q[1];
rz(-pi) q[2];
rz(2.7799928) q[3];
sx q[3];
rz(-2.1146833) q[3];
sx q[3];
rz(1.4023449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3086243) q[2];
sx q[2];
rz(-2.7699351) q[2];
sx q[2];
rz(-0.0011681636) q[2];
rz(2.809281) q[3];
sx q[3];
rz(-1.4609591) q[3];
sx q[3];
rz(0.47962475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7697656) q[0];
sx q[0];
rz(-1.9244939) q[0];
sx q[0];
rz(1.0960854) q[0];
rz(-1.3990043) q[1];
sx q[1];
rz(-1.636248) q[1];
sx q[1];
rz(-0.91986626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77114785) q[0];
sx q[0];
rz(-0.4397529) q[0];
sx q[0];
rz(-1.670776) q[0];
rz(2.8613352) q[2];
sx q[2];
rz(-1.9748903) q[2];
sx q[2];
rz(1.0755444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2640317) q[1];
sx q[1];
rz(-1.4204867) q[1];
sx q[1];
rz(-0.64964215) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9971531) q[3];
sx q[3];
rz(-1.0424992) q[3];
sx q[3];
rz(2.4424255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2190735) q[2];
sx q[2];
rz(-2.8963431) q[2];
sx q[2];
rz(-1.6314487) q[2];
rz(0.63589969) q[3];
sx q[3];
rz(-2.4363775) q[3];
sx q[3];
rz(-0.51457921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.9566327) q[0];
sx q[0];
rz(-0.867046) q[0];
sx q[0];
rz(1.4956723) q[0];
rz(-0.11790672) q[1];
sx q[1];
rz(-1.7472183) q[1];
sx q[1];
rz(0.25679055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1908156) q[0];
sx q[0];
rz(-2.7047892) q[0];
sx q[0];
rz(1.703116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7071827) q[2];
sx q[2];
rz(-2.5108178) q[2];
sx q[2];
rz(1.3303016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49918711) q[1];
sx q[1];
rz(-0.39246628) q[1];
sx q[1];
rz(-0.74340238) q[1];
rz(-0.54739807) q[3];
sx q[3];
rz(-1.8845356) q[3];
sx q[3];
rz(2.1496839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7439277) q[2];
sx q[2];
rz(-0.74349082) q[2];
sx q[2];
rz(2.1916154) q[2];
rz(2.5395565) q[3];
sx q[3];
rz(-1.9046015) q[3];
sx q[3];
rz(0.33815798) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4468741) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(-2.1318204) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(-0.80712423) q[2];
sx q[2];
rz(-2.7257498) q[2];
sx q[2];
rz(2.7422703) q[2];
rz(-0.39691858) q[3];
sx q[3];
rz(-2.3229058) q[3];
sx q[3];
rz(-2.3380052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
