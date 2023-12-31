OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(-0.64557689) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49255532) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(1.9306246) q[0];
x q[1];
rz(2.2013118) q[2];
sx q[2];
rz(-1.5144381) q[2];
sx q[2];
rz(1.2727357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9420535) q[1];
sx q[1];
rz(-0.76738165) q[1];
sx q[1];
rz(-0.3120504) q[1];
rz(-pi) q[2];
rz(-0.50819355) q[3];
sx q[3];
rz(-2.4262706) q[3];
sx q[3];
rz(3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-2.5966068) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(2.3166336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4175407) q[0];
sx q[0];
rz(-1.9736971) q[0];
sx q[0];
rz(-0.20362644) q[0];
rz(-1.5603793) q[2];
sx q[2];
rz(-1.1679107) q[2];
sx q[2];
rz(-1.252623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7823239) q[1];
sx q[1];
rz(-1.4546848) q[1];
sx q[1];
rz(-0.5639204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9460658) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(2.5415681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(3.1030531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8124213) q[0];
sx q[0];
rz(-0.0077795452) q[0];
sx q[0];
rz(-1.6532941) q[0];
rz(-pi) q[1];
rz(0.58263393) q[2];
sx q[2];
rz(-0.27432549) q[2];
sx q[2];
rz(-0.81253101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4290532) q[1];
sx q[1];
rz(-0.4936115) q[1];
sx q[1];
rz(-1.4406956) q[1];
rz(-pi) q[2];
rz(2.3452957) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(0.89087668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12476607) q[0];
sx q[0];
rz(-0.15206465) q[0];
sx q[0];
rz(2.2269339) q[0];
rz(-1.0267369) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.7550215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86707592) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(-2.0348674) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5024662) q[3];
sx q[3];
rz(-1.5838606) q[3];
sx q[3];
rz(-2.5845205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(-1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.224684) q[0];
sx q[0];
rz(-1.2901257) q[0];
sx q[0];
rz(1.0018437) q[0];
x q[1];
rz(2.4079608) q[2];
sx q[2];
rz(-2.1918104) q[2];
sx q[2];
rz(2.6526407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34199076) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(-2.0288543) q[1];
rz(-1.3729172) q[3];
sx q[3];
rz(-2.0298925) q[3];
sx q[3];
rz(2.0898233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75509214) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(2.7639672) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5650425) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(0.25175005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50778163) q[2];
sx q[2];
rz(-2.0785329) q[2];
sx q[2];
rz(0.43916647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5988016) q[1];
sx q[1];
rz(-0.32074499) q[1];
sx q[1];
rz(-2.9221605) q[1];
rz(-pi) q[2];
rz(1.6842501) q[3];
sx q[3];
rz(-2.2621584) q[3];
sx q[3];
rz(1.3464348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6043828) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.4253915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0910949) q[0];
sx q[0];
rz(-1.1384581) q[0];
sx q[0];
rz(1.7763441) q[0];
rz(2.695589) q[2];
sx q[2];
rz(-1.6288695) q[2];
sx q[2];
rz(0.29689483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2196301) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(2.4835543) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1281934) q[3];
sx q[3];
rz(-1.0671167) q[3];
sx q[3];
rz(-0.79813938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.210775) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(0.67529768) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589688) q[0];
sx q[0];
rz(-1.0685182) q[0];
sx q[0];
rz(1.9320022) q[0];
rz(-2.2194527) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(-2.5285072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68554316) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(-0.088940253) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3648871) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.3148274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7816759) q[0];
sx q[0];
rz(-1.3206375) q[0];
sx q[0];
rz(-0.71235384) q[0];
rz(-pi) q[1];
rz(-2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(3.0968015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7278553) q[1];
sx q[1];
rz(-1.9836042) q[1];
sx q[1];
rz(-2.2933943) q[1];
rz(-pi) q[2];
rz(-2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(-1.0178125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(2.749696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6673198) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(-1.5865109) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0060843) q[2];
sx q[2];
rz(-0.91943179) q[2];
sx q[2];
rz(-1.1609921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7614903) q[1];
sx q[1];
rz(-0.83161608) q[1];
sx q[1];
rz(-2.889061) q[1];
rz(-pi) q[2];
rz(-2.6565353) q[3];
sx q[3];
rz(-2.5254446) q[3];
sx q[3];
rz(-2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(2.3618868) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(-0.34655456) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
