OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(-2.1704817) q[0];
sx q[0];
rz(0.71075332) q[0];
rz(-2.1739668) q[1];
sx q[1];
rz(-0.014048014) q[1];
sx q[1];
rz(2.3486121) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7638057) q[0];
sx q[0];
rz(-1.9024693) q[0];
sx q[0];
rz(0.8140696) q[0];
rz(-pi) q[1];
rz(-1.8618199) q[2];
sx q[2];
rz(-3.0918192) q[2];
sx q[2];
rz(-0.39779824) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88411623) q[1];
sx q[1];
rz(-2.3168644) q[1];
sx q[1];
rz(2.5349239) q[1];
x q[2];
rz(-2.1737384) q[3];
sx q[3];
rz(-1.0558898) q[3];
sx q[3];
rz(-0.076857278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27796081) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(-2.5521736) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(2.2404501) q[0];
rz(0.98183739) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(-2.1493886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311148) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(0.94816072) q[0];
rz(-pi) q[1];
rz(2.9770933) q[2];
sx q[2];
rz(-0.82969159) q[2];
sx q[2];
rz(2.6890597) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34582576) q[1];
sx q[1];
rz(-0.92928934) q[1];
sx q[1];
rz(-0.87916763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68282376) q[3];
sx q[3];
rz(-1.0299333) q[3];
sx q[3];
rz(-1.4861432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.09484) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(2.181634) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(2.4175194) q[0];
rz(0.11101668) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(-0.012880005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2483978) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(3.061122) q[0];
rz(0.78978844) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(-2.546026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9835123) q[1];
sx q[1];
rz(-1.8077824) q[1];
sx q[1];
rz(-1.7416341) q[1];
rz(-pi) q[2];
rz(-2.0184085) q[3];
sx q[3];
rz(-0.73799282) q[3];
sx q[3];
rz(0.18692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(0.4970099) q[3];
sx q[3];
rz(-0.23247601) q[3];
sx q[3];
rz(2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67320353) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(-0.47624269) q[0];
rz(-2.6561123) q[1];
sx q[1];
rz(-2.6464033) q[1];
sx q[1];
rz(-2.9229497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784558) q[0];
sx q[0];
rz(-1.7972203) q[0];
sx q[0];
rz(-1.9853659) q[0];
rz(-1.5831469) q[2];
sx q[2];
rz(-2.2943779) q[2];
sx q[2];
rz(-0.45958322) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87590295) q[1];
sx q[1];
rz(-2.4445718) q[1];
sx q[1];
rz(-2.4494996) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3371631) q[3];
sx q[3];
rz(-1.2645742) q[3];
sx q[3];
rz(0.15627565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.789088) q[2];
sx q[2];
rz(-2.3829491) q[2];
sx q[2];
rz(-0.57141203) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-2.5550227) q[3];
sx q[3];
rz(-0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649696) q[0];
sx q[0];
rz(-0.40410703) q[0];
sx q[0];
rz(-0.47082666) q[0];
rz(0.91104031) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-0.43112531) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3456988) q[0];
sx q[0];
rz(-0.78471334) q[0];
sx q[0];
rz(0.037952947) q[0];
x q[1];
rz(-0.25283739) q[2];
sx q[2];
rz(-0.95273861) q[2];
sx q[2];
rz(-1.1065266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70111245) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(-3.1140226) q[1];
rz(-pi) q[2];
rz(-2.9399559) q[3];
sx q[3];
rz(-1.7072225) q[3];
sx q[3];
rz(-1.0659983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-0.63569331) q[3];
sx q[3];
rz(-2.671833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895198) q[0];
sx q[0];
rz(-0.85318035) q[0];
sx q[0];
rz(-3.0186655) q[0];
rz(2.7561103) q[1];
sx q[1];
rz(-0.79817927) q[1];
sx q[1];
rz(0.72365671) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6521344) q[0];
sx q[0];
rz(-1.3969027) q[0];
sx q[0];
rz(0.045920865) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0711818) q[2];
sx q[2];
rz(-2.7636409) q[2];
sx q[2];
rz(2.9674781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78711975) q[1];
sx q[1];
rz(-2.2406089) q[1];
sx q[1];
rz(2.3666827) q[1];
rz(-2.0615929) q[3];
sx q[3];
rz(-2.7088968) q[3];
sx q[3];
rz(-0.015263488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3896997) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(-2.4683118) q[2];
rz(-0.26345396) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.70169705) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(-2.6252966) q[0];
rz(0.75597489) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(-0.36852536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4936016) q[0];
sx q[0];
rz(-0.89009919) q[0];
sx q[0];
rz(2.2590911) q[0];
x q[1];
rz(0.045727878) q[2];
sx q[2];
rz(-0.70529443) q[2];
sx q[2];
rz(-2.4207585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8175537) q[1];
sx q[1];
rz(-1.4020668) q[1];
sx q[1];
rz(0.1076474) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6266699) q[3];
sx q[3];
rz(-2.2201) q[3];
sx q[3];
rz(-2.2442371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64510173) q[2];
sx q[2];
rz(-0.056702159) q[2];
sx q[2];
rz(0.48218316) q[2];
rz(0.21969806) q[3];
sx q[3];
rz(-0.69742656) q[3];
sx q[3];
rz(2.4612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41736233) q[0];
sx q[0];
rz(-2.9927232) q[0];
sx q[0];
rz(2.505488) q[0];
rz(0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(-0.87463921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2518305) q[0];
sx q[0];
rz(-1.5127851) q[0];
sx q[0];
rz(-1.7202352) q[0];
rz(-2.4061255) q[2];
sx q[2];
rz(-0.16263419) q[2];
sx q[2];
rz(-2.5474608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88747178) q[1];
sx q[1];
rz(-1.9560818) q[1];
sx q[1];
rz(-2.0068568) q[1];
rz(0.80031826) q[3];
sx q[3];
rz(-0.85552514) q[3];
sx q[3];
rz(-3.0190938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26802289) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(-0.91656172) q[2];
rz(0.77740866) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(0.53434813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022631835) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-0.62260038) q[1];
sx q[1];
rz(2.5737305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647301) q[0];
sx q[0];
rz(-1.9807388) q[0];
sx q[0];
rz(-2.7565844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90654208) q[2];
sx q[2];
rz(-0.34053206) q[2];
sx q[2];
rz(0.9685002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1901924) q[1];
sx q[1];
rz(-1.4245598) q[1];
sx q[1];
rz(-1.7758572) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80983272) q[3];
sx q[3];
rz(-0.91512915) q[3];
sx q[3];
rz(-1.1717494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(2.589321) q[2];
rz(0.28198379) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(-0.10274398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2118537) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(-0.1317568) q[0];
rz(-2.6051104) q[1];
sx q[1];
rz(-2.9832612) q[1];
sx q[1];
rz(0.90824711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73777481) q[0];
sx q[0];
rz(-1.5409011) q[0];
sx q[0];
rz(2.2653518) q[0];
x q[1];
rz(1.1029673) q[2];
sx q[2];
rz(-1.8059633) q[2];
sx q[2];
rz(0.18593341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40501198) q[1];
sx q[1];
rz(-2.2078276) q[1];
sx q[1];
rz(2.4819961) q[1];
rz(-pi) q[2];
rz(1.796714) q[3];
sx q[3];
rz(-1.572567) q[3];
sx q[3];
rz(-1.2745119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6324255) q[2];
sx q[2];
rz(-2.2278892) q[2];
sx q[2];
rz(2.8636279) q[2];
rz(-2.5748504) q[3];
sx q[3];
rz(-0.20213474) q[3];
sx q[3];
rz(2.22866) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3315898) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(-2.6592061) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-2.161088) q[2];
sx q[2];
rz(-1.7750778) q[2];
sx q[2];
rz(1.9951174) q[2];
rz(3.1101221) q[3];
sx q[3];
rz(-1.0647872) q[3];
sx q[3];
rz(-2.9797462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
