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
rz(1.6990868) q[0];
sx q[0];
rz(-1.8449755) q[0];
sx q[0];
rz(-1.698864) q[0];
rz(2.7891085) q[1];
sx q[1];
rz(3.6678996) q[1];
sx q[1];
rz(7.6511135) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57940021) q[0];
sx q[0];
rz(-0.44792563) q[0];
sx q[0];
rz(0.62904398) q[0];
x q[1];
rz(1.6085195) q[2];
sx q[2];
rz(-1.2153271) q[2];
sx q[2];
rz(-1.0743574) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1795111) q[1];
sx q[1];
rz(-1.700811) q[1];
sx q[1];
rz(2.6165753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.164293) q[3];
sx q[3];
rz(-0.90850368) q[3];
sx q[3];
rz(0.73847258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0679396) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(1.937872) q[2];
rz(-2.1667571) q[3];
sx q[3];
rz(-2.1741512) q[3];
sx q[3];
rz(1.8956641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342733) q[0];
sx q[0];
rz(-1.0618671) q[0];
sx q[0];
rz(-0.88428307) q[0];
rz(-0.70941225) q[1];
sx q[1];
rz(-1.246289) q[1];
sx q[1];
rz(-0.90219227) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220778) q[0];
sx q[0];
rz(-1.7762707) q[0];
sx q[0];
rz(0.46066649) q[0];
x q[1];
rz(-0.96376597) q[2];
sx q[2];
rz(-1.4861408) q[2];
sx q[2];
rz(-0.36106685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.480011) q[1];
sx q[1];
rz(-1.8152092) q[1];
sx q[1];
rz(-0.30765012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7092) q[3];
sx q[3];
rz(-1.6507727) q[3];
sx q[3];
rz(-0.13025912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0238637) q[2];
sx q[2];
rz(-0.29568299) q[2];
sx q[2];
rz(1.6531061) q[2];
rz(1.214341) q[3];
sx q[3];
rz(-1.4597273) q[3];
sx q[3];
rz(-1.8620209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073256) q[0];
sx q[0];
rz(-1.3995582) q[0];
sx q[0];
rz(-0.16635995) q[0];
rz(0.69794377) q[1];
sx q[1];
rz(-0.78015399) q[1];
sx q[1];
rz(1.4770329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.03334) q[0];
sx q[0];
rz(-1.1287924) q[0];
sx q[0];
rz(-3.1382794) q[0];
rz(-2.338411) q[2];
sx q[2];
rz(-0.381261) q[2];
sx q[2];
rz(-1.8800637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0946413) q[1];
sx q[1];
rz(-1.1186386) q[1];
sx q[1];
rz(0.96479123) q[1];
rz(-pi) q[2];
rz(3.1055177) q[3];
sx q[3];
rz(-2.6429308) q[3];
sx q[3];
rz(2.7447678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.6537689) q[2];
sx q[2];
rz(-1.9705801) q[2];
sx q[2];
rz(0.6146532) q[2];
rz(-1.4608308) q[3];
sx q[3];
rz(-1.5771644) q[3];
sx q[3];
rz(3.0998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.3294687) q[0];
sx q[0];
rz(1.6229269) q[0];
rz(-0.80865771) q[1];
sx q[1];
rz(-1.8452019) q[1];
sx q[1];
rz(-0.048695806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2364427) q[0];
sx q[0];
rz(-0.76950162) q[0];
sx q[0];
rz(2.6506846) q[0];
rz(0.38191958) q[2];
sx q[2];
rz(-1.4617698) q[2];
sx q[2];
rz(-2.7545021) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9756226) q[1];
sx q[1];
rz(-0.73117729) q[1];
sx q[1];
rz(2.380409) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1118495) q[3];
sx q[3];
rz(-2.4623639) q[3];
sx q[3];
rz(-1.3857057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.101717) q[2];
sx q[2];
rz(-2.3287435) q[2];
sx q[2];
rz(1.8854878) q[2];
rz(0.060955437) q[3];
sx q[3];
rz(-0.90164369) q[3];
sx q[3];
rz(2.2140293) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98546472) q[0];
sx q[0];
rz(-2.0138854) q[0];
sx q[0];
rz(0.10966478) q[0];
rz(-2.1127286) q[1];
sx q[1];
rz(-1.2867462) q[1];
sx q[1];
rz(-1.5922155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76413918) q[0];
sx q[0];
rz(-1.8528588) q[0];
sx q[0];
rz(-0.46746032) q[0];
rz(-pi) q[1];
rz(-2.0711871) q[2];
sx q[2];
rz(-1.8922046) q[2];
sx q[2];
rz(-1.6381581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.938443) q[1];
sx q[1];
rz(-2.0324273) q[1];
sx q[1];
rz(1.0814971) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6084592) q[3];
sx q[3];
rz(-1.6051014) q[3];
sx q[3];
rz(-1.487285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5900383) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(0.48198286) q[2];
rz(-0.41989741) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(2.4317252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14966203) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(-0.71257198) q[0];
rz(2.271671) q[1];
sx q[1];
rz(-2.7674119) q[1];
sx q[1];
rz(-3.1178927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256692) q[0];
sx q[0];
rz(-2.5953889) q[0];
sx q[0];
rz(-2.8234981) q[0];
x q[1];
rz(-2.4392088) q[2];
sx q[2];
rz(-2.0310035) q[2];
sx q[2];
rz(2.5277939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5795676) q[1];
sx q[1];
rz(-1.7365841) q[1];
sx q[1];
rz(3.1400984) q[1];
rz(-pi) q[2];
rz(-2.7530144) q[3];
sx q[3];
rz(-0.51296189) q[3];
sx q[3];
rz(-1.4298317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80453834) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(1.1402593) q[2];
rz(-1.987847) q[3];
sx q[3];
rz(-1.2146344) q[3];
sx q[3];
rz(-1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425728) q[0];
sx q[0];
rz(-1.3874929) q[0];
sx q[0];
rz(-0.36780372) q[0];
rz(-1.6416719) q[1];
sx q[1];
rz(-0.92253128) q[1];
sx q[1];
rz(-2.0466764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.80847) q[0];
sx q[0];
rz(-2.5980691) q[0];
sx q[0];
rz(-2.031206) q[0];
x q[1];
rz(-1.7000654) q[2];
sx q[2];
rz(-1.0241) q[2];
sx q[2];
rz(1.136029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75083878) q[1];
sx q[1];
rz(-1.3350643) q[1];
sx q[1];
rz(-1.5890934) q[1];
rz(-1.1672375) q[3];
sx q[3];
rz(-1.5903683) q[3];
sx q[3];
rz(0.38364601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0423923) q[2];
sx q[2];
rz(-1.7902057) q[2];
sx q[2];
rz(3.0858827) q[2];
rz(0.67692155) q[3];
sx q[3];
rz(-0.61546314) q[3];
sx q[3];
rz(0.87853471) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95213503) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(-2.1761555) q[0];
rz(1.8769439) q[1];
sx q[1];
rz(-2.3955884) q[1];
sx q[1];
rz(-2.8146578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345987) q[0];
sx q[0];
rz(-0.36834684) q[0];
sx q[0];
rz(-1.8828431) q[0];
rz(-pi) q[1];
rz(-0.44595309) q[2];
sx q[2];
rz(-1.5610245) q[2];
sx q[2];
rz(0.31644145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4591229) q[1];
sx q[1];
rz(-1.7791479) q[1];
sx q[1];
rz(-0.18110593) q[1];
rz(-pi) q[2];
rz(0.14029221) q[3];
sx q[3];
rz(-0.45278087) q[3];
sx q[3];
rz(1.0990395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7586729) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(-2.3455589) q[2];
rz(3.054079) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(-0.93761888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336808) q[0];
sx q[0];
rz(-1.7681363) q[0];
sx q[0];
rz(2.8117836) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-0.70711702) q[1];
sx q[1];
rz(-2.0475533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6716794) q[0];
sx q[0];
rz(-1.6080583) q[0];
sx q[0];
rz(0.47084634) q[0];
x q[1];
rz(-1.173063) q[2];
sx q[2];
rz(-1.5088044) q[2];
sx q[2];
rz(0.50794377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0713404) q[1];
sx q[1];
rz(-1.0254745) q[1];
sx q[1];
rz(1.7252183) q[1];
rz(-pi) q[2];
rz(-2.002203) q[3];
sx q[3];
rz(-1.5346705) q[3];
sx q[3];
rz(-0.86078913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0988934) q[2];
sx q[2];
rz(-0.77184474) q[2];
sx q[2];
rz(1.137255) q[2];
rz(-0.62816652) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(-2.1342733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657418) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(-2.0487336) q[0];
rz(1.2650371) q[1];
sx q[1];
rz(-1.3053514) q[1];
sx q[1];
rz(1.0337894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7395679) q[0];
sx q[0];
rz(-1.0453859) q[0];
sx q[0];
rz(2.5891586) q[0];
x q[1];
rz(-0.7598147) q[2];
sx q[2];
rz(-1.3926448) q[2];
sx q[2];
rz(2.9896328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6792843) q[1];
sx q[1];
rz(-2.1888715) q[1];
sx q[1];
rz(-2.0831773) q[1];
x q[2];
rz(1.4321253) q[3];
sx q[3];
rz(-1.6059173) q[3];
sx q[3];
rz(-0.39681527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8826302) q[2];
sx q[2];
rz(-2.4179103) q[2];
sx q[2];
rz(-0.21931973) q[2];
rz(2.3389881) q[3];
sx q[3];
rz(-1.8290627) q[3];
sx q[3];
rz(-0.32138166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5103067) q[0];
sx q[0];
rz(-1.8230556) q[0];
sx q[0];
rz(-2.0429116) q[0];
rz(2.9639099) q[1];
sx q[1];
rz(-1.0164574) q[1];
sx q[1];
rz(-0.16998092) q[1];
rz(0.06391025) q[2];
sx q[2];
rz(-2.2722681) q[2];
sx q[2];
rz(2.0427335) q[2];
rz(2.1228409) q[3];
sx q[3];
rz(-2.9186747) q[3];
sx q[3];
rz(-1.0001343) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
