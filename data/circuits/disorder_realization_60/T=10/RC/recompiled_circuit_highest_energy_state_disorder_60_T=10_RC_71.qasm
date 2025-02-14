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
rz(0.76437104) q[0];
sx q[0];
rz(1.8005014) q[0];
sx q[0];
rz(10.330248) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(1.5009872) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78202263) q[0];
sx q[0];
rz(-1.3728314) q[0];
sx q[0];
rz(-0.85967608) q[0];
rz(-pi) q[1];
rz(-2.0101648) q[2];
sx q[2];
rz(-1.7305534) q[2];
sx q[2];
rz(-1.8607832) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5792306) q[1];
sx q[1];
rz(-3.0931614) q[1];
sx q[1];
rz(-2.6391451) q[1];
rz(-pi) q[2];
rz(1.0085868) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(-2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(-1.9979075) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2442653) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(-0.51189297) q[0];
rz(1.5533252) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(-1.1221251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0034135) q[0];
sx q[0];
rz(-1.8818568) q[0];
sx q[0];
rz(-1.6775234) q[0];
x q[1];
rz(1.4830515) q[2];
sx q[2];
rz(-0.43102396) q[2];
sx q[2];
rz(-1.6602248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2832239) q[1];
sx q[1];
rz(-2.5852381) q[1];
sx q[1];
rz(-2.2409641) q[1];
rz(-pi) q[2];
rz(0.53776017) q[3];
sx q[3];
rz(-2.3590292) q[3];
sx q[3];
rz(-1.3937221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3731709) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-2.6785417) q[2];
rz(2.566346) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082044) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(-0.43011618) q[0];
rz(-2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(2.5045085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642536) q[0];
sx q[0];
rz(-1.5298944) q[0];
sx q[0];
rz(-0.014057191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78539679) q[2];
sx q[2];
rz(-1.5893418) q[2];
sx q[2];
rz(-2.3634499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8608537) q[1];
sx q[1];
rz(-0.85608427) q[1];
sx q[1];
rz(1.2398941) q[1];
rz(-pi) q[2];
rz(-0.62406059) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(-2.000119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-0.72506881) q[2];
rz(-1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8722039) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(1.4440906) q[0];
rz(0.048642453) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(0.28894249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699997) q[0];
sx q[0];
rz(-0.44737383) q[0];
sx q[0];
rz(-2.1567221) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.028325105) q[2];
sx q[2];
rz(-1.7212756) q[2];
sx q[2];
rz(-0.17720824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0045341) q[1];
sx q[1];
rz(-1.6797795) q[1];
sx q[1];
rz(1.7458899) q[1];
rz(0.66529243) q[3];
sx q[3];
rz(-1.6136618) q[3];
sx q[3];
rz(0.24584082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.284953) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(2.7613769) q[2];
rz(-2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(-1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3087092) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(-1.0943476) q[0];
rz(2.265918) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(1.099115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581061) q[0];
sx q[0];
rz(-1.1855584) q[0];
sx q[0];
rz(1.7190821) q[0];
x q[1];
rz(0.78754707) q[2];
sx q[2];
rz(-0.468245) q[2];
sx q[2];
rz(2.5494273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5531986) q[1];
sx q[1];
rz(-1.0244245) q[1];
sx q[1];
rz(-2.9562034) q[1];
rz(-pi) q[2];
rz(0.91644561) q[3];
sx q[3];
rz(-1.4438085) q[3];
sx q[3];
rz(2.6892852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(-0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(-1.2716768) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.6082825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23841183) q[0];
sx q[0];
rz(-3.0738214) q[0];
sx q[0];
rz(1.7525903) q[0];
rz(-1.4727696) q[2];
sx q[2];
rz(-1.5614126) q[2];
sx q[2];
rz(-2.2414152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.530513) q[1];
sx q[1];
rz(-2.8332498) q[1];
sx q[1];
rz(-0.18402305) q[1];
x q[2];
rz(-0.23171429) q[3];
sx q[3];
rz(-2.3149256) q[3];
sx q[3];
rz(0.8207013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61234683) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(-0.38522729) q[2];
rz(0.084608229) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(-0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(-2.7440942) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(0.0040815512) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57848141) q[0];
sx q[0];
rz(-2.3627776) q[0];
sx q[0];
rz(1.7527197) q[0];
x q[1];
rz(2.053431) q[2];
sx q[2];
rz(-0.74649278) q[2];
sx q[2];
rz(1.6255524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34985182) q[1];
sx q[1];
rz(-0.71261969) q[1];
sx q[1];
rz(0.60421677) q[1];
rz(2.2479731) q[3];
sx q[3];
rz(-1.0056408) q[3];
sx q[3];
rz(-0.86097417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6323382) q[2];
sx q[2];
rz(-1.961144) q[2];
sx q[2];
rz(2.9285367) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(-2.842438) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(-2.1655653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19071391) q[0];
sx q[0];
rz(-1.6813206) q[0];
sx q[0];
rz(0.24675639) q[0];
x q[1];
rz(0.86324228) q[2];
sx q[2];
rz(-2.9094271) q[2];
sx q[2];
rz(-1.6492998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8638226) q[1];
sx q[1];
rz(-2.7105717) q[1];
sx q[1];
rz(2.1716539) q[1];
rz(1.8194514) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(-1.4506884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(-1.3236375) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-2.8470993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2015304) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(-1.6424204) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(-1.7971719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7758862) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(-1.9802753) q[0];
rz(2.1189519) q[2];
sx q[2];
rz(-1.022199) q[2];
sx q[2];
rz(2.5238069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0186179) q[1];
sx q[1];
rz(-0.6280762) q[1];
sx q[1];
rz(-1.1372034) q[1];
rz(0.88560652) q[3];
sx q[3];
rz(-2.2937498) q[3];
sx q[3];
rz(0.9566488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.4538291) q[2];
rz(2.7135571) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(-1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(1.6499299) q[0];
rz(-2.5904169) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(0.59250441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4649749) q[0];
sx q[0];
rz(-1.3335557) q[0];
sx q[0];
rz(0.28740164) q[0];
rz(-pi) q[1];
rz(2.1001753) q[2];
sx q[2];
rz(-1.0905078) q[2];
sx q[2];
rz(-0.65435556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69018248) q[1];
sx q[1];
rz(-1.3082778) q[1];
sx q[1];
rz(-2.7464944) q[1];
rz(-pi) q[2];
rz(-1.5207401) q[3];
sx q[3];
rz(-2.069469) q[3];
sx q[3];
rz(0.071144516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42980117) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(3.0832624) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2881099) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(1.7755605) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-0.38717196) q[2];
sx q[2];
rz(-0.6091112) q[2];
sx q[2];
rz(2.1412639) q[2];
rz(0.80647918) q[3];
sx q[3];
rz(-0.69533336) q[3];
sx q[3];
rz(1.7736377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
