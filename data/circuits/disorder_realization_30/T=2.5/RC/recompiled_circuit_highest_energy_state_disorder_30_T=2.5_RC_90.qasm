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
rz(0.27424115) q[0];
sx q[0];
rz(3.5477717) q[0];
sx q[0];
rz(9.6532333) q[0];
rz(-2.7574381) q[1];
sx q[1];
rz(-1.6698807) q[1];
sx q[1];
rz(2.6928112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327553) q[0];
sx q[0];
rz(-1.7102557) q[0];
sx q[0];
rz(2.7701401) q[0];
rz(-0.12334646) q[2];
sx q[2];
rz(-1.984086) q[2];
sx q[2];
rz(1.9537499) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.083605974) q[1];
sx q[1];
rz(-0.068002105) q[1];
sx q[1];
rz(1.2108832) q[1];
rz(-2.9308546) q[3];
sx q[3];
rz(-0.02861261) q[3];
sx q[3];
rz(-0.60625854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3936798) q[2];
sx q[2];
rz(-2.9069275) q[2];
sx q[2];
rz(-2.5133384) q[2];
rz(-1.0660394) q[3];
sx q[3];
rz(-0.99848905) q[3];
sx q[3];
rz(2.0121241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5265441) q[0];
sx q[0];
rz(-2.5874309) q[0];
sx q[0];
rz(-0.91914415) q[0];
rz(0.17658405) q[1];
sx q[1];
rz(-1.6932026) q[1];
sx q[1];
rz(-0.44763705) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.823671) q[0];
sx q[0];
rz(-1.5783675) q[0];
sx q[0];
rz(1.5764144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8488919) q[2];
sx q[2];
rz(-1.9534885) q[2];
sx q[2];
rz(-1.1130321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3474975) q[1];
sx q[1];
rz(-1.7298051) q[1];
sx q[1];
rz(0.68557941) q[1];
rz(-pi) q[2];
rz(1.0466879) q[3];
sx q[3];
rz(-1.1231622) q[3];
sx q[3];
rz(1.4967048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8944019) q[2];
sx q[2];
rz(-1.2195419) q[2];
sx q[2];
rz(2.2578237) q[2];
rz(2.4023458) q[3];
sx q[3];
rz(-1.2284307) q[3];
sx q[3];
rz(-1.8993186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833375) q[0];
sx q[0];
rz(-2.1012335) q[0];
sx q[0];
rz(-3.1357159) q[0];
rz(2.8351496) q[1];
sx q[1];
rz(-1.076661) q[1];
sx q[1];
rz(-1.8773016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731544) q[0];
sx q[0];
rz(-1.1958164) q[0];
sx q[0];
rz(1.2545579) q[0];
rz(1.4219079) q[2];
sx q[2];
rz(-1.5369986) q[2];
sx q[2];
rz(-1.4567018) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78759277) q[1];
sx q[1];
rz(-1.3718604) q[1];
sx q[1];
rz(-2.5721694) q[1];
rz(-0.98250947) q[3];
sx q[3];
rz(-2.4949565) q[3];
sx q[3];
rz(-1.2842632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3197202) q[2];
sx q[2];
rz(-1.3786653) q[2];
sx q[2];
rz(-2.482448) q[2];
rz(1.2705196) q[3];
sx q[3];
rz(-0.29265413) q[3];
sx q[3];
rz(1.6856153) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6419411) q[0];
sx q[0];
rz(-0.15436509) q[0];
sx q[0];
rz(0.9374215) q[0];
rz(-1.591466) q[1];
sx q[1];
rz(-2.4731686) q[1];
sx q[1];
rz(-2.392427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2523152) q[0];
sx q[0];
rz(-1.2752644) q[0];
sx q[0];
rz(-0.79177971) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5116483) q[2];
sx q[2];
rz(-2.1376746) q[2];
sx q[2];
rz(0.34467372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8038102) q[1];
sx q[1];
rz(-2.446108) q[1];
sx q[1];
rz(1.8799604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9177298) q[3];
sx q[3];
rz(-0.77896691) q[3];
sx q[3];
rz(1.1457486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59430993) q[2];
sx q[2];
rz(-2.3197033) q[2];
sx q[2];
rz(1.8853356) q[2];
rz(-2.3516583) q[3];
sx q[3];
rz(-0.94090763) q[3];
sx q[3];
rz(-1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.054852) q[0];
sx q[0];
rz(-2.5331443) q[0];
sx q[0];
rz(0.21859455) q[0];
rz(-2.356148) q[1];
sx q[1];
rz(-1.1050858) q[1];
sx q[1];
rz(-1.8484263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620813) q[0];
sx q[0];
rz(-1.1089227) q[0];
sx q[0];
rz(2.1681468) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3082383) q[2];
sx q[2];
rz(-1.6505604) q[2];
sx q[2];
rz(-2.3924116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47673273) q[1];
sx q[1];
rz(-0.97366758) q[1];
sx q[1];
rz(0.67311153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3985211) q[3];
sx q[3];
rz(-0.2187116) q[3];
sx q[3];
rz(-0.89994535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.834632) q[2];
sx q[2];
rz(-2.0768276) q[2];
sx q[2];
rz(1.5131697) q[2];
rz(0.73897922) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(-1.6966604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4479248) q[0];
sx q[0];
rz(-1.6532927) q[0];
sx q[0];
rz(-2.779575) q[0];
rz(-2.0118878) q[1];
sx q[1];
rz(-1.874322) q[1];
sx q[1];
rz(-0.76090181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1429962) q[0];
sx q[0];
rz(-0.79299295) q[0];
sx q[0];
rz(-2.8727813) q[0];
x q[1];
rz(-1.1298112) q[2];
sx q[2];
rz(-0.79409696) q[2];
sx q[2];
rz(1.5979021) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2471091) q[1];
sx q[1];
rz(-1.037408) q[1];
sx q[1];
rz(-0.24786149) q[1];
rz(-pi) q[2];
rz(2.6057556) q[3];
sx q[3];
rz(-0.69147516) q[3];
sx q[3];
rz(-2.7590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0249944) q[2];
sx q[2];
rz(-0.28433871) q[2];
sx q[2];
rz(-1.1762168) q[2];
rz(-1.0233915) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(-0.37337506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244331) q[0];
sx q[0];
rz(-2.7130782) q[0];
sx q[0];
rz(1.3065216) q[0];
rz(-0.10082968) q[1];
sx q[1];
rz(-1.3918624) q[1];
sx q[1];
rz(-0.94591013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27594179) q[0];
sx q[0];
rz(-1.2667333) q[0];
sx q[0];
rz(-1.7070387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77581279) q[2];
sx q[2];
rz(-0.6686506) q[2];
sx q[2];
rz(0.60912998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7335763) q[1];
sx q[1];
rz(-2.9936572) q[1];
sx q[1];
rz(-2.7285517) q[1];
rz(-2.2846388) q[3];
sx q[3];
rz(-1.3009239) q[3];
sx q[3];
rz(-1.2601579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0239608) q[2];
sx q[2];
rz(-1.6098216) q[2];
sx q[2];
rz(-2.7839933) q[2];
rz(0.050497342) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(-2.7549506) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015775) q[0];
sx q[0];
rz(-1.3192663) q[0];
sx q[0];
rz(2.5460119) q[0];
rz(-1.2254084) q[1];
sx q[1];
rz(-1.775368) q[1];
sx q[1];
rz(2.7092686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29994592) q[0];
sx q[0];
rz(-0.52021813) q[0];
sx q[0];
rz(0.33824091) q[0];
rz(-pi) q[1];
rz(1.7244699) q[2];
sx q[2];
rz(-1.6343079) q[2];
sx q[2];
rz(-0.85086289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22307587) q[1];
sx q[1];
rz(-2.1380414) q[1];
sx q[1];
rz(1.5722797) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0995977) q[3];
sx q[3];
rz(-0.75529683) q[3];
sx q[3];
rz(0.55213682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7633742) q[2];
sx q[2];
rz(-1.32064) q[2];
sx q[2];
rz(-1.0279559) q[2];
rz(1.6592735) q[3];
sx q[3];
rz(-1.6927398) q[3];
sx q[3];
rz(2.3894374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1886275) q[0];
sx q[0];
rz(-1.9099706) q[0];
sx q[0];
rz(0.46284437) q[0];
rz(-1.9075182) q[1];
sx q[1];
rz(-1.2214829) q[1];
sx q[1];
rz(5*pi/8) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5355204) q[0];
sx q[0];
rz(-1.3957481) q[0];
sx q[0];
rz(0.025786215) q[0];
x q[1];
rz(-0.47204702) q[2];
sx q[2];
rz(-0.85523048) q[2];
sx q[2];
rz(-0.89821076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8527216) q[1];
sx q[1];
rz(-0.72119547) q[1];
sx q[1];
rz(-1.3836963) q[1];
rz(1.6912446) q[3];
sx q[3];
rz(-1.9102194) q[3];
sx q[3];
rz(1.898693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50596333) q[2];
sx q[2];
rz(-1.3869945) q[2];
sx q[2];
rz(-2.3734234) q[2];
rz(0.65028894) q[3];
sx q[3];
rz(-0.30935973) q[3];
sx q[3];
rz(0.55997854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77785093) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(1.1458696) q[0];
rz(-2.5762985) q[1];
sx q[1];
rz(-0.71563131) q[1];
sx q[1];
rz(-0.40564767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7956896) q[0];
sx q[0];
rz(-2.3064605) q[0];
sx q[0];
rz(2.664734) q[0];
rz(-pi) q[1];
rz(-0.94735165) q[2];
sx q[2];
rz(-2.6877786) q[2];
sx q[2];
rz(2.712534) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1310869) q[1];
sx q[1];
rz(-1.0684135) q[1];
sx q[1];
rz(-1.2817863) q[1];
rz(-pi) q[2];
rz(-0.43934699) q[3];
sx q[3];
rz(-1.4733107) q[3];
sx q[3];
rz(2.0078703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.046935) q[2];
sx q[2];
rz(-1.9104702) q[2];
sx q[2];
rz(-2.6222353) q[2];
rz(0.46118593) q[3];
sx q[3];
rz(-1.0591732) q[3];
sx q[3];
rz(0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984579) q[0];
sx q[0];
rz(-1.5512137) q[0];
sx q[0];
rz(2.2358538) q[0];
rz(-1.657919) q[1];
sx q[1];
rz(-1.1653733) q[1];
sx q[1];
rz(1.9001874) q[1];
rz(-2.9979669) q[2];
sx q[2];
rz(-0.96114071) q[2];
sx q[2];
rz(-3.0428934) q[2];
rz(2.3158793) q[3];
sx q[3];
rz(-2.5057247) q[3];
sx q[3];
rz(2.417939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
