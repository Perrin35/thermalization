OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(2.2428089) q[0];
sx q[0];
rz(8.1221683) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(-1.5136493) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2500656) q[0];
sx q[0];
rz(-2.6882049) q[0];
sx q[0];
rz(2.2946623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3816125) q[2];
sx q[2];
rz(-0.70356762) q[2];
sx q[2];
rz(2.4576996) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4931655) q[1];
sx q[1];
rz(-1.2305639) q[1];
sx q[1];
rz(3.0837584) q[1];
rz(-pi) q[2];
rz(-2.1964425) q[3];
sx q[3];
rz(-1.8813475) q[3];
sx q[3];
rz(0.20598447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0836432) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(0.38457125) q[2];
rz(-2.7390076) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(-0.80818498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.758051) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(-0.59355271) q[0];
rz(-1.3213762) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(-0.039332565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101321) q[0];
sx q[0];
rz(-1.1173741) q[0];
sx q[0];
rz(2.7260145) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63644511) q[2];
sx q[2];
rz(-1.6628254) q[2];
sx q[2];
rz(-0.68816371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96728863) q[1];
sx q[1];
rz(-1.5212465) q[1];
sx q[1];
rz(0.51731717) q[1];
rz(1.5402628) q[3];
sx q[3];
rz(-1.021592) q[3];
sx q[3];
rz(0.98992482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3010657) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(0.62704101) q[2];
rz(2.7035233) q[3];
sx q[3];
rz(-1.6431243) q[3];
sx q[3];
rz(1.1367249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7196734) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(-1.3713974) q[0];
rz(0.60945359) q[1];
sx q[1];
rz(-0.89027673) q[1];
sx q[1];
rz(1.9166463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44642513) q[0];
sx q[0];
rz(-1.7621627) q[0];
sx q[0];
rz(-1.6686977) q[0];
x q[1];
rz(1.0880427) q[2];
sx q[2];
rz(-0.73763371) q[2];
sx q[2];
rz(-0.18407735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6869192) q[1];
sx q[1];
rz(-2.6489566) q[1];
sx q[1];
rz(-0.83322816) q[1];
rz(-pi) q[2];
rz(-0.87525778) q[3];
sx q[3];
rz(-1.7376383) q[3];
sx q[3];
rz(2.5341906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3043392) q[2];
sx q[2];
rz(-2.2746268) q[2];
sx q[2];
rz(0.85401094) q[2];
rz(-0.52418661) q[3];
sx q[3];
rz(-1.6322522) q[3];
sx q[3];
rz(-0.9700276) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22685856) q[0];
sx q[0];
rz(-0.15758841) q[0];
sx q[0];
rz(-2.3478813) q[0];
rz(-3.1397676) q[1];
sx q[1];
rz(-2.5853214) q[1];
sx q[1];
rz(2.3307641) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49634078) q[0];
sx q[0];
rz(-0.84028572) q[0];
sx q[0];
rz(-2.8964983) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0433719) q[2];
sx q[2];
rz(-0.7458936) q[2];
sx q[2];
rz(-0.30442849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6135679) q[1];
sx q[1];
rz(-0.9903637) q[1];
sx q[1];
rz(0.40596647) q[1];
rz(0.14489095) q[3];
sx q[3];
rz(-1.1456744) q[3];
sx q[3];
rz(1.8505877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4578555) q[2];
sx q[2];
rz(-0.33317864) q[2];
sx q[2];
rz(-1.4972756) q[2];
rz(1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(1.3076521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7233647) q[0];
sx q[0];
rz(-1.4606553) q[0];
sx q[0];
rz(2.8826868) q[0];
rz(-3.018191) q[1];
sx q[1];
rz(-0.63107189) q[1];
sx q[1];
rz(-0.48615989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82993342) q[0];
sx q[0];
rz(-1.0898542) q[0];
sx q[0];
rz(-0.89998683) q[0];
rz(-pi) q[1];
rz(-0.27433594) q[2];
sx q[2];
rz(-2.6736709) q[2];
sx q[2];
rz(-2.5578424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3868931) q[1];
sx q[1];
rz(-0.4122977) q[1];
sx q[1];
rz(-1.0846433) q[1];
rz(-pi) q[2];
rz(0.93323243) q[3];
sx q[3];
rz(-2.3321816) q[3];
sx q[3];
rz(0.903086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2316124) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.3952338) q[2];
rz(2.1868271) q[3];
sx q[3];
rz(-1.3937817) q[3];
sx q[3];
rz(-2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2024277) q[0];
sx q[0];
rz(-1.2347777) q[0];
sx q[0];
rz(-0.760461) q[0];
rz(-0.83433926) q[1];
sx q[1];
rz(-1.1015247) q[1];
sx q[1];
rz(1.1044501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4882979) q[0];
sx q[0];
rz(-2.1135931) q[0];
sx q[0];
rz(0.25570583) q[0];
rz(-2.274674) q[2];
sx q[2];
rz(-0.50627497) q[2];
sx q[2];
rz(-2.3828854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.901746) q[1];
sx q[1];
rz(-0.86565986) q[1];
sx q[1];
rz(0.88259952) q[1];
rz(1.1105632) q[3];
sx q[3];
rz(-0.62556872) q[3];
sx q[3];
rz(0.21504185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96222782) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(-0.041570138) q[2];
rz(-2.9465594) q[3];
sx q[3];
rz(-0.2427559) q[3];
sx q[3];
rz(-0.78708762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389909) q[0];
sx q[0];
rz(-2.2881303) q[0];
sx q[0];
rz(-2.3882197) q[0];
rz(-2.0808749) q[1];
sx q[1];
rz(-2.3371688) q[1];
sx q[1];
rz(-2.9341968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78847892) q[0];
sx q[0];
rz(-0.99035701) q[0];
sx q[0];
rz(-0.81932318) q[0];
rz(-pi) q[1];
rz(-1.8600402) q[2];
sx q[2];
rz(-1.1221894) q[2];
sx q[2];
rz(-1.9236652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2826536) q[1];
sx q[1];
rz(-2.5368854) q[1];
sx q[1];
rz(2.7207123) q[1];
rz(-1.5630521) q[3];
sx q[3];
rz(-1.6098579) q[3];
sx q[3];
rz(-2.9970053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14219001) q[2];
sx q[2];
rz(-1.157607) q[2];
sx q[2];
rz(-3.0863777) q[2];
rz(2.2986872) q[3];
sx q[3];
rz(-0.75777811) q[3];
sx q[3];
rz(0.88111669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.76999369) q[0];
sx q[0];
rz(-2.433625) q[0];
sx q[0];
rz(-1.9644894) q[0];
rz(-2.2616995) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(0.060861977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5755479) q[0];
sx q[0];
rz(-1.6948997) q[0];
sx q[0];
rz(3.0461765) q[0];
rz(2.6598236) q[2];
sx q[2];
rz(-1.2854115) q[2];
sx q[2];
rz(0.18295675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5338075) q[1];
sx q[1];
rz(-2.5609511) q[1];
sx q[1];
rz(-2.9279079) q[1];
x q[2];
rz(1.8350321) q[3];
sx q[3];
rz(-1.4995534) q[3];
sx q[3];
rz(2.3087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0941144) q[2];
sx q[2];
rz(-2.5884509) q[2];
sx q[2];
rz(-1.4609569) q[2];
rz(2.285752) q[3];
sx q[3];
rz(-2.1688921) q[3];
sx q[3];
rz(-2.262825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0395373) q[0];
sx q[0];
rz(-2.2886031) q[0];
sx q[0];
rz(-0.0041740388) q[0];
rz(1.8904842) q[1];
sx q[1];
rz(-0.1336385) q[1];
sx q[1];
rz(0.68152308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96723667) q[0];
sx q[0];
rz(-0.78034329) q[0];
sx q[0];
rz(-1.1514949) q[0];
rz(-0.60891843) q[2];
sx q[2];
rz(-1.0002197) q[2];
sx q[2];
rz(-1.0341687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79922593) q[1];
sx q[1];
rz(-1.7569506) q[1];
sx q[1];
rz(0.28719814) q[1];
rz(2.8130346) q[3];
sx q[3];
rz(-0.78379455) q[3];
sx q[3];
rz(-2.4155922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.556813) q[2];
sx q[2];
rz(-0.4759554) q[2];
sx q[2];
rz(-1.7468096) q[2];
rz(-2.7252588) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(0.39286119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3861179) q[0];
sx q[0];
rz(-2.8135354) q[0];
sx q[0];
rz(2.1395444) q[0];
rz(-2.877032) q[1];
sx q[1];
rz(-1.325565) q[1];
sx q[1];
rz(1.6511668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39367657) q[0];
sx q[0];
rz(-1.0320469) q[0];
sx q[0];
rz(1.7231621) q[0];
rz(-pi) q[1];
rz(1.7298431) q[2];
sx q[2];
rz(-2.3727086) q[2];
sx q[2];
rz(-1.3231708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7539053) q[1];
sx q[1];
rz(-1.6847982) q[1];
sx q[1];
rz(3.0579159) q[1];
x q[2];
rz(-2.5633859) q[3];
sx q[3];
rz(-3.0319408) q[3];
sx q[3];
rz(0.46441098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.010290535) q[2];
sx q[2];
rz(-1.7346069) q[2];
sx q[2];
rz(-1.7245801) q[2];
rz(-0.33306444) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(1.7634348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604816) q[0];
sx q[0];
rz(-1.2662553) q[0];
sx q[0];
rz(1.8929831) q[0];
rz(0.026451182) q[1];
sx q[1];
rz(-1.6117663) q[1];
sx q[1];
rz(1.5996006) q[1];
rz(-1.1556861) q[2];
sx q[2];
rz(-0.90183707) q[2];
sx q[2];
rz(3.1169008) q[2];
rz(0.61458434) q[3];
sx q[3];
rz(-1.8329704) q[3];
sx q[3];
rz(-0.80154673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
