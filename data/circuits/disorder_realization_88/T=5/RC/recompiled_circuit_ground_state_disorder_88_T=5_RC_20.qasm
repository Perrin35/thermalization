OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.705536) q[0];
sx q[0];
rz(-2.2826865) q[0];
sx q[0];
rz(-0.98281759) q[0];
rz(1.1801899) q[1];
sx q[1];
rz(-2.4153914) q[1];
sx q[1];
rz(2.243431) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58098823) q[0];
sx q[0];
rz(-0.33412558) q[0];
sx q[0];
rz(0.32291205) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9299081) q[2];
sx q[2];
rz(-1.058033) q[2];
sx q[2];
rz(2.7732234) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8232728) q[1];
sx q[1];
rz(-1.987769) q[1];
sx q[1];
rz(0.61064536) q[1];
rz(-pi) q[2];
rz(-1.6706561) q[3];
sx q[3];
rz(-1.7152399) q[3];
sx q[3];
rz(0.58510548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1232221) q[2];
sx q[2];
rz(-1.8223338) q[2];
sx q[2];
rz(-2.3940864) q[2];
rz(1.2627164) q[3];
sx q[3];
rz(-1.5371753) q[3];
sx q[3];
rz(-0.023716299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338407) q[0];
sx q[0];
rz(-0.045578651) q[0];
sx q[0];
rz(-2.0763092) q[0];
rz(0.65159687) q[1];
sx q[1];
rz(-2.7862796) q[1];
sx q[1];
rz(1.5078872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3443073) q[0];
sx q[0];
rz(-2.4977614) q[0];
sx q[0];
rz(2.0876398) q[0];
x q[1];
rz(-2.5152757) q[2];
sx q[2];
rz(-0.56506116) q[2];
sx q[2];
rz(2.9941699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44017556) q[1];
sx q[1];
rz(-1.9721627) q[1];
sx q[1];
rz(1.0640776) q[1];
rz(1.6117068) q[3];
sx q[3];
rz(-0.61167756) q[3];
sx q[3];
rz(-3.1168695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4773341) q[2];
sx q[2];
rz(-1.6691672) q[2];
sx q[2];
rz(0.092546917) q[2];
rz(-0.44068286) q[3];
sx q[3];
rz(-0.40658545) q[3];
sx q[3];
rz(1.1455166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6630845) q[0];
sx q[0];
rz(-2.9496851) q[0];
sx q[0];
rz(-1.9051911) q[0];
rz(2.9036486) q[1];
sx q[1];
rz(-1.4711719) q[1];
sx q[1];
rz(2.1740289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2264573) q[0];
sx q[0];
rz(-2.1564556) q[0];
sx q[0];
rz(2.4047635) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5406303) q[2];
sx q[2];
rz(-1.8701359) q[2];
sx q[2];
rz(2.9082852) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1246207) q[1];
sx q[1];
rz(-1.3836814) q[1];
sx q[1];
rz(1.8797423) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4400515) q[3];
sx q[3];
rz(-0.41660158) q[3];
sx q[3];
rz(2.1611913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48270109) q[2];
sx q[2];
rz(-2.7693558) q[2];
sx q[2];
rz(2.2230395) q[2];
rz(2.3679768) q[3];
sx q[3];
rz(-0.88671237) q[3];
sx q[3];
rz(-0.95019597) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2660265) q[0];
sx q[0];
rz(-2.4022864) q[0];
sx q[0];
rz(2.8493122) q[0];
rz(2.267011) q[1];
sx q[1];
rz(-2.5129109) q[1];
sx q[1];
rz(0.94327092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90664141) q[0];
sx q[0];
rz(-0.72777339) q[0];
sx q[0];
rz(2.0208738) q[0];
x q[1];
rz(-2.9246632) q[2];
sx q[2];
rz(-1.1116905) q[2];
sx q[2];
rz(2.2709059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5041793) q[1];
sx q[1];
rz(-2.4054745) q[1];
sx q[1];
rz(2.0448684) q[1];
rz(-pi) q[2];
rz(-0.28369035) q[3];
sx q[3];
rz(-2.1296888) q[3];
sx q[3];
rz(-2.7070466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6167407) q[2];
sx q[2];
rz(-2.0276232) q[2];
sx q[2];
rz(0.88391602) q[2];
rz(2.0058477) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(-0.69042027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32490548) q[0];
sx q[0];
rz(-0.08520928) q[0];
sx q[0];
rz(2.9499522) q[0];
rz(1.2958255) q[1];
sx q[1];
rz(-1.192966) q[1];
sx q[1];
rz(-2.6909018) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0518351) q[0];
sx q[0];
rz(-2.20443) q[0];
sx q[0];
rz(3.0283338) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8692372) q[2];
sx q[2];
rz(-1.5030454) q[2];
sx q[2];
rz(0.94558796) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5366392) q[1];
sx q[1];
rz(-0.99483314) q[1];
sx q[1];
rz(0.085809068) q[1];
x q[2];
rz(-2.8487318) q[3];
sx q[3];
rz(-1.5080875) q[3];
sx q[3];
rz(0.69578275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4638764) q[2];
sx q[2];
rz(-1.8666942) q[2];
sx q[2];
rz(2.1837168) q[2];
rz(3.0887582) q[3];
sx q[3];
rz(-2.5962679) q[3];
sx q[3];
rz(0.42832819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65474725) q[0];
sx q[0];
rz(-2.0992278) q[0];
sx q[0];
rz(0.47252193) q[0];
rz(0.76332244) q[1];
sx q[1];
rz(-0.7205874) q[1];
sx q[1];
rz(2.8514013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4241353) q[0];
sx q[0];
rz(-2.4032434) q[0];
sx q[0];
rz(-2.8108869) q[0];
x q[1];
rz(1.0792208) q[2];
sx q[2];
rz(-0.60814684) q[2];
sx q[2];
rz(1.3754932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.063263254) q[1];
sx q[1];
rz(-1.4681879) q[1];
sx q[1];
rz(2.4797477) q[1];
rz(-pi) q[2];
rz(0.17866023) q[3];
sx q[3];
rz(-1.7669441) q[3];
sx q[3];
rz(2.2771108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7758238) q[2];
sx q[2];
rz(-0.86161986) q[2];
sx q[2];
rz(-0.59930581) q[2];
rz(1.722909) q[3];
sx q[3];
rz(-1.8214046) q[3];
sx q[3];
rz(0.9532218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80239427) q[0];
sx q[0];
rz(-1.9559487) q[0];
sx q[0];
rz(1.8307357) q[0];
rz(-1.5015548) q[1];
sx q[1];
rz(-2.7412667) q[1];
sx q[1];
rz(-0.60360533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045403584) q[0];
sx q[0];
rz(-1.9902162) q[0];
sx q[0];
rz(2.9950401) q[0];
rz(-2.7947803) q[2];
sx q[2];
rz(-0.57935994) q[2];
sx q[2];
rz(2.8413136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1511076) q[1];
sx q[1];
rz(-1.8112887) q[1];
sx q[1];
rz(-0.069995306) q[1];
rz(-pi) q[2];
rz(2.9333618) q[3];
sx q[3];
rz(-0.75246598) q[3];
sx q[3];
rz(2.3104639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8169868) q[2];
sx q[2];
rz(-2.2548455) q[2];
sx q[2];
rz(2.9465607) q[2];
rz(2.4845691) q[3];
sx q[3];
rz(-1.421509) q[3];
sx q[3];
rz(0.11369625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103771) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(-3.069416) q[0];
rz(-2.3585034) q[1];
sx q[1];
rz(-2.5091645) q[1];
sx q[1];
rz(2.2236688) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1872301) q[0];
sx q[0];
rz(-1.7779121) q[0];
sx q[0];
rz(2.6595479) q[0];
rz(-pi) q[1];
rz(-1.8391219) q[2];
sx q[2];
rz(-0.96728281) q[2];
sx q[2];
rz(2.8203815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0705441) q[1];
sx q[1];
rz(-1.8080825) q[1];
sx q[1];
rz(-0.29233934) q[1];
x q[2];
rz(1.9959171) q[3];
sx q[3];
rz(-1.3367869) q[3];
sx q[3];
rz(1.0581512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5607295) q[2];
sx q[2];
rz(-1.9274638) q[2];
sx q[2];
rz(1.6403991) q[2];
rz(-0.8864657) q[3];
sx q[3];
rz(-1.3366046) q[3];
sx q[3];
rz(-3.0361573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5378872) q[0];
sx q[0];
rz(-2.7925346) q[0];
sx q[0];
rz(0.50759298) q[0];
rz(2.4137068) q[1];
sx q[1];
rz(-1.6517086) q[1];
sx q[1];
rz(-1.1061888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78680404) q[0];
sx q[0];
rz(-1.3193325) q[0];
sx q[0];
rz(-1.2261816) q[0];
rz(-pi) q[1];
rz(-0.35850766) q[2];
sx q[2];
rz(-1.7353183) q[2];
sx q[2];
rz(1.2929971) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.008279) q[1];
sx q[1];
rz(-2.6402557) q[1];
sx q[1];
rz(0.8322844) q[1];
rz(-pi) q[2];
rz(-0.9000128) q[3];
sx q[3];
rz(-2.4379895) q[3];
sx q[3];
rz(0.37868628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20405208) q[2];
sx q[2];
rz(-1.1037408) q[2];
sx q[2];
rz(0.44720116) q[2];
rz(-1.7043381) q[3];
sx q[3];
rz(-2.3279326) q[3];
sx q[3];
rz(-1.6566431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5337885) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(-0.5150038) q[0];
rz(-1.1732514) q[1];
sx q[1];
rz(-1.2898022) q[1];
sx q[1];
rz(-2.692093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052348) q[0];
sx q[0];
rz(-0.32669386) q[0];
sx q[0];
rz(1.5004083) q[0];
rz(-pi) q[1];
rz(-0.88926824) q[2];
sx q[2];
rz(-0.73822108) q[2];
sx q[2];
rz(-0.29730931) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.044552) q[1];
sx q[1];
rz(-0.97787338) q[1];
sx q[1];
rz(-2.135599) q[1];
rz(-2.2902255) q[3];
sx q[3];
rz(-1.3847305) q[3];
sx q[3];
rz(2.1894941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0939193) q[2];
sx q[2];
rz(-0.27457044) q[2];
sx q[2];
rz(-0.74335113) q[2];
rz(-0.36176935) q[3];
sx q[3];
rz(-2.1105364) q[3];
sx q[3];
rz(0.37775347) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7985228) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(0.15380225) q[1];
sx q[1];
rz(-1.8228795) q[1];
sx q[1];
rz(0.22620329) q[1];
rz(-3.0085887) q[2];
sx q[2];
rz(-1.8756744) q[2];
sx q[2];
rz(-1.6185624) q[2];
rz(-0.50688898) q[3];
sx q[3];
rz(-2.4573368) q[3];
sx q[3];
rz(-1.5507207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
