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
rz(0.38415456) q[1];
sx q[1];
rz(-1.471712) q[1];
sx q[1];
rz(-2.6928112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81918994) q[0];
sx q[0];
rz(-0.39562851) q[0];
sx q[0];
rz(-0.36901335) q[0];
x q[1];
rz(3.0182462) q[2];
sx q[2];
rz(-1.1575067) q[2];
sx q[2];
rz(1.1878428) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0579867) q[1];
sx q[1];
rz(-0.068002105) q[1];
sx q[1];
rz(1.9307095) q[1];
rz(-pi) q[2];
rz(-0.027979942) q[3];
sx q[3];
rz(-1.5648119) q[3];
sx q[3];
rz(2.3877092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7479129) q[2];
sx q[2];
rz(-0.23466514) q[2];
sx q[2];
rz(2.5133384) q[2];
rz(-2.0755532) q[3];
sx q[3];
rz(-2.1431036) q[3];
sx q[3];
rz(-1.1294686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5265441) q[0];
sx q[0];
rz(-0.55416179) q[0];
sx q[0];
rz(-0.91914415) q[0];
rz(-2.9650086) q[1];
sx q[1];
rz(-1.6932026) q[1];
sx q[1];
rz(-0.44763705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31792163) q[0];
sx q[0];
rz(-1.5783675) q[0];
sx q[0];
rz(-1.5651783) q[0];
rz(-pi) q[1];
rz(1.8488919) q[2];
sx q[2];
rz(-1.1881042) q[2];
sx q[2];
rz(-2.0285605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7270738) q[1];
sx q[1];
rz(-2.440732) q[1];
sx q[1];
rz(-2.8935211) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0949048) q[3];
sx q[3];
rz(-1.1231622) q[3];
sx q[3];
rz(1.6448878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8944019) q[2];
sx q[2];
rz(-1.9220507) q[2];
sx q[2];
rz(2.2578237) q[2];
rz(0.73924685) q[3];
sx q[3];
rz(-1.2284307) q[3];
sx q[3];
rz(-1.242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7582551) q[0];
sx q[0];
rz(-1.0403591) q[0];
sx q[0];
rz(3.1357159) q[0];
rz(-2.8351496) q[1];
sx q[1];
rz(-1.076661) q[1];
sx q[1];
rz(-1.2642911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684382) q[0];
sx q[0];
rz(-1.9457762) q[0];
sx q[0];
rz(1.8870348) q[0];
x q[1];
rz(3.1074171) q[2];
sx q[2];
rz(-1.4219936) q[2];
sx q[2];
rz(3.0224295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90904655) q[1];
sx q[1];
rz(-2.1276444) q[1];
sx q[1];
rz(1.3358447) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7448768) q[3];
sx q[3];
rz(-1.0457888) q[3];
sx q[3];
rz(-0.58806149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3197202) q[2];
sx q[2];
rz(-1.7629273) q[2];
sx q[2];
rz(2.482448) q[2];
rz(-1.871073) q[3];
sx q[3];
rz(-0.29265413) q[3];
sx q[3];
rz(1.6856153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4996516) q[0];
sx q[0];
rz(-2.9872276) q[0];
sx q[0];
rz(2.2041712) q[0];
rz(-1.591466) q[1];
sx q[1];
rz(-2.4731686) q[1];
sx q[1];
rz(-2.392427) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8892775) q[0];
sx q[0];
rz(-1.8663283) q[0];
sx q[0];
rz(2.3498129) q[0];
x q[1];
rz(2.5116483) q[2];
sx q[2];
rz(-1.0039181) q[2];
sx q[2];
rz(-0.34467372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73203304) q[1];
sx q[1];
rz(-0.91425843) q[1];
sx q[1];
rz(-0.2486807) q[1];
x q[2];
rz(2.375256) q[3];
sx q[3];
rz(-1.727399) q[3];
sx q[3];
rz(-0.58567724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5472827) q[2];
sx q[2];
rz(-0.82188934) q[2];
sx q[2];
rz(1.8853356) q[2];
rz(-2.3516583) q[3];
sx q[3];
rz(-0.94090763) q[3];
sx q[3];
rz(1.9739715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.054852) q[0];
sx q[0];
rz(-2.5331443) q[0];
sx q[0];
rz(-2.9229981) q[0];
rz(2.356148) q[1];
sx q[1];
rz(-1.1050858) q[1];
sx q[1];
rz(1.8484263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3446418) q[0];
sx q[0];
rz(-2.0984637) q[0];
sx q[0];
rz(-0.54191858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3082383) q[2];
sx q[2];
rz(-1.6505604) q[2];
sx q[2];
rz(0.74918109) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.433328) q[1];
sx q[1];
rz(-0.86769968) q[1];
sx q[1];
rz(-0.82872699) q[1];
x q[2];
rz(3.1035091) q[3];
sx q[3];
rz(-1.786219) q[3];
sx q[3];
rz(1.076339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30696067) q[2];
sx q[2];
rz(-2.0768276) q[2];
sx q[2];
rz(-1.5131697) q[2];
rz(-0.73897922) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(1.6966604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4479248) q[0];
sx q[0];
rz(-1.6532927) q[0];
sx q[0];
rz(0.36201763) q[0];
rz(1.1297049) q[1];
sx q[1];
rz(-1.874322) q[1];
sx q[1];
rz(2.3806908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7603455) q[0];
sx q[0];
rz(-1.3804304) q[0];
sx q[0];
rz(-2.3668853) q[0];
x q[1];
rz(2.731853) q[2];
sx q[2];
rz(-0.86977661) q[2];
sx q[2];
rz(1.0052731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80428159) q[1];
sx q[1];
rz(-1.7836522) q[1];
sx q[1];
rz(1.0237184) q[1];
rz(-pi) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-2.1511249) q[3];
sx q[3];
rz(2.1023242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1165983) q[2];
sx q[2];
rz(-2.8572539) q[2];
sx q[2];
rz(1.9653758) q[2];
rz(-2.1182012) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(-2.7682176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1171595) q[0];
sx q[0];
rz(-0.42851448) q[0];
sx q[0];
rz(-1.3065216) q[0];
rz(-3.040763) q[1];
sx q[1];
rz(-1.3918624) q[1];
sx q[1];
rz(0.94591013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3358766) q[0];
sx q[0];
rz(-1.44084) q[0];
sx q[0];
rz(-2.8348599) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51351764) q[2];
sx q[2];
rz(-1.1217204) q[2];
sx q[2];
rz(-2.8359536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1506699) q[1];
sx q[1];
rz(-1.7062108) q[1];
sx q[1];
rz(1.5110498) q[1];
rz(-pi) q[2];
rz(-2.2846388) q[3];
sx q[3];
rz(-1.3009239) q[3];
sx q[3];
rz(-1.2601579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1176318) q[2];
sx q[2];
rz(-1.5317711) q[2];
sx q[2];
rz(-2.7839933) q[2];
rz(-3.0910953) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(-2.7549506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9015775) q[0];
sx q[0];
rz(-1.8223263) q[0];
sx q[0];
rz(2.5460119) q[0];
rz(-1.9161842) q[1];
sx q[1];
rz(-1.3662246) q[1];
sx q[1];
rz(2.7092686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29994592) q[0];
sx q[0];
rz(-2.6213745) q[0];
sx q[0];
rz(-0.33824091) q[0];
rz(-0.064266845) q[2];
sx q[2];
rz(-1.7241577) q[2];
sx q[2];
rz(0.72976412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22031517) q[1];
sx q[1];
rz(-2.5743458) q[1];
sx q[1];
rz(-0.0023283255) q[1];
rz(2.6981399) q[3];
sx q[3];
rz(-0.93741527) q[3];
sx q[3];
rz(-1.2284281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7633742) q[2];
sx q[2];
rz(-1.32064) q[2];
sx q[2];
rz(-1.0279559) q[2];
rz(-1.6592735) q[3];
sx q[3];
rz(-1.4488528) q[3];
sx q[3];
rz(-0.75215522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1886275) q[0];
sx q[0];
rz(-1.231622) q[0];
sx q[0];
rz(-2.6787483) q[0];
rz(-1.2340744) q[1];
sx q[1];
rz(-1.9201098) q[1];
sx q[1];
rz(-3*pi/8) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60607227) q[0];
sx q[0];
rz(-1.7458445) q[0];
sx q[0];
rz(-0.025786215) q[0];
rz(-pi) q[1];
rz(-1.0888363) q[2];
sx q[2];
rz(-0.83363998) q[2];
sx q[2];
rz(2.9046975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14069362) q[1];
sx q[1];
rz(-1.447666) q[1];
sx q[1];
rz(0.85832125) q[1];
rz(-pi) q[2];
rz(-0.34171243) q[3];
sx q[3];
rz(-1.6843421) q[3];
sx q[3];
rz(0.36817238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50596333) q[2];
sx q[2];
rz(-1.3869945) q[2];
sx q[2];
rz(-2.3734234) q[2];
rz(0.65028894) q[3];
sx q[3];
rz(-2.8322329) q[3];
sx q[3];
rz(2.5816141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77785093) q[0];
sx q[0];
rz(-2.7010481) q[0];
sx q[0];
rz(-1.995723) q[0];
rz(2.5762985) q[1];
sx q[1];
rz(-2.4259613) q[1];
sx q[1];
rz(-0.40564767) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7956896) q[0];
sx q[0];
rz(-2.3064605) q[0];
sx q[0];
rz(-0.47685868) q[0];
rz(-1.1937386) q[2];
sx q[2];
rz(-1.8296281) q[2];
sx q[2];
rz(-1.4260615) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7024955) q[1];
sx q[1];
rz(-1.8232574) q[1];
sx q[1];
rz(0.52045937) q[1];
x q[2];
rz(0.43934699) q[3];
sx q[3];
rz(-1.6682819) q[3];
sx q[3];
rz(-1.1337224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.046935) q[2];
sx q[2];
rz(-1.9104702) q[2];
sx q[2];
rz(2.6222353) q[2];
rz(0.46118593) q[3];
sx q[3];
rz(-2.0824194) q[3];
sx q[3];
rz(-0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2431348) q[0];
sx q[0];
rz(-1.5903789) q[0];
sx q[0];
rz(-0.90573885) q[0];
rz(1.657919) q[1];
sx q[1];
rz(-1.9762194) q[1];
sx q[1];
rz(-1.2414052) q[1];
rz(1.3686539) q[2];
sx q[2];
rz(-0.62424651) q[2];
sx q[2];
rz(-0.14870283) q[2];
rz(-2.6775581) q[3];
sx q[3];
rz(-1.1190718) q[3];
sx q[3];
rz(-3.0116871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
