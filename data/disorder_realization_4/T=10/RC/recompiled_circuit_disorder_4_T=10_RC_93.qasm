OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.714158) q[0];
sx q[0];
rz(-2.7058869) q[0];
sx q[0];
rz(-0.92619196) q[0];
rz(-1.1821094) q[1];
sx q[1];
rz(3.8745772) q[1];
sx q[1];
rz(12.193845) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6719138) q[0];
sx q[0];
rz(-0.9979453) q[0];
sx q[0];
rz(-1.8125305) q[0];
rz(-pi) q[1];
rz(-1.6765321) q[2];
sx q[2];
rz(-2.1991962) q[2];
sx q[2];
rz(0.65199967) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8358742) q[1];
sx q[1];
rz(-0.22138518) q[1];
sx q[1];
rz(-2.1165044) q[1];
rz(1.3543234) q[3];
sx q[3];
rz(-1.4604005) q[3];
sx q[3];
rz(0.9339827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42840502) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(2.2170846) q[2];
rz(1.6690147) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532042) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(2.9470782) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40614906) q[0];
sx q[0];
rz(-0.14649728) q[0];
sx q[0];
rz(2.2551401) q[0];
rz(-0.61972159) q[2];
sx q[2];
rz(-1.1148858) q[2];
sx q[2];
rz(-1.0085307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3000852) q[1];
sx q[1];
rz(-1.2282787) q[1];
sx q[1];
rz(3.1373346) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8781434) q[3];
sx q[3];
rz(-0.69034319) q[3];
sx q[3];
rz(-3.1171947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7028971) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(-0.99059087) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.3593339) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(-0.24770501) q[1];
sx q[1];
rz(-1.1876371) q[1];
sx q[1];
rz(1.8050271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7860387) q[0];
sx q[0];
rz(-0.22909129) q[0];
sx q[0];
rz(-1.3130207) q[0];
rz(0.58972085) q[2];
sx q[2];
rz(-1.9460742) q[2];
sx q[2];
rz(-0.42082618) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6429236) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(0.81865262) q[1];
x q[2];
rz(0.41575723) q[3];
sx q[3];
rz(-0.81167479) q[3];
sx q[3];
rz(0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(2.1170199) q[2];
rz(-1.2767977) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(-1.4039325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0899635) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(-0.35332638) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68217512) q[2];
sx q[2];
rz(-2.1795159) q[2];
sx q[2];
rz(-1.9145554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6744068) q[1];
sx q[1];
rz(-1.8605839) q[1];
sx q[1];
rz(-3.0249216) q[1];
x q[2];
rz(-1.4675667) q[3];
sx q[3];
rz(-0.69403115) q[3];
sx q[3];
rz(-0.3895143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(-2.7385353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-2.1257341) q[0];
rz(-1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(-3.1075081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344454) q[0];
sx q[0];
rz(-2.7140518) q[0];
sx q[0];
rz(-2.7464675) q[0];
rz(0.39880619) q[2];
sx q[2];
rz(-0.6119234) q[2];
sx q[2];
rz(-1.2955701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9220229) q[1];
sx q[1];
rz(-1.8551738) q[1];
sx q[1];
rz(-0.0071830458) q[1];
rz(-2.2250416) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(-1.786701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74026996) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(1.3327538) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.4532582) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(0.43584287) q[0];
rz(-1.1054989) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(1.9357392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57731956) q[0];
sx q[0];
rz(-0.84007971) q[0];
sx q[0];
rz(0.3198091) q[0];
rz(-1.4044936) q[2];
sx q[2];
rz(-3.0450833) q[2];
sx q[2];
rz(-1.4788747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44991325) q[1];
sx q[1];
rz(-2.0813585) q[1];
sx q[1];
rz(-1.6010067) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1084656) q[3];
sx q[3];
rz(-1.1713542) q[3];
sx q[3];
rz(1.2101733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(0.40766019) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(-2.6962962) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62149858) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(-2.4682585) q[0];
rz(0.78397059) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(2.6046682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250701) q[0];
sx q[0];
rz(-1.636583) q[0];
sx q[0];
rz(1.487088) q[0];
rz(1.827048) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(-1.8585376) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1310136) q[1];
sx q[1];
rz(-0.77242935) q[1];
sx q[1];
rz(-1.7060682) q[1];
rz(-pi) q[2];
rz(1.4173996) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(-0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15988222) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(-0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(0.38761815) q[0];
rz(0.11501137) q[1];
sx q[1];
rz(-1.7128046) q[1];
sx q[1];
rz(-2.8040335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885376) q[0];
sx q[0];
rz(-2.6153784) q[0];
sx q[0];
rz(-0.40286581) q[0];
rz(-1.6667716) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(-2.1542187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90606373) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(-0.99793418) q[1];
rz(-pi) q[2];
rz(0.75546219) q[3];
sx q[3];
rz(-1.5243693) q[3];
sx q[3];
rz(-1.4481627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(-1.8743275) q[2];
rz(-0.77504843) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(-0.22098456) q[0];
rz(2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-1.0029213) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5054277) q[0];
sx q[0];
rz(-1.3190077) q[0];
sx q[0];
rz(-1.4279143) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9858176) q[2];
sx q[2];
rz(-2.0009811) q[2];
sx q[2];
rz(2.1625724) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3993966) q[1];
sx q[1];
rz(-2.4376215) q[1];
sx q[1];
rz(-0.75615551) q[1];
x q[2];
rz(2.8691611) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(1.2215134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(-0.15850244) q[2];
rz(-0.45378271) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52699387) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(0.19113834) q[0];
rz(-0.29516164) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(0.89231649) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96517262) q[0];
sx q[0];
rz(-1.7473548) q[0];
sx q[0];
rz(0.82133349) q[0];
x q[1];
rz(3.0478165) q[2];
sx q[2];
rz(-0.9320335) q[2];
sx q[2];
rz(0.057657777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.416623) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(-0.23857393) q[1];
x q[2];
rz(0.51384135) q[3];
sx q[3];
rz(-2.6548879) q[3];
sx q[3];
rz(2.0972308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3698547) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-2.2820293) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(-1.6620811) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(-3.1142278) q[2];
sx q[2];
rz(-1.8623427) q[2];
sx q[2];
rz(1.3775415) q[2];
rz(2.0966093) q[3];
sx q[3];
rz(-3.0621739) q[3];
sx q[3];
rz(-0.0041181507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];