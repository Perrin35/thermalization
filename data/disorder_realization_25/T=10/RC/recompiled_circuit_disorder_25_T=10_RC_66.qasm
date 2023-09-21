OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1782921) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(0.51142366) q[0];
rz(-2.7152039) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(-0.79977712) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3275798) q[1];
sx q[1];
rz(-0.50427932) q[1];
sx q[1];
rz(-3.1370647) q[1];
x q[2];
rz(-1.5929234) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(-2.4165137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(2.5722356) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5434718) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(0.55364451) q[0];
rz(1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935826) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(2.4539024) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4670062) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(-1.9763725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1439647) q[1];
sx q[1];
rz(-1.3647172) q[1];
sx q[1];
rz(-3.1254083) q[1];
x q[2];
rz(-1.3319098) q[3];
sx q[3];
rz(-2.2242862) q[3];
sx q[3];
rz(-1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(-2.3705204) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-3.070389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6547346) q[0];
sx q[0];
rz(-1.808597) q[0];
sx q[0];
rz(-2.9391857) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0729012) q[2];
sx q[2];
rz(-2.6325668) q[2];
sx q[2];
rz(-0.44391649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5350266) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(-1.4024629) q[1];
rz(-0.86665385) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(1.9748953) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3281876) q[0];
sx q[0];
rz(-1.7338599) q[0];
sx q[0];
rz(-1.3379315) q[0];
rz(-pi) q[1];
rz(-0.0091323098) q[2];
sx q[2];
rz(-1.3924686) q[2];
sx q[2];
rz(0.4724617) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64926681) q[1];
sx q[1];
rz(-0.71486231) q[1];
sx q[1];
rz(1.8148755) q[1];
rz(-pi) q[2];
rz(2.5656384) q[3];
sx q[3];
rz(-2.7154657) q[3];
sx q[3];
rz(1.061561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.2129983) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0096668) q[0];
sx q[0];
rz(-1.8312694) q[0];
sx q[0];
rz(2.7564604) q[0];
rz(2.2239387) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(1.7078924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4124406) q[1];
sx q[1];
rz(-1.8262595) q[1];
sx q[1];
rz(1.6925807) q[1];
rz(-pi) q[2];
rz(-0.17739399) q[3];
sx q[3];
rz(-1.8284441) q[3];
sx q[3];
rz(-0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-0.58475959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8557381) q[0];
sx q[0];
rz(-1.8830839) q[0];
sx q[0];
rz(-0.62494846) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61442394) q[2];
sx q[2];
rz(-1.5922058) q[2];
sx q[2];
rz(1.1504088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11440052) q[1];
sx q[1];
rz(-2.235734) q[1];
sx q[1];
rz(2.7726638) q[1];
rz(2.8486409) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(-2.0041182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.6284846) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-0.54876304) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(-2.5351977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8618146) q[0];
sx q[0];
rz(-1.8068411) q[0];
sx q[0];
rz(1.0311014) q[0];
rz(-pi) q[1];
rz(0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(0.19336685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.079271) q[1];
sx q[1];
rz(-2.3289526) q[1];
sx q[1];
rz(-0.15516849) q[1];
rz(2.2783979) q[3];
sx q[3];
rz(-0.51897012) q[3];
sx q[3];
rz(0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(0.50841224) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4906625) q[0];
sx q[0];
rz(-2.7224446) q[0];
sx q[0];
rz(-1.3952257) q[0];
rz(-1.9367427) q[2];
sx q[2];
rz(-1.9908675) q[2];
sx q[2];
rz(0.36047381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7714872) q[1];
sx q[1];
rz(-2.5574554) q[1];
sx q[1];
rz(0.41569709) q[1];
x q[2];
rz(-1.6119192) q[3];
sx q[3];
rz(-1.0121945) q[3];
sx q[3];
rz(-2.5839644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.078538744) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(-2.3030346) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(1.1484336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6156857) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(0.023948897) q[0];
rz(-pi) q[1];
rz(-2.1512595) q[2];
sx q[2];
rz(-0.60539421) q[2];
sx q[2];
rz(2.6661172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.941274) q[1];
sx q[1];
rz(-2.4793844) q[1];
sx q[1];
rz(2.3108285) q[1];
rz(-0.095454772) q[3];
sx q[3];
rz(-2.7881552) q[3];
sx q[3];
rz(-0.27775882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(2.443312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822599) q[0];
sx q[0];
rz(-1.2685304) q[0];
sx q[0];
rz(-1.4955826) q[0];
rz(0.024784879) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(0.63400808) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2472898) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(1.1029907) q[1];
rz(-0.24153696) q[3];
sx q[3];
rz(-2.1250484) q[3];
sx q[3];
rz(2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(0.60224709) q[2];
sx q[2];
rz(-1.2802274) q[2];
sx q[2];
rz(1.0614492) q[2];
rz(2.2610353) q[3];
sx q[3];
rz(-1.4546483) q[3];
sx q[3];
rz(1.6402257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
