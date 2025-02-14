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
rz(-0.8987838) q[0];
sx q[0];
rz(-1.838983) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(1.6279434) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0078137579) q[0];
sx q[0];
rz(-1.8651175) q[0];
sx q[0];
rz(-1.9208287) q[0];
rz(-pi) q[1];
rz(0.75998016) q[2];
sx q[2];
rz(-2.438025) q[2];
sx q[2];
rz(0.68389308) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.096949654) q[1];
sx q[1];
rz(-1.5162807) q[1];
sx q[1];
rz(1.9115555) q[1];
rz(-2.1964425) q[3];
sx q[3];
rz(-1.8813475) q[3];
sx q[3];
rz(-2.9356082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0579494) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(2.7570214) q[2];
rz(0.40258506) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3835417) q[0];
sx q[0];
rz(-1.6076247) q[0];
sx q[0];
rz(0.59355271) q[0];
rz(-1.3213762) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(-0.039332565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7205374) q[0];
sx q[0];
rz(-2.5364784) q[0];
sx q[0];
rz(2.2626816) q[0];
rz(-pi) q[1];
rz(-1.4565384) q[2];
sx q[2];
rz(-2.2041177) q[2];
sx q[2];
rz(0.81482023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96728863) q[1];
sx q[1];
rz(-1.6203462) q[1];
sx q[1];
rz(0.51731717) q[1];
rz(0.049841471) q[3];
sx q[3];
rz(-2.5916272) q[3];
sx q[3];
rz(2.0932239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3010657) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(-0.62704101) q[2];
rz(0.4380694) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(-2.0048678) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7196734) q[0];
sx q[0];
rz(-2.0543583) q[0];
sx q[0];
rz(1.3713974) q[0];
rz(-0.60945359) q[1];
sx q[1];
rz(-2.2513159) q[1];
sx q[1];
rz(-1.2249464) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44642513) q[0];
sx q[0];
rz(-1.37943) q[0];
sx q[0];
rz(1.4728949) q[0];
rz(-2.0535499) q[2];
sx q[2];
rz(-2.4039589) q[2];
sx q[2];
rz(0.18407735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65388262) q[1];
sx q[1];
rz(-1.9284029) q[1];
sx q[1];
rz(-0.34642152) q[1];
x q[2];
rz(-1.8277934) q[3];
sx q[3];
rz(-0.71201268) q[3];
sx q[3];
rz(0.76699996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3043392) q[2];
sx q[2];
rz(-0.86696583) q[2];
sx q[2];
rz(-0.85401094) q[2];
rz(-0.52418661) q[3];
sx q[3];
rz(-1.6322522) q[3];
sx q[3];
rz(-0.9700276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22685856) q[0];
sx q[0];
rz(-0.15758841) q[0];
sx q[0];
rz(2.3478813) q[0];
rz(-0.0018250068) q[1];
sx q[1];
rz(-2.5853214) q[1];
sx q[1];
rz(0.8108286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13769606) q[0];
sx q[0];
rz(-0.76327774) q[0];
sx q[0];
rz(1.3062817) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89703441) q[2];
sx q[2];
rz(-1.9193726) q[2];
sx q[2];
rz(-0.86212117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6135679) q[1];
sx q[1];
rz(-2.151229) q[1];
sx q[1];
rz(2.7356262) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9998897) q[3];
sx q[3];
rz(-1.4388814) q[3];
sx q[3];
rz(-2.9219081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4578555) q[2];
sx q[2];
rz(-2.808414) q[2];
sx q[2];
rz(1.4972756) q[2];
rz(1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(1.3076521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.4182279) q[0];
sx q[0];
rz(-1.6809373) q[0];
sx q[0];
rz(-0.25890589) q[0];
rz(0.12340165) q[1];
sx q[1];
rz(-2.5105208) q[1];
sx q[1];
rz(-2.6554328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3116592) q[0];
sx q[0];
rz(-1.0898542) q[0];
sx q[0];
rz(0.89998683) q[0];
x q[1];
rz(-1.4347378) q[2];
sx q[2];
rz(-1.1216838) q[2];
sx q[2];
rz(2.8633008) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9100719) q[1];
sx q[1];
rz(-1.2086444) q[1];
sx q[1];
rz(-2.9400184) q[1];
rz(0.87032373) q[3];
sx q[3];
rz(-2.016267) q[3];
sx q[3];
rz(-1.1401759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90998021) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(-0.95476556) q[3];
sx q[3];
rz(-1.747811) q[3];
sx q[3];
rz(2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93916494) q[0];
sx q[0];
rz(-1.9068149) q[0];
sx q[0];
rz(0.760461) q[0];
rz(-0.83433926) q[1];
sx q[1];
rz(-2.0400679) q[1];
sx q[1];
rz(2.0371425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6532947) q[0];
sx q[0];
rz(-1.0279995) q[0];
sx q[0];
rz(-0.25570583) q[0];
rz(-0.8669187) q[2];
sx q[2];
rz(-0.50627497) q[2];
sx q[2];
rz(2.3828854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1431434) q[1];
sx q[1];
rz(-2.1996212) q[1];
sx q[1];
rz(2.5004908) q[1];
rz(-0.99640019) q[3];
sx q[3];
rz(-1.8339012) q[3];
sx q[3];
rz(-0.97364473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1793648) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(0.041570138) q[2];
rz(-2.9465594) q[3];
sx q[3];
rz(-0.2427559) q[3];
sx q[3];
rz(-0.78708762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70260173) q[0];
sx q[0];
rz(-0.85346237) q[0];
sx q[0];
rz(0.75337291) q[0];
rz(2.0808749) q[1];
sx q[1];
rz(-0.80442387) q[1];
sx q[1];
rz(0.20739584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30880659) q[0];
sx q[0];
rz(-2.1784003) q[0];
sx q[0];
rz(0.73143801) q[0];
rz(-1.2815525) q[2];
sx q[2];
rz(-1.1221894) q[2];
sx q[2];
rz(1.9236652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.064697178) q[1];
sx q[1];
rz(-1.3363779) q[1];
sx q[1];
rz(-0.56273048) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5785406) q[3];
sx q[3];
rz(-1.6098579) q[3];
sx q[3];
rz(-2.9970053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14219001) q[2];
sx q[2];
rz(-1.9839857) q[2];
sx q[2];
rz(3.0863777) q[2];
rz(2.2986872) q[3];
sx q[3];
rz(-2.3838145) q[3];
sx q[3];
rz(-0.88111669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371599) q[0];
sx q[0];
rz(-0.70796767) q[0];
sx q[0];
rz(1.1771033) q[0];
rz(-2.2616995) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(0.060861977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5755479) q[0];
sx q[0];
rz(-1.6948997) q[0];
sx q[0];
rz(-0.095416165) q[0];
x q[1];
rz(1.8905156) q[2];
sx q[2];
rz(-2.0315315) q[2];
sx q[2];
rz(1.2416897) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7877373) q[1];
sx q[1];
rz(-1.0050049) q[1];
sx q[1];
rz(-1.4325527) q[1];
rz(-1.8350321) q[3];
sx q[3];
rz(-1.6420393) q[3];
sx q[3];
rz(2.3087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.047478288) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1020553) q[0];
sx q[0];
rz(-0.85298959) q[0];
sx q[0];
rz(0.0041740388) q[0];
rz(-1.8904842) q[1];
sx q[1];
rz(-3.0079542) q[1];
sx q[1];
rz(-2.4600696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96723667) q[0];
sx q[0];
rz(-2.3612494) q[0];
sx q[0];
rz(-1.1514949) q[0];
rz(-pi) q[1];
rz(-0.84284346) q[2];
sx q[2];
rz(-2.3327391) q[2];
sx q[2];
rz(0.12241546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8261823) q[1];
sx q[1];
rz(-1.8528954) q[1];
sx q[1];
rz(-1.3768857) q[1];
x q[2];
rz(0.32855804) q[3];
sx q[3];
rz(-0.78379455) q[3];
sx q[3];
rz(-0.72600049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5847797) q[2];
sx q[2];
rz(-0.4759554) q[2];
sx q[2];
rz(1.7468096) q[2];
rz(-0.41633385) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(-0.39286119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3861179) q[0];
sx q[0];
rz(-0.32805726) q[0];
sx q[0];
rz(-2.1395444) q[0];
rz(-0.26456061) q[1];
sx q[1];
rz(-1.325565) q[1];
sx q[1];
rz(1.4904259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39367657) q[0];
sx q[0];
rz(-1.0320469) q[0];
sx q[0];
rz(1.7231621) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4117496) q[2];
sx q[2];
rz(-2.3727086) q[2];
sx q[2];
rz(-1.8184219) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24764168) q[1];
sx q[1];
rz(-0.14130768) q[1];
sx q[1];
rz(2.2013478) q[1];
rz(-pi) q[2];
rz(-1.6308925) q[3];
sx q[3];
rz(-1.6625685) q[3];
sx q[3];
rz(-2.0962146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1313021) q[2];
sx q[2];
rz(-1.7346069) q[2];
sx q[2];
rz(1.4170125) q[2];
rz(2.8085282) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(1.7634348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604816) q[0];
sx q[0];
rz(-1.2662553) q[0];
sx q[0];
rz(1.8929831) q[0];
rz(-0.026451182) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(2.6698723) q[2];
sx q[2];
rz(-0.77009554) q[2];
sx q[2];
rz(-0.64252616) q[2];
rz(2.706017) q[3];
sx q[3];
rz(-0.66146225) q[3];
sx q[3];
rz(0.41730455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
