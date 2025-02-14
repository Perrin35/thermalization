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
rz(2.4647291) q[0];
sx q[0];
rz(-1.530175) q[0];
sx q[0];
rz(-0.29729182) q[0];
rz(-5.5041432) q[1];
sx q[1];
rz(6.1225523) q[1];
sx q[1];
rz(10.860722) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4477618) q[0];
sx q[0];
rz(-1.0658456) q[0];
sx q[0];
rz(-2.2044066) q[0];
rz(3.0129635) q[2];
sx q[2];
rz(-1.8278215) q[2];
sx q[2];
rz(2.2976053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0432192) q[1];
sx q[1];
rz(-2.3465152) q[1];
sx q[1];
rz(-0.45795346) q[1];
rz(-pi) q[2];
rz(2.8332356) q[3];
sx q[3];
rz(-1.7614109) q[3];
sx q[3];
rz(-2.7558143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5433189) q[2];
sx q[2];
rz(-2.7136549) q[2];
sx q[2];
rz(0.16779009) q[2];
rz(-1.1281891) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096864916) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(2.0665533) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(-1.3137821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15720651) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(2.8048672) q[0];
rz(2.9552835) q[2];
sx q[2];
rz(-1.4042451) q[2];
sx q[2];
rz(-3.0465517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35006501) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(-1.4552874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1978549) q[3];
sx q[3];
rz(-0.35344703) q[3];
sx q[3];
rz(-1.8693876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4904867) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(2.7172078) q[2];
rz(-0.73244798) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(-1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(-0.5027698) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(2.3263993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.316534) q[0];
sx q[0];
rz(-1.382917) q[0];
sx q[0];
rz(1.618439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4645459) q[2];
sx q[2];
rz(-1.4521404) q[2];
sx q[2];
rz(-3.075171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3383662) q[1];
sx q[1];
rz(-2.7684281) q[1];
sx q[1];
rz(0.79014059) q[1];
rz(-pi) q[2];
rz(0.87343862) q[3];
sx q[3];
rz(-2.0202352) q[3];
sx q[3];
rz(-3.0118503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0649071) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(-2.0150851) q[2];
rz(-1.3569776) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(3.0016541) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821238) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(0.35991392) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(2.7684033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554337) q[0];
sx q[0];
rz(-2.9669274) q[0];
sx q[0];
rz(0.96439488) q[0];
rz(-2.7794851) q[2];
sx q[2];
rz(-0.59382861) q[2];
sx q[2];
rz(0.98021843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14228257) q[1];
sx q[1];
rz(-1.9110084) q[1];
sx q[1];
rz(-0.55026023) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3803102) q[3];
sx q[3];
rz(-2.0107186) q[3];
sx q[3];
rz(2.2504077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.084426247) q[2];
sx q[2];
rz(-1.4320222) q[2];
sx q[2];
rz(-0.90844321) q[2];
rz(-2.0716136) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(0.98328868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95803607) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(2.2285484) q[0];
rz(-2.4519582) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(2.3675809) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15439437) q[0];
sx q[0];
rz(-0.78332213) q[0];
sx q[0];
rz(-2.4670967) q[0];
rz(-1.5390603) q[2];
sx q[2];
rz(-0.77546111) q[2];
sx q[2];
rz(-2.8808978) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9076428) q[1];
sx q[1];
rz(-1.7617258) q[1];
sx q[1];
rz(-1.1348866) q[1];
rz(0.048581913) q[3];
sx q[3];
rz(-0.54837275) q[3];
sx q[3];
rz(-0.89357943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.245605) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(-2.4020933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5189811) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(1.125289) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3430905) q[0];
sx q[0];
rz(-2.6406086) q[0];
sx q[0];
rz(0.72879099) q[0];
x q[1];
rz(-2.7281458) q[2];
sx q[2];
rz(-1.3020497) q[2];
sx q[2];
rz(0.15971249) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9003657) q[1];
sx q[1];
rz(-1.3123742) q[1];
sx q[1];
rz(-1.4421806) q[1];
rz(-0.31611021) q[3];
sx q[3];
rz(-2.5694642) q[3];
sx q[3];
rz(-0.24219777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3681727) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(1.0991905) q[2];
rz(2.1558971) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2987591) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(0.66566268) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(-0.98141247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2546651) q[0];
sx q[0];
rz(-1.5136295) q[0];
sx q[0];
rz(-2.1421823) q[0];
rz(-pi) q[1];
rz(1.8265884) q[2];
sx q[2];
rz(-1.0793387) q[2];
sx q[2];
rz(-1.0813528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8725207) q[1];
sx q[1];
rz(-1.4351294) q[1];
sx q[1];
rz(-1.1024464) q[1];
x q[2];
rz(2.2432765) q[3];
sx q[3];
rz(-1.6784462) q[3];
sx q[3];
rz(-0.09569351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0770646) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(-0.36901739) q[2];
rz(1.7391694) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(1.3212475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3082665) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(0.33018026) q[0];
rz(-0.45686832) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(-0.59828573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5766524) q[0];
sx q[0];
rz(-2.1974753) q[0];
sx q[0];
rz(-2.1184854) q[0];
rz(0.070934709) q[2];
sx q[2];
rz(-1.1026376) q[2];
sx q[2];
rz(0.93229655) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6574994) q[1];
sx q[1];
rz(-1.7998905) q[1];
sx q[1];
rz(-2.799396) q[1];
x q[2];
rz(-2.4870858) q[3];
sx q[3];
rz(-1.084436) q[3];
sx q[3];
rz(2.8974722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8699708) q[2];
sx q[2];
rz(-0.87646708) q[2];
sx q[2];
rz(-0.064621933) q[2];
rz(-1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86218631) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(-2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(-1.0313787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169075) q[0];
sx q[0];
rz(-2.3992736) q[0];
sx q[0];
rz(-0.10984212) q[0];
x q[1];
rz(-1.4508444) q[2];
sx q[2];
rz(-1.2799096) q[2];
sx q[2];
rz(-1.5121258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045193521) q[1];
sx q[1];
rz(-0.61122433) q[1];
sx q[1];
rz(-1.0525714) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1886173) q[3];
sx q[3];
rz(-2.6572356) q[3];
sx q[3];
rz(-0.25121197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3970268) q[2];
sx q[2];
rz(-2.4137745) q[2];
sx q[2];
rz(-0.087372027) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(-0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(-2.0024306) q[0];
rz(2.6423404) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(-1.8081236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3402417) q[0];
sx q[0];
rz(-1.6610067) q[0];
sx q[0];
rz(-3.0349681) q[0];
rz(1.75976) q[2];
sx q[2];
rz(-0.25936959) q[2];
sx q[2];
rz(-1.237902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7391679) q[1];
sx q[1];
rz(-0.3350733) q[1];
sx q[1];
rz(0.70783028) q[1];
rz(-pi) q[2];
rz(0.74204294) q[3];
sx q[3];
rz(-0.84845209) q[3];
sx q[3];
rz(-2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.8583543) q[2];
rz(-1.2011652) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82861154) q[0];
sx q[0];
rz(-0.77684488) q[0];
sx q[0];
rz(-2.0364398) q[0];
rz(-3.037187) q[1];
sx q[1];
rz(-1.5900292) q[1];
sx q[1];
rz(1.1463696) q[1];
rz(0.5055867) q[2];
sx q[2];
rz(-1.3973631) q[2];
sx q[2];
rz(0.33968795) q[2];
rz(2.6324568) q[3];
sx q[3];
rz(-2.4066299) q[3];
sx q[3];
rz(0.015711333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
