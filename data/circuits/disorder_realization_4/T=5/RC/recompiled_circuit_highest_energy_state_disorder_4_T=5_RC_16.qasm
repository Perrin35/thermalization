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
rz(0.29851222) q[0];
sx q[0];
rz(4.4311509) q[0];
sx q[0];
rz(9.4480954) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(-1.4758543) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97853547) q[0];
sx q[0];
rz(-3.0594303) q[0];
sx q[0];
rz(-1.0121679) q[0];
rz(-pi) q[1];
rz(1.453773) q[2];
sx q[2];
rz(-0.3323148) q[2];
sx q[2];
rz(1.8154643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45881328) q[1];
sx q[1];
rz(-0.87535697) q[1];
sx q[1];
rz(-2.0942445) q[1];
rz(-0.9528927) q[3];
sx q[3];
rz(-2.3639819) q[3];
sx q[3];
rz(-2.4617406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4618571) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(0.30065817) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41475007) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-2.9657189) q[0];
rz(2.580592) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(-0.19634136) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4787653) q[0];
sx q[0];
rz(-1.2365554) q[0];
sx q[0];
rz(-1.1007376) q[0];
rz(-pi) q[1];
rz(2.3679737) q[2];
sx q[2];
rz(-0.16628312) q[2];
sx q[2];
rz(2.6321049) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93028101) q[1];
sx q[1];
rz(-2.086211) q[1];
sx q[1];
rz(-3.1392908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7815468) q[3];
sx q[3];
rz(-2.776675) q[3];
sx q[3];
rz(-1.5860189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.089513) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(-1.4698131) q[2];
rz(2.6164264) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(2.4770881) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0199652) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(2.5110733) q[0];
rz(-1.9969253) q[1];
sx q[1];
rz(-1.2455218) q[1];
sx q[1];
rz(0.05680457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5543514) q[0];
sx q[0];
rz(-1.9375174) q[0];
sx q[0];
rz(-1.1318867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0089097) q[2];
sx q[2];
rz(-2.3724634) q[2];
sx q[2];
rz(3.1388856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2395798) q[1];
sx q[1];
rz(-2.6551592) q[1];
sx q[1];
rz(0.3063267) q[1];
rz(1.1435236) q[3];
sx q[3];
rz(-1.3766854) q[3];
sx q[3];
rz(0.33188348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8593665) q[2];
sx q[2];
rz(-1.6702009) q[2];
sx q[2];
rz(-2.7024506) q[2];
rz(1.3965083) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948931) q[0];
sx q[0];
rz(-0.72143227) q[0];
sx q[0];
rz(-1.0149581) q[0];
rz(-3.0789442) q[1];
sx q[1];
rz(-2.1868314) q[1];
sx q[1];
rz(0.28422022) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260176) q[0];
sx q[0];
rz(-1.9997445) q[0];
sx q[0];
rz(1.9563849) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8880803) q[2];
sx q[2];
rz(-0.51157839) q[2];
sx q[2];
rz(1.6588841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1831467) q[1];
sx q[1];
rz(-1.5827521) q[1];
sx q[1];
rz(-0.35838106) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6246068) q[3];
sx q[3];
rz(-2.887886) q[3];
sx q[3];
rz(-2.1068609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8404428) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-1.1287929) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(0.40529761) q[0];
rz(0.48794508) q[1];
sx q[1];
rz(-2.0351724) q[1];
sx q[1];
rz(-0.36283666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3317446) q[0];
sx q[0];
rz(-2.6001626) q[0];
sx q[0];
rz(0.62744139) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58541457) q[2];
sx q[2];
rz(-1.8903132) q[2];
sx q[2];
rz(-1.8210757) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0367485) q[1];
sx q[1];
rz(-1.7916153) q[1];
sx q[1];
rz(1.5209496) q[1];
rz(-pi) q[2];
rz(0.69546206) q[3];
sx q[3];
rz(-1.3331279) q[3];
sx q[3];
rz(-2.9836751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(-2.4190767) q[2];
rz(-0.51269382) q[3];
sx q[3];
rz(-0.86944681) q[3];
sx q[3];
rz(-1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(-2.784506) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(-1.9836609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5362893) q[0];
sx q[0];
rz(-0.9035631) q[0];
sx q[0];
rz(0.31179223) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91041301) q[2];
sx q[2];
rz(-2.0636301) q[2];
sx q[2];
rz(-0.016066859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.14640644) q[1];
sx q[1];
rz(-2.1651092) q[1];
sx q[1];
rz(-1.8176778) q[1];
rz(-pi) q[2];
rz(2.255106) q[3];
sx q[3];
rz(-0.89596701) q[3];
sx q[3];
rz(-1.0533028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.57924119) q[2];
sx q[2];
rz(-1.3902105) q[2];
sx q[2];
rz(-0.58970279) q[2];
rz(-1.3127182) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(0.5635128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967459) q[0];
sx q[0];
rz(-0.99269301) q[0];
sx q[0];
rz(-2.8712811) q[0];
rz(0.42701328) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(-0.40831533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664551) q[0];
sx q[0];
rz(-1.9388559) q[0];
sx q[0];
rz(1.6987263) q[0];
x q[1];
rz(0.10119844) q[2];
sx q[2];
rz(-2.2658341) q[2];
sx q[2];
rz(-1.1388495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.22173126) q[1];
sx q[1];
rz(-1.1581005) q[1];
sx q[1];
rz(0.92642118) q[1];
rz(-pi) q[2];
rz(0.72260489) q[3];
sx q[3];
rz(-1.6480903) q[3];
sx q[3];
rz(-0.80394905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.938505) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(-0.92424029) q[2];
rz(0.83485323) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95558178) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(-2.7741449) q[0];
rz(-2.5732749) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(-2.7265991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5552705) q[0];
sx q[0];
rz(-1.30997) q[0];
sx q[0];
rz(1.1010948) q[0];
rz(-pi) q[1];
rz(-2.3363386) q[2];
sx q[2];
rz(-2.4404011) q[2];
sx q[2];
rz(-2.5561577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4492682) q[1];
sx q[1];
rz(-2.0776761) q[1];
sx q[1];
rz(-0.70517069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5852195) q[3];
sx q[3];
rz(-2.7957749) q[3];
sx q[3];
rz(2.7399969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(2.3084194) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(-1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9678765) q[0];
sx q[0];
rz(-0.40206566) q[0];
sx q[0];
rz(-0.014884431) q[0];
rz(2.0170276) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(3.0344149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10725645) q[0];
sx q[0];
rz(-1.7198623) q[0];
sx q[0];
rz(2.3019522) q[0];
rz(-pi) q[1];
rz(0.82294269) q[2];
sx q[2];
rz(-1.2185893) q[2];
sx q[2];
rz(-0.65300452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6139504) q[1];
sx q[1];
rz(-2.6384934) q[1];
sx q[1];
rz(-0.31219074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8802683) q[3];
sx q[3];
rz(-1.7423986) q[3];
sx q[3];
rz(-1.8971407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5772446) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(-2.2819819) q[2];
rz(0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49941007) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(1.1580178) q[0];
rz(-2.6462789) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(-0.29702979) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.934108) q[0];
sx q[0];
rz(-1.5667324) q[0];
sx q[0];
rz(-1.8898176) q[0];
rz(-pi) q[1];
rz(-3.0115602) q[2];
sx q[2];
rz(-2.8388925) q[2];
sx q[2];
rz(-1.5720194) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7630939) q[1];
sx q[1];
rz(-1.8619672) q[1];
sx q[1];
rz(-1.4062455) q[1];
rz(-pi) q[2];
x q[2];
rz(0.016344109) q[3];
sx q[3];
rz(-0.90829231) q[3];
sx q[3];
rz(1.3682883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.872252) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(2.2807109) q[2];
rz(0.46368972) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(-2.004682) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1562445) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.8078177) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(-2.4598224) q[2];
sx q[2];
rz(-2.3538156) q[2];
sx q[2];
rz(0.71629477) q[2];
rz(-2.648223) q[3];
sx q[3];
rz(-1.2262288) q[3];
sx q[3];
rz(-0.83362383) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
