OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(8.0170595) q[0];
sx q[0];
rz(9.4693139) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(-1.2244616) q[1];
sx q[1];
rz(-0.78723025) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48101857) q[0];
sx q[0];
rz(-0.74587599) q[0];
sx q[0];
rz(-2.3700299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2455315) q[2];
sx q[2];
rz(-1.6084387) q[2];
sx q[2];
rz(-1.4311439) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41809717) q[1];
sx q[1];
rz(-0.58291328) q[1];
sx q[1];
rz(-1.6969796) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9028893) q[3];
sx q[3];
rz(-1.9145962) q[3];
sx q[3];
rz(-1.7834922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74701509) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(-0.75418312) q[2];
rz(-1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(-2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542145) q[0];
sx q[0];
rz(-0.35254756) q[0];
sx q[0];
rz(-1.3341599) q[0];
rz(-1.2913903) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(-0.87563595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9459045) q[0];
sx q[0];
rz(-0.92600805) q[0];
sx q[0];
rz(-2.2553205) q[0];
rz(-2.0673758) q[2];
sx q[2];
rz(-1.6255709) q[2];
sx q[2];
rz(-0.77067256) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8671682) q[1];
sx q[1];
rz(-0.53571415) q[1];
sx q[1];
rz(-0.97656004) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3631212) q[3];
sx q[3];
rz(-1.8559578) q[3];
sx q[3];
rz(-0.36237291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.69677532) q[2];
sx q[2];
rz(-0.70187086) q[2];
sx q[2];
rz(-0.87265054) q[2];
rz(2.3526092) q[3];
sx q[3];
rz(-2.240447) q[3];
sx q[3];
rz(0.22855973) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92662421) q[0];
sx q[0];
rz(-1.3525532) q[0];
sx q[0];
rz(-0.0080000814) q[0];
rz(-0.81832111) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(1.2368894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87419009) q[0];
sx q[0];
rz(-1.2498901) q[0];
sx q[0];
rz(-0.85808922) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34864595) q[2];
sx q[2];
rz(-1.1599891) q[2];
sx q[2];
rz(0.88627671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9540959) q[1];
sx q[1];
rz(-1.264123) q[1];
sx q[1];
rz(2.5727954) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.453492) q[3];
sx q[3];
rz(-1.8059732) q[3];
sx q[3];
rz(1.4405574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41039738) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(0.94089874) q[2];
rz(-2.0375552) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(1.1436536) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10821548) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-2.9040842) q[0];
rz(-0.6595276) q[1];
sx q[1];
rz(-2.567465) q[1];
sx q[1];
rz(0.017017078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68786808) q[0];
sx q[0];
rz(-2.9659418) q[0];
sx q[0];
rz(1.8499299) q[0];
rz(-pi) q[1];
rz(-0.41827664) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(1.0470225) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.58910546) q[1];
sx q[1];
rz(-2.1334339) q[1];
sx q[1];
rz(-1.6943114) q[1];
x q[2];
rz(2.3364725) q[3];
sx q[3];
rz(-2.1305354) q[3];
sx q[3];
rz(2.2162958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.235405) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(0.91442937) q[2];
rz(-0.36214456) q[3];
sx q[3];
rz(-0.80949628) q[3];
sx q[3];
rz(-2.5563498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.888716) q[0];
sx q[0];
rz(-1.3739561) q[0];
sx q[0];
rz(-0.1121029) q[0];
rz(1.0653227) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(2.0723453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4010942) q[0];
sx q[0];
rz(-1.2477269) q[0];
sx q[0];
rz(-1.857398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.094355) q[2];
sx q[2];
rz(-2.0531056) q[2];
sx q[2];
rz(0.87076887) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4878975) q[1];
sx q[1];
rz(-1.2733165) q[1];
sx q[1];
rz(-1.9031699) q[1];
rz(0.97158708) q[3];
sx q[3];
rz(-1.4607444) q[3];
sx q[3];
rz(0.19796619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7142882) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(0.53421268) q[2];
rz(0.30321768) q[3];
sx q[3];
rz(-1.5847619) q[3];
sx q[3];
rz(1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7233906) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(2.2499625) q[0];
rz(-1.760969) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(2.2183529) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6398481) q[0];
sx q[0];
rz(-0.11127936) q[0];
sx q[0];
rz(-0.48868816) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60439749) q[2];
sx q[2];
rz(-2.6816419) q[2];
sx q[2];
rz(-2.9876171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90369836) q[1];
sx q[1];
rz(-2.782202) q[1];
sx q[1];
rz(0.87025799) q[1];
rz(-pi) q[2];
rz(-1.1679064) q[3];
sx q[3];
rz(-2.1729762) q[3];
sx q[3];
rz(-0.17113479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0455857) q[2];
sx q[2];
rz(-2.7092689) q[2];
sx q[2];
rz(1.879479) q[2];
rz(-1.9803842) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(-1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68814174) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(-1.5334817) q[0];
rz(2.7104132) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(-1.7327488) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20074318) q[0];
sx q[0];
rz(-2.4398167) q[0];
sx q[0];
rz(-2.6883342) q[0];
rz(1.7869912) q[2];
sx q[2];
rz(-0.39297418) q[2];
sx q[2];
rz(-0.68902389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6366439) q[1];
sx q[1];
rz(-1.726999) q[1];
sx q[1];
rz(-1.6976995) q[1];
rz(-pi) q[2];
x q[2];
rz(0.124229) q[3];
sx q[3];
rz(-2.0696739) q[3];
sx q[3];
rz(1.4829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2048753) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(1.3196866) q[2];
rz(-3.1070869) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.3694793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40783229) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.2169417) q[0];
rz(2.2159684) q[1];
sx q[1];
rz(-1.7763014) q[1];
sx q[1];
rz(1.2006753) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9618098) q[0];
sx q[0];
rz(-2.2968074) q[0];
sx q[0];
rz(2.2711193) q[0];
rz(-pi) q[1];
rz(-0.56270069) q[2];
sx q[2];
rz(-0.68018736) q[2];
sx q[2];
rz(-2.0923751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0266704) q[1];
sx q[1];
rz(-1.611913) q[1];
sx q[1];
rz(2.369472) q[1];
rz(2.2766477) q[3];
sx q[3];
rz(-2.6078913) q[3];
sx q[3];
rz(0.362901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.919148) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-1.0486802) q[2];
rz(-2.9564296) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(-0.10704253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.01934) q[0];
sx q[0];
rz(-2.4801319) q[0];
sx q[0];
rz(-2.3308603) q[0];
rz(2.4108389) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(2.1810541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385242) q[0];
sx q[0];
rz(-1.5674066) q[0];
sx q[0];
rz(1.5698213) q[0];
rz(-0.87709092) q[2];
sx q[2];
rz(-0.6232647) q[2];
sx q[2];
rz(0.40568144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78923046) q[1];
sx q[1];
rz(-1.8327729) q[1];
sx q[1];
rz(1.0220703) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0960363) q[3];
sx q[3];
rz(-0.48431319) q[3];
sx q[3];
rz(0.37176311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(1.6020927) q[2];
rz(-2.0942073) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(1.2654977) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78028107) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(-0.58746946) q[0];
rz(-2.7777708) q[1];
sx q[1];
rz(-1.9963341) q[1];
sx q[1];
rz(2.2344373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5575378) q[0];
sx q[0];
rz(-1.4901925) q[0];
sx q[0];
rz(0.35146468) q[0];
rz(2.467359) q[2];
sx q[2];
rz(-0.57949726) q[2];
sx q[2];
rz(-2.5052239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2550019) q[1];
sx q[1];
rz(-2.0200517) q[1];
sx q[1];
rz(-1.8043065) q[1];
x q[2];
rz(-1.7916405) q[3];
sx q[3];
rz(-2.7691602) q[3];
sx q[3];
rz(0.34162921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63984799) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(2.7196344) q[2];
rz(-0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(-2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.52787732) q[0];
sx q[0];
rz(-1.2953225) q[0];
sx q[0];
rz(-3.0184826) q[0];
rz(0.70115024) q[1];
sx q[1];
rz(-1.9155365) q[1];
sx q[1];
rz(1.6446) q[1];
rz(-1.20303) q[2];
sx q[2];
rz(-0.90819539) q[2];
sx q[2];
rz(2.09203) q[2];
rz(-1.1604068) q[3];
sx q[3];
rz(-1.3153362) q[3];
sx q[3];
rz(1.5462331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
